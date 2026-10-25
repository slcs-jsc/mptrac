#!/usr/bin/env bash
set -euo pipefail

# GPU-DD NetCDF scaling study. Dry-run is the default; DD_IO_RUN=1 executes.
script_dir=$(cd "$(dirname "$0")" && pwd)
mptrac_dir=$(cd "$script_dir/../.." && pwd)
work_root=${DD_IO_WORK_DIR:-$PWD/dd_io_scaling_137}
resolutions=${DD_IO_RESOLUTIONS:-1,0.5,0.25,0.125}
export DD_BENCH_LAYOUTS=${DD_BENCH_LAYOUTS:-1x4,2x2,4x1}
export DD_BENCH_COUNTS=${DD_BENCH_COUNTS:-10000000}
export DD_BENCH_MODES=${DD_BENCH_MODES:-1}
export DD_BENCH_REPEATS=${DD_BENCH_REPEATS:-3}
export DD_BENCH_RANKS=${DD_BENCH_RANKS:-4}
export DD_BENCH_NODES=${DD_BENCH_NODES:-1}
export DD_BENCH_RANKS_PER_NODE=${DD_BENCH_RANKS_PER_NODE:-4}
export DD_BENCH_CANONICAL_INIT=1
export DD_SIMULATION_HOURS=${DD_SIMULATION_HOURS:-6}
export DD_SORT_DT=${DD_SORT_DT:-1800}
export DD_DIFFUSION=${DD_DIFFUSION:-0}
export DD_USE_GPU=${DD_USE_GPU:-1}
export DD_IO_GPU_SAMPLE=${DD_IO_GPU_SAMPLE:-1}
export DD_IO_CHECK_FINITE=1
export DD_IO_PROFILE=${DD_IO_PROFILE:-1}
export DD_BENCH_MAX_DISK_GIB=${DD_BENCH_MAX_DISK_GIB:-32}
export DD_MET_NZ=137 DD_MET_EP=140

case ",${resolutions}," in
    *,,*) echo "[ERROR] Empty resolution in DD_IO_RESOLUTIONS" >&2; exit 1 ;;
esac
mkdir -p "$work_root"
work_root=$(cd "$work_root" && pwd)
report_results() {
    if [[ -f $work_root/plan.csv ]]; then
        python3 "$script_dir/io_scaling_report.py" "$work_root" \
            --files-per-run "$((DD_SIMULATION_HOURS + 1))" \
            --repeats "$DD_BENCH_REPEATS" --ranks "$DD_BENCH_RANKS" \
            --gpus "$((DD_USE_GPU * DD_BENCH_RANKS))" || true
    fi
}
trap report_results EXIT
if [[ ${DD_IO_RUN:-0} == 1 && $DD_IO_PROFILE != 1 && -z ${DD_IO_TRACE_LIB:-} ]]; then
    export DD_IO_TRACE_LIB=$work_root/libnetcdf_read_trace.so
    "${DD_MPICC:-mpicc}" -O2 -fPIC -shared \
        "$script_dir/netcdf_read_trace.c" -o "$DD_IO_TRACE_LIB" \
        -I "$mptrac_dir/libs/build/include" \
        -L "$mptrac_dir/libs/build/lib" -lnetcdf -ldl
fi
if [[ ${DD_IO_RUN:-0} == 1 && -z ${DD_MET_VALIDATE_EXEC:-} ]]; then
    export DD_MET_VALIDATE_EXEC=$work_root/validate_wind_grid
    "${DD_MPICC:-mpicc}" -O2 "$script_dir/validate_wind_grid.c" \
        -o "$DD_MET_VALIDATE_EXEC" \
        -I "$mptrac_dir/libs/build/include" \
        -L "$mptrac_dir/libs/build/lib" -lnetcdf -lm
fi
plan=$work_root/plan.csv
printf '%s\n' 'resolution,nx,ny,nz,ex_global,ey_global,ep,wind_met_gib,wind_buffer_gib,largest_local_ex,largest_local_ey,two_local_met_gib,estimated_gpu_gib,wind_files_gib' > "$plan"

IFS=, read -r -a resolution_list <<< "$resolutions"
for resolution in "${resolution_list[@]}"; do
    case $resolution in
        1) DD_MET_NX=360; DD_MET_NY=181 ;;
        0.5) DD_MET_NX=720; DD_MET_NY=361 ;;
        0.25) DD_MET_NX=1440; DD_MET_NY=721 ;;
        0.125) DD_MET_NX=2880; DD_MET_NY=1441 ;;
        *) echo "[ERROR] Unsupported resolution: $resolution" >&2; exit 1 ;;
    esac
    export DD_MET_NX DD_MET_NY
    export DD_BENCH_WORK_DIR=$work_root/res_${resolution//./p}
    mkdir -p "$DD_BENCH_WORK_DIR"
    python3 - "$resolution" "$DD_MET_NX" "$DD_MET_NY" "$plan" \
        "${DD_IO_GPU_MEMORY_GIB:-40}" "${DD_IO_HOST_MEMORY_GIB:-503}" \
        "${DD_IO_AVAILABLE_DISK_GIB:-400}" "$work_root" \
        "$((DD_SIMULATION_HOURS + 1))" <<'PY'
import csv
import math
import shutil
import sys

resolution, nx, ny, path, gpu_limit, host_limit, disk_limit, root, nfiles = sys.argv[1:]
nx, ny = int(nx), int(ny)
nfiles = int(nfiles)
gpu_limit, host_limit, disk_limit = map(float, (gpu_limit, host_limit, disk_limit))
ep = 140
gib = 2**30
met_bytes = lambda ex, ey: 76 * ex * ey * ep + 96 * ex * ey + 32768
global_bytes = met_bytes(nx + 4, ny + 5)
local_sizes = []
for nlon, nlat in ((1, 4), (2, 2), (4, 1)):
    ex = nx // nlon + 10
    ey = ny // nlat + 10
    local_sizes.append((met_bytes(ex, ey), ex, ey))
local_bytes, ex, ey = max(local_sizes)
wind_buffer = 4 * (nx + 4) * (ny + 5) * ep
# Particle/cache/DD storage and CUDA/OpenACC runtime receive 4 GiB reserve.
estimated_gpu = (2 * local_bytes + 4 * gib) / gib
file_bytes = 4 * nx * ny * (11 * 137 + 24) + 8 * (nx + ny + 137)
wind_files = nfiles * file_bytes / gib
particle_disk = 16  # conservative per-resolution binary output allowance
if estimated_gpu > 0.9 * gpu_limit:
    raise SystemExit(f"[ERROR] GPU preflight failed: {estimated_gpu:.1f} GiB > 90% of {gpu_limit} GiB")
if (global_bytes + wind_buffer > 0.8 * host_limit * gib):
    raise SystemExit("[ERROR] Host preflight failed for global wind generator")
if wind_files + particle_disk > disk_limit:
    raise SystemExit("[ERROR] Disk preflight failed for this resolution")
if wind_files + particle_disk > 0.8 * shutil.disk_usage(root).free / gib:
    raise SystemExit("[ERROR] Insufficient free file-system space")
with open(path, "a", newline="") as stream:
    csv.writer(stream).writerow((resolution, nx, ny, 137, nx + 4, ny + 5, ep,
        f"{global_bytes/gib:.3f}", f"{wind_buffer/gib:.3f}", ex, ey,
        f"{2*local_bytes/gib:.3f}", f"{estimated_gpu:.3f}", f"{wind_files:.3f}"))
print(f"[INFO] {resolution}°: wind host {global_bytes/gib:.2f} GiB + "
      f"{wind_buffer/gib:.2f} GiB buffer; max GPU estimate "
      f"{estimated_gpu:.2f} GiB; {nfiles} wind files ~{wind_files:.2f} GiB")
PY
    if [[ ${DD_IO_RUN:-0} == 1 ]]; then
        if [[ -z ${SLURM_JOB_ID:-} ]]; then
            echo "[ERROR] Execute the GPU matrix inside a Slurm allocation" >&2
            exit 1
        fi
        if [[ -e $DD_BENCH_WORK_DIR/summary.csv ]]; then
            echo "[ERROR] Existing results: $DD_BENCH_WORK_DIR" >&2
            exit 1
        fi
        "$script_dir/benchmark.sh" > "$DD_BENCH_WORK_DIR/benchmark.log" 2>&1
        python3 - "$DD_BENCH_WORK_DIR" "$resolution" \
            "$((DD_SIMULATION_HOURS + 1))" <<'PY'
import csv
import pathlib
import sys
root, resolution, expected = pathlib.Path(sys.argv[1]), sys.argv[2], int(sys.argv[3])
files = sorted((root / "wind_data").glob("wind_*.nc"))
if len(files) != expected:
    raise SystemExit(f"[ERROR] Expected {expected} wind files, found {len(files)}")
with (root / "wind_files.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(("resolution", "filename", "bytes"))
    writer.writerows((resolution, path.name, path.stat().st_size) for path in files)
print(f"[INFO] {resolution}° actual wind disk: "
      f"{sum(p.stat().st_size for p in files)/2**30:.3f} GiB")
PY
    else
        DD_BENCH_DRY_RUN=1 "$script_dir/benchmark.sh" \
            > "$DD_BENCH_WORK_DIR/dry_run.log" 2>&1
    fi
done
echo "[INFO] Memory/disk plan: $plan"
