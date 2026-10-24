#!/bin/bash
set -euo pipefail

# Run inside a JUWELS Booster allocation, or set DD_USE_GPU=0 for a local
# CPU smoke test. All output lives below DD_BENCH_WORK_DIR.
script_dir=$(cd "$(dirname "$0")" && pwd)
mptrac_dir=$(cd "$script_dir/../.." && pwd)
export LD_LIBRARY_PATH="$mptrac_dir/libs/build/lib:${LD_LIBRARY_PATH:-}"
work_root=${DD_BENCH_WORK_DIR:-$PWD/dd_benchmark}
mkdir -p "$work_root"
work_root=$(cd "$work_root" && pwd)

export DD_BENCHMARK=1
export DD_USE_GPU=${DD_USE_GPU:-1}
export DD_DOMAINS_LON=${DD_DOMAINS_LON:-2}
export DD_DOMAINS_LAT=${DD_DOMAINS_LAT:-2}
export DD_SIMULATION_HOURS=${DD_SIMULATION_HOURS:-6}
export DD_SORT_DT=${DD_SORT_DT:-1800}
export DD_DIFFUSION=${DD_DIFFUSION:-0}
export DD_ATM_DT_OUT=${DD_ATM_DT_OUT:-86400}
export DD_NP_MARGIN_PERCENT=${DD_NP_MARGIN_PERCENT:-25}
export DD_THRUST=${DD_THRUST:-$([[ $DD_USE_GPU == 1 ]] && echo 1 || echo 0)}
export DD_COMPILER=${DD_COMPILER:-$([[ $DD_USE_GPU == 1 ]] && echo nvc || echo gcc)}
export DD_MPICC=${DD_MPICC:-mpicc}
export DD_MPIRUN=${DD_MPIRUN:-mpirun}
export DD_WIND_ALPHA=${DD_WIND_ALPHA:-90}
export DD_WIND_PERT_U=${DD_WIND_PERT_U:-0}
export DD_WIND_PERT_V=${DD_WIND_PERT_V:-0}
export DD_WIND_PERT_W=${DD_WIND_PERT_W:-0}
export DD_WIND_PERT_K=${DD_WIND_PERT_K:-24}

counts=${DD_BENCH_COUNTS:-100000,1000000,10000000}
modes=${DD_BENCH_MODES:-${DD_EXCHANGE_MODE:-0,1}}
repeats=${DD_BENCH_REPEATS:-1}
disk_limit=${DD_BENCH_MAX_DISK_GIB:-64}
memory_limit=${DD_BENCH_MAX_PARTICLE_MEMORY_GIB:-24}
counts=${counts//,/ }
modes=${modes//,/ }
read -r -a count_list <<< "$counts"
read -r -a mode_list <<< "$modes"
if (( ${#count_list[@]} == 0 || ${#mode_list[@]} == 0 )); then
    echo "[ERROR] Empty benchmark particle-count or exchange-mode matrix" >&2
    exit 1
fi
for value in "$DD_DOMAINS_LON" "$DD_DOMAINS_LAT" "$DD_SIMULATION_HOURS" \
    "$DD_ATM_DT_OUT" "$DD_NP_MARGIN_PERCENT" "$repeats" "$disk_limit" \
    "$memory_limit" \
    "${count_list[@]}"; do
    if [[ ! $value =~ ^[0-9]+$ ]]; then
        echo "[ERROR] Benchmark sizes, hours, repeats, and disk limit must be integers" >&2
        exit 1
    fi
done
for mode in "${mode_list[@]}"; do
    if [[ $mode != 0 && $mode != 1 ]]; then
        echo "[ERROR] DD_BENCH_MODES accepts only 0 and 1" >&2
        exit 1
    fi
done
if [[ $DD_USE_GPU != 0 && $DD_USE_GPU != 1 ]]; then
    echo "[ERROR] DD_USE_GPU must be 0 or 1" >&2
    exit 1
fi
if [[ -n ${DD_NP_MAX:-} && ! ${DD_NP_MAX} =~ ^[0-9]+$ ]]; then
    echo "[ERROR] DD_NP_MAX must be an integer" >&2
    exit 1
fi
ranks=$((DD_DOMAINS_LON * DD_DOMAINS_LAT))
if (( ranks < 2 || repeats < 1 || DD_SIMULATION_HOURS < 1 ||
      DD_ATM_DT_OUT < 600 || DD_ATM_DT_OUT % 600 != 0 )); then
    echo "[ERROR] Use at least two ranks, positive repeats/hours, and an output interval divisible by 600 s" >&2
    exit 1
fi
t0=707443200  # 2022-06-02 00:00 UTC, as in run.sh.
if (( t0 % DD_ATM_DT_OUT != 0 )); then
    echo "[ERROR] DD_ATM_DT_OUT must divide the initial time to retain the first snapshot" >&2
    exit 1
fi
tstop=$((t0 + DD_SIMULATION_HOURS * 3600))
noutputs=$((tstop / DD_ATM_DT_OUT - (t0 - 1) / DD_ATM_DT_OUT))
if (( tstop % DD_ATM_DT_OUT != 0 )); then noutputs=$((noutputs + 1)); fi
nq=5
if [[ ${DD_ISOSURF:-0} == 3 ]]; then nq=6; fi

# Derive nearly square particle grids for each requested approximate count.
# Explicit lon/lat strides override this calculation for a single-count run.
if [[ -n ${DD_PARTICLES_LON_STRIDE:-} || -n ${DD_PARTICLES_LAT_STRIDE:-} ]]; then
    if [[ -z ${DD_PARTICLES_LON_STRIDE:-} || -z ${DD_PARTICLES_LAT_STRIDE:-} ||
          ${#count_list[@]} -ne 1 ]]; then
        echo "[ERROR] Set both DD_PARTICLES_*_STRIDE values with one DD_BENCH_COUNTS entry" >&2
        exit 1
    fi
fi
declare -a lon_list lat_list total_list capacity_list
estimated_bytes=0
for target in "${count_list[@]}"; do
    if (( target < ranks )); then
        echo "[ERROR] Target particle count must be at least the number of MPI ranks" >&2
        exit 1
    fi
    if [[ -n ${DD_PARTICLES_LON_STRIDE:-} ]]; then
        lon=$DD_PARTICLES_LON_STRIDE
        lat=$DD_PARTICLES_LAT_STRIDE
    else
        read -r lon lat < <(python3 -c 'import math,sys
n=int(sys.argv[1]); ranks=int(sys.argv[2]); per=math.ceil(n/ranks)
lon=max(1,round(math.sqrt(2*per))); lat=math.ceil(per/lon)
print(lon,lat)' "$target" "$ranks")
    fi
    if [[ ! $lon =~ ^[0-9]+$ || ! $lat =~ ^[0-9]+$ ]] ||
       (( lon < 1 || lat < 1 )); then
        echo "[ERROR] Particle strides must be positive integers" >&2
        exit 1
    fi
    per_rank=$((lon * lat))
    total=$((per_rank * ranks))
    capacity=${DD_NP_MAX:-$((per_rank + (per_rank * DD_NP_MARGIN_PERCENT + 99) / 100 + 1))}
    if (( capacity <= per_rank || capacity > 2147483647 || total > 2147483647 )); then
        echo "[ERROR] NP capacity or particle IDs exceed the supported signed-int range" >&2
        exit 1
    fi
    if (( capacity * (8 * (4 + nq) + 80) > memory_limit * 1073741824 )); then
        echo "[ERROR] Estimated per-rank particle arrays exceed DD_BENCH_MAX_PARTICLE_MEMORY_GIB=$memory_limit" >&2
        exit 1
    fi
    echo "[INFO] Target $target: actual $total, initial/rank $per_rank, NP/rank $capacity, particle arrays/rank ~$((capacity * (8 * (4 + nq) + 80) / 1048576)) MiB"
    lon_list+=("$lon")
    lat_list+=("$lat")
    total_list+=("$total")
    capacity_list+=("$capacity")
    # One hard-linked initial state per size; all cases have noutputs files.
    estimated_bytes=$((estimated_bytes + total * 8 * (4 + nq) *
        (1 + noutputs * ${#mode_list[@]} * repeats)))
done
if (( estimated_bytes > disk_limit * 1073741824 )); then
    echo "[ERROR] Estimated particle files exceed DD_BENCH_MAX_DISK_GIB=$disk_limit" >&2
    exit 1
fi
echo "[INFO] Estimated particle-file storage: $((estimated_bytes / 1073741824)) GiB"
echo "[INFO] Particle snapshots per case: $noutputs"
if [[ ${DD_BENCH_DRY_RUN:-0} == 1 ]]; then
    echo "[INFO] Dry run complete; no particle or wind files were generated"
    exit 0
fi

wind_external=0
if [[ -n ${DD_WIND_DIR:-} ]]; then wind_external=1; fi
wind_dir=${DD_WIND_DIR:-$work_root/wind_data}
mkdir -p "$wind_dir"
wind_dir=$(cd "$wind_dir" && pwd)
export DD_WIND_DIR=$wind_dir
wind_key="hours=$DD_SIMULATION_HOURS alpha=$DD_WIND_ALPHA u=$DD_WIND_PERT_U v=$DD_WIND_PERT_V w=$DD_WIND_PERT_W k=$DD_WIND_PERT_K"
generate_wind=0
if (( ! wind_external )); then
    if [[ -f $wind_dir/.benchmark_config ]]; then
        if [[ $(cat "$wind_dir/.benchmark_config") != "$wind_key" ]]; then
            echo "[ERROR] Cached wind settings differ; use a new DD_BENCH_WORK_DIR" >&2
            exit 1
        fi
    else
        generate_wind=1
    fi
fi
for ((hour=0; hour<=DD_SIMULATION_HOURS; hour++)); do
    filename=$(date -u -d "2022-06-02 00:00:00 UTC + $hour hours" +wind_%Y_%m_%d_%H.nc)
    if [[ ! -f $wind_dir/$filename ]]; then
        if (( ! wind_external )); then
            generate_wind=1
        else
            echo "[ERROR] Missing external meteorological file: $wind_dir/$filename" >&2
            exit 1
        fi
    fi
done
if (( generate_wind )); then
    if [[ -n ${WIND_EXEC:-} ]]; then
        if [[ ! -x $WIND_EXEC ]]; then
            echo "[ERROR] WIND_EXEC is not executable: $WIND_EXEC" >&2
            exit 1
        fi
    else
        echo "[INFO] Building the global-grid wind generator"
        first_capacity=${capacity_list[0]}
        defs="-DNP=$first_capacity -DEX=364 -DEY=186 -DEP=90 -DNQ=$nq"
        build_opt=()
        if [[ $DD_COMPILER == nvc ]]; then
            build_opt=(OPT="${DD_BENCH_OPT:--O3 -Mlarge_arrays}")
        elif [[ -n ${DD_BENCH_OPT:-} ]]; then
            build_opt=(OPT="$DD_BENCH_OPT")
        fi
        make -C "$mptrac_dir/src" clean > "$work_root/build_wind.log" 2>&1
        make -C "$mptrac_dir/src" -j2 wind STATIC=0 GPU=0 \
            COMPILER="$DD_COMPILER" DEFINES="$defs" "${build_opt[@]}" \
            >> "$work_root/build_wind.log" 2>&1
        cp "$mptrac_dir/src/wind" "$work_root/wind_global"
        export WIND_EXEC=$work_root/wind_global
    fi
fi

summary=$work_root/summary.csv
compiled_capacity=-1
for index in "${!count_list[@]}"; do
    lon=${lon_list[$index]}
    lat=${lat_list[$index]}
    total=${total_list[$index]}
    capacity=${capacity_list[$index]}
    per_rank=$((lon * lat))
    init_source=""
    for ((repeat=1; repeat<=repeats; repeat++)); do
        for mode in "${mode_list[@]}"; do
            case_dir=$work_root/n${total}/repeat${repeat}/mode${mode}
            if [[ -e $case_dir/data ]]; then
                echo "[ERROR] Existing benchmark case: $case_dir; use a fresh work directory" >&2
                exit 1
            fi
            mkdir -p "$case_dir"
            export DD_PARTICLES_LON_STRIDE=$lon DD_PARTICLES_LAT_STRIDE=$lat
            export DD_NP_MAX=$capacity DD_EXCHANGE_MODE=$mode
            setup_args=(skip_run)
            if (( compiled_capacity != capacity )); then
                setup_args=(compile "${setup_args[@]}")
                compiled_capacity=$capacity
            fi
            if (( generate_wind )); then
                generate_wind=0
            else
                setup_args+=(skip_wind)
            fi
            if [[ -n $init_source ]]; then export DD_INIT_SOURCE=$init_source
            else unset DD_INIT_SOURCE || true
            fi
            echo "[INFO] Setup particles=$total, mode=$mode, repeat=$repeat, NP=$capacity"
            (cd "$case_dir" && "$script_dir/run.sh" "${setup_args[@]}" \
                > setup.log 2>&1)
            python3 - "$case_dir/data" "$ranks" "$per_rank" <<'PY'
import pathlib
import struct
import sys
root, ranks, expected = pathlib.Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
for rank in range(ranks):
    path = root / f"data.{rank}" / "init.bin"
    with path.open("rb") as stream:
        version, count = struct.unpack("=ii", stream.read(8))
    if version != 100 or count != expected:
        raise SystemExit(f"[ERROR] Unexpected initialization count in {path}: {count}")
PY
            if [[ ! -f $wind_dir/.benchmark_config && $wind_external == 0 ]]; then
                printf '%s\n' "$wind_key" > "$wind_dir/.benchmark_config"
            fi
            if [[ -z $init_source ]]; then init_source=$case_dir/data; fi

            echo "[INFO] Run particles=$total, mode=$mode, repeat=$repeat"
            start=$(date +%s.%N)
            if [[ -n ${SLURM_JOB_ID:-} ]]; then
                if [[ $DD_USE_GPU == 1 ]]; then
                    (cd "$case_dir/data" && srun --ntasks="$ranks" \
                        --gpus-per-task=1 --gpu-bind=single:1 --cpu-bind=cores \
                        "$mptrac_dir/src/trac" dirlist.tab config.ctl init.bin ATM_BASENAME atm) \
                        > "$case_dir/trac.log" 2>&1
                else
                    (cd "$case_dir/data" && srun --ntasks="$ranks" \
                        "$mptrac_dir/src/trac" dirlist.tab config.ctl init.bin ATM_BASENAME atm) \
                        > "$case_dir/trac.log" 2>&1
                fi
            else
                (cd "$case_dir/data" && "$DD_MPIRUN" --oversubscribe -np "$ranks" \
                    "$mptrac_dir/src/trac" dirlist.tab config.ctl init.bin ATM_BASENAME atm) \
                    > "$case_dir/trac.log" 2>&1
            fi
            end=$(date +%s.%N)
            python3 -c 'import sys; print(float(sys.argv[2])-float(sys.argv[1]))' \
                "$start" "$end" > "$case_dir/launcher_wall_seconds.txt"
            python3 "$script_dir/benchmark_report.py" --case "$case_dir" \
                --log "$case_dir/trac.log" --ranks "$ranks" \
                --gpus "$((DD_USE_GPU * ranks))" \
                --lon-domains "$DD_DOMAINS_LON" --lat-domains "$DD_DOMAINS_LAT" \
                --per-rank "$per_rank" --id-stride "${DD_INIT_IDX_STRIDE:-$per_rank}" \
                --capacity "$capacity" --nq "$nq" --mode "$mode" \
                --sort-dt "$DD_SORT_DT" --hours "$DD_SIMULATION_HOURS" \
                --repeat "$repeat" --header > "$case_dir/result.csv" \
                2> "$case_dir/report.log"
            if [[ ! -f $summary ]]; then
                cp "$case_dir/result.csv" "$summary"
            else
                tail -n 1 "$case_dir/result.csv" >> "$summary"
            fi
            cat "$case_dir/report.log"
        done
    done
done
echo "[INFO] Summary: $summary"
