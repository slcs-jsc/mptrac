#!/bin/bash
#
# Run one (scenario, particle-count) cell of the exchange-mode I/O campaign.
# Runs benchmark.sh for layouts 1x4,2x2,4x1 x modes 0,1 x 3 repeats with the
# locked campaign configuration, then merges the case summary into
# $DD_CAMPAIGN_ROOT/campaign_summary.csv.
#
# Submit from the campaign root directory (job logs land in the submit dir):
#   sbatch --export=ALL,DD_CAMPAIGN_ROOT=$PWD,DD_CAMPAIGN_SCENARIO=A,DD_CAMPAIGN_COUNT=10000000 \
#       /path/to/mptrac/tests/dd_test/campaign_exchange_run.sbatch
#
# Environment (all required unless noted):
#   DD_CAMPAIGN_ROOT     campaign root prepared by campaign_exchange_prep.sh
#   DD_CAMPAIGN_SCENARIO A (zonal transport) or B (perturbations)
#   DD_CAMPAIGN_COUNT    global particle count of the prepared init cell
#   DD_CAMPAIGN_SOURCE   mptrac checkout (default: repo containing this script)
#
# Walltime guidance: JUWELS Booster accepts at most 24 hours per job. The
# default is 12 hours; the 100M-particle cell may need sbatch --time=24:00:00.
#SBATCH --job-name=dd_exch_run
#SBATCH --account=slmet
#SBATCH --partition=booster
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=12
#SBATCH --hint=nomultithread
#SBATCH --time=12:00:00
set -euo pipefail

script_path=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
default_source=$(cd "$script_path/../.." && pwd)

root=${DD_CAMPAIGN_ROOT:?Set DD_CAMPAIGN_ROOT}
scenario=${DD_CAMPAIGN_SCENARIO:?Set DD_CAMPAIGN_SCENARIO (A or B)}
count=${DD_CAMPAIGN_COUNT:?Set DD_CAMPAIGN_COUNT}
mptrac_dir=${DD_CAMPAIGN_SOURCE:-$default_source}

ml nvidia-compilers ParaStationMPI CMake netCDF/4.9.3

wind_dir=$root/wind/$scenario
init_dir=$root/init/c$count
work_dir=$root/runs/$scenario/c$count

for dep in "$wind_dir" "$init_dir"; do
    if [[ ! -d $dep ]]; then
        echo "[ERROR] Missing campaign input: $dep (run campaign_exchange_prep.sh)" >&2
        exit 1
    fi
done
if [[ -e $work_dir ]]; then
    echo "[ERROR] Work directory already exists: $work_dir" >&2
    echo "[ERROR] Remove it or choose a different DD_CAMPAIGN_COUNT/SCENARIO." >&2
    exit 1
fi

case $scenario in
    A)
        export DD_WIND_ALPHA=90 DD_WIND_SPEED=50
        export DD_WIND_PERT_U=0 DD_WIND_PERT_V=0 DD_WIND_PERT_W=0
        ;;
    B)
        export DD_WIND_ALPHA=90 DD_WIND_SPEED=0
        export DD_WIND_PERT_U=3 DD_WIND_PERT_V=3 DD_WIND_PERT_W=0.03 DD_WIND_PERT_K=24
        ;;
    *)
        echo "[ERROR] Unknown DD_CAMPAIGN_SCENARIO: $scenario" >&2
        exit 1
        ;;
esac

export OMP_NUM_THREADS=12
export OMP_PLACES=cores OMP_PROC_BIND=close
export DD_BENCH_CPUS_PER_TASK=12
export DD_COMPILER=nvc DD_MPICC=mpicc
export DD_USE_GPU=1 DD_THRUST=1
export DD_MET_NX=2880 DD_MET_NY=1441 DD_MET_NZ=137 DD_MET_EP=140
export DD_MET_NC_READ_MODE=1
export DD_SIMULATION_HOURS=6 DD_SORT_DT=1800 DD_DIFFUSION=0 DD_ATM_DT_OUT=86400
export DD_BENCH_COUNTS=$count
export DD_BENCH_LAYOUTS=${DD_BENCH_LAYOUTS:-1x4,2x2,4x1}
export DD_BENCH_MODES=${DD_BENCH_MODES:-0,1}
export DD_BENCH_REPEATS=${DD_BENCH_REPEATS:-3}
n_layouts=$(tr , \  <<< "$DD_BENCH_LAYOUTS" | wc -w)
n_modes=$(tr , \  <<< "$DD_BENCH_MODES" | wc -w)
export DD_WIND_DIR=$wind_dir
export DD_BENCH_INIT_SOURCE=$init_dir
export DD_BENCH_MAX_DISK_GIB=$(python3 -c 'import sys
count, nq, layouts, outputs, modes, repeats = map(int, sys.argv[1:7])
gib = count * 8 * (4 + nq) * layouts * (1 + outputs * modes * repeats) / 2**30
print(int(gib * 1.1) + 200)' "$count" 5 "$n_layouts" 7 "$n_modes" "$DD_BENCH_REPEATS")

echo "[INFO] Campaign cell: scenario=$scenario count=$count"
echo "[INFO] Wind: $wind_dir | Init: $init_dir | Work: $work_dir"
echo "[INFO] Disk guard: ${DD_BENCH_MAX_DISK_GIB} GiB"

DD_BENCH_WORK_DIR=$work_dir bash "$mptrac_dir/tests/dd_test/benchmark.sh"

python3 - "$work_dir/summary.csv" "$root/campaign_summary.csv" "$scenario" "$count" <<'PY'
import csv, sys
from pathlib import Path
src, dst, scenario, count = sys.argv[1], Path(sys.argv[2]), sys.argv[3], sys.argv[4]
with open(src, newline="") as stream:
    rows = list(csv.DictReader(stream))
if not rows:
    raise SystemExit(f"[ERROR] Empty summary: {src}")
extra = {"scenario": scenario, "campaign_count": count}
fields = list(extra) + list(rows[0])
write_header = not dst.exists()
with open(dst, "a", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=fields)
    if write_header:
        writer.writeheader()
    for row in rows:
        writer.writerow({**extra, **row})
print(f"[INFO] Appended {len(rows)} rows to {dst}")
PY

echo "[INFO] Cell complete: $work_dir"
