#!/usr/bin/env bash
#
# Prepare inputs for the exchange-mode I/O campaign (see README.md).
#
# Creates, under $DD_CAMPAIGN_ROOT:
#   wind/<scenario>/        wind_YYYY_MM_DD_HH.nc files per wind scenario
#   init/c<count>/          global init.bin plus per-layout rank partitions
#                           (init_layout<L>/data/data.<i>/init.bin)
#   source_revision.txt     repository revision used for the campaign
#
# Usage (run on a cluster login/compute node with build tools, or locally):
#   DD_CAMPAIGN_ROOT=/path/to/campaign ./campaign_exchange_prep.sh
#
# Environment:
#   DD_CAMPAIGN_COUNTS     comma-separated global particle counts
#                          (default: 100000,1000000,10000000,100000000)
#   DD_CAMPAIGN_LAYOUTS    comma-separated layouts (default: 1x4,2x2,4x1)
#   DD_CAMPAIGN_SCENARIOS  comma-separated wind scenarios (default: A,B)
#   DD_MET_NX/NY/NZ        met grid (default: 2880 1441 137)
#   DD_MET_EP              met pressure levels (default: 140)
#   DD_SIMULATION_HOURS    simulation hours (default: 6)
#   DD_PREP_COMPILER       compiler for wind/time2jsec (default: gcc)
#   DD_PREP_MAKE_J         -j level for prep builds (default: 8)
#   DD_INIT_*              see init_global.sh (seed, box, compiler, ...)
#
# The script is idempotent only for wind files (existing non-empty files are
# kept). Initialized particle files are always regenerated.
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mptrac_dir=$(cd "$script_dir/../.." && pwd)
export LD_LIBRARY_PATH="$mptrac_dir/libs/build/lib:${LD_LIBRARY_PATH:-}"

root=${DD_CAMPAIGN_ROOT:?Set DD_CAMPAIGN_ROOT to the campaign output directory}
IFS=, read -r -a counts <<< "${DD_CAMPAIGN_COUNTS:-100000,1000000,10000000,100000000}"
IFS=, read -r -a layouts <<< "${DD_CAMPAIGN_LAYOUTS:-1x4,2x2,4x1}"
IFS=, read -r -a scenarios <<< "${DD_CAMPAIGN_SCENARIOS:-A,B}"
met_nx=${DD_MET_NX:-2880}
met_ny=${DD_MET_NY:-1441}
met_nz=${DD_MET_NZ:-137}
met_ep=${DD_MET_EP:-140}
hours=${DD_SIMULATION_HOURS:-6}
compiler=${DD_PREP_COMPILER:-gcc}
make_j=${DD_PREP_MAKE_J:-8}
dt_met=3600
ntimesteps=$((hours * 3600 / dt_met + 1))

mkdir -p "$root"
root=$(cd "$root" && pwd)
git -C "$mptrac_dir" rev-parse HEAD > "$root/source_revision.txt" 2>/dev/null || \
    printf 'unknown\n' > "$root/source_revision.txt"

echo "[INFO] Campaign root: $root"
echo "[INFO] Counts: ${counts[*]}"
echo "[INFO] Layouts: ${layouts[*]}"
echo "[INFO] Scenarios: ${scenarios[*]}"

#################################################################################################################
# Wind scenarios (generated once, shared by all counts, modes, layouts)
#################################################################################################################

scenario_alpha() { case $1 in A) echo 90;; B) echo 90;; *) echo "[ERROR] Unknown scenario $1" >&2; exit 1;; esac; }
scenario_speed() { case $1 in A) echo 50;; B) echo 0;; *) echo "[ERROR] Unknown scenario $1" >&2; exit 1;; esac; }
scenario_pertu() { case $1 in A) echo 0;; B) echo 3;; *) echo "[ERROR] Unknown scenario $1" >&2; exit 1;; esac; }
scenario_pertv() { case $1 in A) echo 0;; B) echo 3;; *) echo "[ERROR] Unknown scenario $1" >&2; exit 1;; esac; }
scenario_pertw() { case $1 in A) echo 0;; B) echo 0.03;; *) echo "[ERROR] Unknown scenario $1" >&2; exit 1;; esac; }
scenario_k()    { case $1 in A) echo 24;; B) echo 24;; *) echo "[ERROR] Unknown scenario $1" >&2; exit 1;; esac; }

if [[ ${#scenarios[@]} -gt 0 ]]; then
    cd "$mptrac_dir/src"
    make clean > "$root/build_prep.log" 2>&1
    make -j "$make_j" wind time2jsec STATIC=0 GPU=0 THRUST=0 COMPILER="$compiler" \
        DEFINES="-DEX=$((met_nx + 4)) -DEY=$((met_ny + 5)) -DEP=$met_ep" \
        >> "$root/build_prep.log" 2>&1
    t0=$("$mptrac_dir/src/time2jsec" 2022 6 2 0 0 0 0)
    for scenario in "${scenarios[@]}"; do
        wind_dir=$root/wind/$scenario
        mkdir -p "$wind_dir"
        cat > "$wind_dir/wind.ctl" <<EOF
WIND_T0 = $t0
WIND_NX = $met_nx
WIND_NY = $met_ny
WIND_NZ = $met_nz
GRID_NX = $met_nx
GRID_NY = $met_ny
WIND_Z0 = 0
WIND_Z1 = 60
WIND_U0 = $(scenario_speed "$scenario")
WIND_U1 = $(scenario_speed "$scenario")
WIND_W0 = 0
WIND_PERT_U = $(scenario_pertu "$scenario")
WIND_PERT_V = $(scenario_pertv "$scenario")
WIND_PERT_W = $(scenario_pertw "$scenario")
WIND_PERT_K = $(scenario_k "$scenario")
WIND_ALPHA = $(scenario_alpha "$scenario")
WIND_LAT_REVERSE = 1
EOF
        echo "[INFO] Generating wind scenario $scenario in $wind_dir"
        for ((i = 0; i < ntimesteps; i++)); do
            new_t0=$(echo "scale=6; $t0 + $i * $dt_met" | bc)
            sed -i "s/^WIND_T0 = .*/WIND_T0 = $new_t0/" "$wind_dir/wind.ctl"
            date -u -d "1970-01-01 + $new_t0 seconds" +"wind_%Y_%m_%d_%H.nc" \
                > "$wind_dir/.expected_name"
            expected=$root/wind/$scenario/$(cat "$wind_dir/.expected_name")
            if [[ -s $expected ]]; then
                echo "[INFO] Keeping existing $expected"
                continue
            fi
            (cd "$wind_dir" && "$mptrac_dir/src/wind" wind.ctl wind) \
                >> "$root/wind_generation.log" 2>&1
        done
        rm -f "$wind_dir/.expected_name"
        n=$(find "$wind_dir" -name 'wind_*.nc' | wc -l)
        if [[ $n -ne $ntimesteps ]]; then
            echo "[ERROR] Scenario $scenario has $n wind files, expected $ntimesteps" >&2
            exit 1
        fi
    done
fi

#################################################################################################################
# Global initial conditions, partitioned per layout
#################################################################################################################

lcm() { python3 -c 'import sys, math; v=1
for x in sys.argv[1:]: v = v * int(x) // math.gcd(v, int(x))
print(v)' "$@"; }

step_lon=1
step_lat=1
for layout in "${layouts[@]}"; do
    dlon=${layout%x*}
    dlat=${layout#*x}
    step_lon=$(lcm "$step_lon" "$dlon")
    step_lat=$(lcm "$step_lat" "$dlat")
done

for count in "${counts[@]}"; do
    read -r global_lon global_lat < <(
        python3 "$script_dir/benchmark_grid.py" "$count" "$step_lon" "$step_lat")
    init_dir=$root/init/c$count
    mkdir -p "$init_dir"
    echo "[INFO] Generating global init for count=$count grid=${global_lon}x${global_lat}"
    DD_INIT_COMPILER="${DD_INIT_COMPILER:-$compiler}" \
        "$script_dir/init_global.sh" "$init_dir/init_global.bin" \
        "$global_lon" "$global_lat"
    python3 "$script_dir/init_partition.py" "$init_dir/init_global.bin" \
        "$init_dir" "$met_nx" "$met_ny" "$(IFS=,; echo "${layouts[*]}")"
    python3 - "$init_dir" <<'PY'
import struct, sys
from pathlib import Path
root = Path(sys.argv[1])
with open(root / "init_global.bin", "rb") as stream:
    _, count = struct.unpack("=ii", stream.read(8))
for layout_dir in sorted(root.glob("init_layout*")):
    total = 0
    for init in sorted(layout_dir.glob("data/data.*/init.bin")):
        with open(init, "rb") as stream:
            _, n = struct.unpack("=ii", stream.read(8))
        total += n
    if total != count:
        raise SystemExit(f"[ERROR] {layout_dir}: partition total {total} != {count}")
    print(f"[INFO] {layout_dir.name}: {total} particles across "
          f"{len(list(layout_dir.glob('data/data.*/init.bin')))} ranks")
PY
done

echo "[INFO] Campaign inputs ready under $root"
