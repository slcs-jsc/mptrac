#!/usr/bin/env bash
set -euo pipefail

# Compare complete local met_t fields from both DD NetCDF readers. This is
# deliberately a small CPU test; the ordinary DD test checks trajectories.
script_dir=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$script_dir/../.." && pwd)
work=${DD_MET_COMPARE_WORK_DIR:-$(mktemp -d)}
mpi_cc=${DD_MPICC:-mpicc}
if [[ -z ${DD_MPICC:-} && -x /usr/bin/mpicc ]]; then
    mpi_cc=/usr/bin/mpicc
fi
mkdir -p "$work"
echo "[INFO] Meteorological reader test directory: $work"

(
    cd "$work"
    DD_DOMAINS_LON=2 DD_DOMAINS_LAT=2 \
    DD_PARTICLES_LON_STRIDE=2 DD_PARTICLES_LAT_STRIDE=2 \
    DD_MET_NX=20 DD_MET_NY=11 DD_MET_NZ=8 \
    DD_SIMULATION_HOURS=1 DD_THRUST=0 DD_USE_GPU=0 \
    DD_MPICC="$mpi_cc" \
    bash "$script_dir/run.sh" compile skip_run skip_compare > setup.log 2>&1
)

defs=(-DMPI -DDD -DNP=8 -DEX_GLOB=24 -DEY_GLOB=16
      -DEP=90 -DEX=20 -DEY=15 -DNQ=5)
"$mpi_cc" -O2 -Werror -Wall -Wextra -fopenmp \
    -I "$repo/src" -I "$repo/libs/build/include" "${defs[@]}" \
    "$script_dir/met_reader_compare.c" "$repo/src/mptrac.o" \
    -L "$repo/libs/build/lib" -lnetcdf -lgsl -lgslcblas -lm \
    -o "$work/met_reader_compare"
"${DD_CC:-cc}" -O2 -Werror -Wall -Wextra \
    "$script_dir/met_reader_encode.c" -lnetcdf -lm \
    -o "$work/met_reader_encode"

export LD_LIBRARY_PATH="$repo/libs/build/lib:${LD_LIBRARY_PATH:-}"
wind=$work/data/wind_data/wind_2022_06_02_00.nc
mkdir -p "$work/encoded"
cp "$wind" "$work/encoded/wind_2022_06_02_00.nc"
"$work/met_reader_encode" "$work/encoded/wind_2022_06_02_00.nc"

for layout in 1x1 1x4 2x2 4x1; do
    nx=${layout%x*}
    ny=${layout#*x}
    "${DD_MPIRUN:-mpirun}" --oversubscribe -np "$((nx * ny))" \
        "$work/met_reader_compare" "$work/data/config.ctl" \
        "$work/encoded/wind_2022_06_02_00.nc" "$nx" "$ny" \
        > "$work/compare_${layout}.log" 2>&1
    count=$(grep -c "different_bytes=0" "$work/compare_${layout}.log")
    if (( count != nx * ny )); then
        echo "[ERROR] Incomplete $layout comparison" >&2
        exit 1
    fi
    echo "[PASS] $layout: $count complete bytewise met_t comparisons"
done

# The existing reader supports a duplicated 0/360 degree endpoint when the
# DD layout has no zonal split. Compare that mapping separately.
mkdir -p "$work/endpoint"
cp "$wind" "$work/endpoint/wind_2022_06_02_00.nc"
"$work/met_reader_encode" "$work/endpoint/wind_2022_06_02_00.nc" \
    --duplicate-endpoint
for layout in 1x1 1x4; do
    nx=${layout%x*}
    ny=${layout#*x}
    "${DD_MPIRUN:-mpirun}" --oversubscribe -np "$((nx * ny))" \
        "$work/met_reader_compare" "$work/data/config.ctl" \
        "$work/endpoint/wind_2022_06_02_00.nc" "$nx" "$ny" \
        > "$work/compare_endpoint_${layout}.log" 2>&1
    count=$(grep -c "different_bytes=0" \
        "$work/compare_endpoint_${layout}.log")
    if (( count != nx * ny )); then
        echo "[ERROR] Incomplete duplicate-endpoint $layout comparison" >&2
        exit 1
    fi
    echo "[PASS] duplicate endpoint $layout: $count comparisons"
done
