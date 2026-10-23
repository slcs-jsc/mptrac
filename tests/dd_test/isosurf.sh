#!/bin/bash
set -euo pipefail

# Run from any directory. Keep the temporary files when a check fails.
script_dir=$(cd "$(dirname "$0")" && pwd)
mptrac_dir=$(cd "$script_dir/../.." && pwd)
work_dir=${DD_ISOSURF_WORK_DIR:-$(mktemp -d /tmp/mptrac_dd_isosurf.XXXXXX)}
wind_dir=${DD_WIND_DIR:-$work_dir/wind_data}
mkdir -p "$work_dir"/{mode0,mode1,serial}
echo "[INFO] Isosurface test directory: $work_dir"

# The wind generator needs the global grid dimensions; the DD build uses
# smaller local arrays and therefore cannot generate these fields itself.
wind_arg=()
if [[ -f "$wind_dir/wind_2022_06_02_06.nc" ]]; then
    wind_arg=(skip_wind)
else
    make -C "$mptrac_dir/src" -B wind DD=0 MPI=0 GPU=0 STATIC=0 \
        DEFINES='-DEX=364 -DEY=186 -DEP=90 -DNQ=6' \
        > "$work_dir/wind-build.log" 2>&1
    cp "$mptrac_dir/src/wind" "$work_dir/wind"
fi

for mode in 0 1; do
    cd "$work_dir/mode$mode"
    run_wind_arg=("${wind_arg[@]}")
    if [[ $mode -eq 1 ]]; then run_wind_arg=(skip_wind); fi
    DD_MPICC=${DD_MPICC:-/usr/bin/mpicc} \
    DD_MPIRUN=${DD_MPIRUN:-/usr/bin/mpirun} \
    WIND_EXEC="$work_dir/wind" \
    DD_WIND_DIR="$wind_dir" \
    DD_DOMAINS_LON=2 DD_DOMAINS_LAT=1 \
    DD_PARTICLES_LON_STRIDE=12 DD_PARTICLES_LAT_STRIDE=1 \
    DD_INIT_Z0=8 DD_INIT_Z1=12 DD_INIT_DZ=4 \
    DD_INIT_IDX_STRIDE=24 DD_WIND_ALPHA=0 \
    DD_NP_MAX=100 DD_THRUST=0 DD_ISOSURF=3 DD_EXCHANGE_MODE=$mode \
        "$script_dir/run.sh" compile "${run_wind_arg[@]}" skip_compare \
        > "$work_dir/mode$mode.log" 2>&1
done

# Use the same initial particles and full meteorological grid on one rank.
cd "$mptrac_dir/src"
make -B trac DD=0 MPI=0 THRUST=0 STATIC=0 GPU=0 \
    DEFINES='-DNP=100 -DEX=364 -DEY=186 -DEP=90 -DNQ=6' \
    > "$work_dir/serial-build.log" 2>&1
for rank in 0 1; do
    serial_dir="$work_dir/serial/data.$rank"
    mkdir -p "$serial_dir"
    cp "$work_dir/mode0/data/data.$rank/init.nc" "$serial_dir/init.nc"
    printf '.\n' > "$serial_dir/dirlist.tab"
    sed -e 's/^DD_SUBDOMAINS_MERIDIONAL = .*/DD_SUBDOMAINS_MERIDIONAL = 1/' \
        -e 's/^DD_SUBDOMAINS_ZONAL = .*/DD_SUBDOMAINS_ZONAL = 1/' \
        -e 's/^GRID_NX = .*/GRID_NX = 360/' \
        -e 's/^GRID_NY = .*/GRID_NY = 181/' \
        "$work_dir/mode0/data/data.$rank/config.ctl" > "$serial_dir/config.ctl"
    (cd "$serial_dir" && "$mptrac_dir/src/trac" dirlist.tab config.ctl init.nc \
        ATM_BASENAME atm > "$work_dir/serial-rank$rank.log" 2>&1)
done

sed -e 's/^NQ = 6$/NQ = 5/' -e '/^QNT_NAME\[5\] = iso_var$/d' \
    "$work_dir/serial/data.0/config.ctl" > "$work_dir/serial/data.0/missing_iso.ctl"
if (cd "$work_dir/serial/data.0" && "$mptrac_dir/src/trac" dirlist.tab missing_iso.ctl \
    init.nc > "$work_dir/missing-iso.log" 2>&1); then
    echo "[ERROR] Missing iso_var quantity was accepted" >&2
    exit 1
fi
grep -q 'Add quantity iso_var for ISOSURF=1, 2, or 3' \
    "$work_dir/missing-iso.log"

python3 "$script_dir/check_isosurf.py" \
    "$work_dir/mode0/data" "$work_dir/mode1/data" "$work_dir/serial"
