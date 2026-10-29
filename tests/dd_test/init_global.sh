#!/usr/bin/env bash
#
# init_global.sh - generate a global air parcel initialization file for the
# domain-decomposition benchmark using the repository's own atm_init binary.
#
# The file reproduces the run.sh INIT_EVENLY=1 setup on a 1x1 global grid:
# node-centered lon/lat lattice over the full globe with cosine-weighted
# latitudes, uniform jitter, and a single pressure level (zeta).  Particle
# identities are deterministic through GSL_RNG_SEED.  Downstream, split the
# global file into per-subdomain init.bin files with init_partition.py and
# feed them to benchmark.sh via DD_BENCH_INIT_SOURCE.
#
# Usage:
#   ./init_global.sh <out_file> <lon_stride> <lat_stride>
#
# The strides are the per-subdomain lattice spacing of a 1x1 global domain
# (i.e. the total grid), so <lon_stride>*<lat_stride> is the exact global
# particle count.  Any positive integers are accepted; benchmark.sh already
# derives its grid as a multiple of every layout dimension.
#
# Environment overrides:
#   DD_INIT_SEED        GSL_RNG_SEED for reproducibility (default: 20261029)
#   DD_INIT_Z0/Z1       pressure level range [hPa] (defaults: 30.0/30.0)
#   DD_INIT_DZ          pressure step [hPa] (default: 1.0)
#   DD_INIT_UZ          vertical jitter half-width [hPa] (default: 60)
#   DD_INIT_ULON        longitude jitter half-width [deg] (default: 360)
#   DD_INIT_ULAT        latitude jitter half-width [deg] (default: 180)
#   DD_INIT_T0          start time [s]; default: 2022-06-02 00:00 UTC
#   DD_INIT_COMPILER    compiler for atm_init (default: gcc)
#   DD_INIT_OPT         optimization flags (default: -O2)
#   DD_INIT_MAKE_J      parallel make jobs (default: 4)
#
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <out_file> <lon_stride> <lat_stride>" >&2
    exit 1
fi

out_file=$1
lon_stride=$2
lat_stride=$3

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mptrac_dir=$(cd "$script_dir/../.." && pwd)
export LD_LIBRARY_PATH="$mptrac_dir/libs/build/lib:${LD_LIBRARY_PATH:-}"

if [[ ! $lon_stride =~ ^[0-9]+$ ]] || [[ ! $lat_stride =~ ^[0-9]+$ ]] ||
   (( lon_stride < 1 || lat_stride < 1 )); then
    echo "[ERROR] Strides must be positive integers (got $lon_stride $lat_stride)" >&2
    exit 1
fi

total=$((lon_stride * lat_stride))
seed=${DD_INIT_SEED:-20261029}
z0=${DD_INIT_Z0:-30.0}
z1=${DD_INIT_Z1:-30.0}
dz=${DD_INIT_DZ:-1.0}
uz=${DD_INIT_UZ:-60}
ulon=${DD_INIT_ULON:-360}
ulat=${DD_INIT_ULAT:-180}
compiler=${DD_INIT_COMPILER:-gcc}
opt=${DD_INIT_OPT:--O2}
make_j=${DD_INIT_MAKE_J:-4}

# Replicate the run.sh initialization box for a 1x1 global domain:
# node-centered lattice with a quarter-step tail padding so that the decimal
# step truncation cannot silently omit a whole particle row.
dlon=$(echo "scale=6; 360.0 / $lon_stride" | bc)
dlat=$(echo "scale=6; 180.0 / $lat_stride" | bc)
lon_offset=$(echo "scale=6; $dlon / 2" | bc)
lat_offset=$(echo "scale=6; $dlat / 2" | bc)
init_lon0=$(echo "scale=6; 0 + $lon_offset" | bc)
init_lon1=$(echo "scale=6; 360 - $lon_offset + $dlon / 4" | bc)
init_lat0=$(echo "scale=6; -90 + $lat_offset" | bc)
init_lat1=$(echo "scale=6; 90 - $lat_offset + $dlat / 4" | bc)

echo "[INFO] Global init: $lon_stride x $lat_stride = $total particles"
echo "[INFO]   Lon box: [$init_lon0, $init_lon1] step $dlon"
echo "[INFO]   Lat box: [$init_lat0, $init_lat1] step $dlat"
echo "[INFO]   zeta: [$z0, $z1] jitter uz=$uz ulon=$ulon ulat=$ulat seed=$seed"

build_dir=$(mktemp -d "${TMPDIR:-/tmp}/init_global_build.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT

met_nx=${DD_MET_NX:-2880}
met_ny=${DD_MET_NY:-1441}
met_nz=${DD_MET_NZ:-137}
met_ep=${DD_MET_EP:-140}

echo "[INFO] Building atm_init with NP=$total"
make -C "$mptrac_dir/src" clean > "$build_dir/build.log" 2>&1
make -C "$mptrac_dir/src" -j"$make_j" atm_init time2jsec STATIC=0 GPU=0 \
    COMPILER="$compiler" OPT="$opt" \
    DEFINES="-DNP=$total -DEX=$((met_nx + 4)) -DEY=$((met_ny + 5)) -DEP=$met_ep" \
    >> "$build_dir/build.log" 2>&1 || {
        cat "$build_dir/build.log" >&2
        echo "[ERROR] Failed to build atm_init" >&2
        exit 1
    }

t0=${DD_INIT_T0:-$("$mptrac_dir"/src/time2jsec 2022 6 2 0 0 0 0)}

cat > "$build_dir/config.ctl" <<EOF
MET_PRESS_LEVEL_DEF = -1
MET_CLAMS = 0
MET_VERT_COORD = 0
MET_DT_OUT = 3600
MET_CONVENTION = 0
MET_NC_READ_MODE = 0
ATM_TYPE = 1
ATM_TYPE_OUT = 1
ADVECT = 2
ADVECT_VERT_COORD = 0
DIFFUSION = 0
TURB_DX_TROP = 0
TURB_DX_STRAT = 0
TURB_DZ_TROP = 0
TURB_DZ_STRAT = 0.0
TURB_MESOX = 0.0
TURB_MESOZ = 0.0
DIRECTION = 1
METBASE = wind
MET_TROPO = 3
TDEC_TROP = 259200
TDEC_STRAT = 259200
DT_MOD = 600
DT_MET = 3600
T_START = $t0
T_STOP = $t0
ATM_DT_OUT = 86400
GRID_DT_OUT = 25920000
GRID_LON0 = 0
GRID_LON1 = 360
GRID_LAT0 = -90
GRID_LAT1 = 90
GRID_NX = 370
GRID_NY = 190
ENSEMBLE = 1
NQ = 5
QNT_NAME[0] = idx
QNT_NAME[1] = zeta
QNT_NAME[2] = m
QNT_NAME[3] = current_subdomain
QNT_NAME[4] = target_subdomain
ISOSURF = 0
DD_HALOS_SIZE = 2
DD_SUBDOMAINS_MERIDIONAL = 1
DD_SUBDOMAINS_ZONAL = 1
DD_EXCHANGE_MODE = 0
DD_SORT_DT = 1800
NAP_MAX = $total
EOF

mkdir -p "$(dirname "$out_file")"

echo "[INFO] Generating $out_file"
GSL_RNG_SEED=$seed "$mptrac_dir"/src/atm_init "$build_dir"/config.ctl "$out_file" \
    INIT_T0 "$t0" INIT_T1 "$t0" \
    INIT_Z0 "$z0" INIT_Z1 "$z1" INIT_DZ "$dz" \
    INIT_UZ "$uz" INIT_ULON "$ulon" INIT_ULAT "$ulat" \
    INIT_LON0 "$init_lon0" INIT_LON1 "$init_lon1" INIT_DLON "$dlon" \
    INIT_LAT0 "$init_lat0" INIT_LAT1 "$init_lat1" INIT_DLAT "$dlat" \
    INIT_EVENLY 1 INIT_IDX_OFFSET 0 \
    > "$build_dir/atm_init.log" 2>&1 || {
        cat "$build_dir/atm_init.log" >&2
        echo "[ERROR] atm_init failed" >&2
        exit 1
    }

# The binary header carries the produced particle count; it must match the
# requested global count exactly.
produced=$(python3 - "$out_file" <<'PY'
import struct
import sys
with open(sys.argv[1], "rb") as stream:
    version, count = struct.unpack("=ii", stream.read(8))
if version != 100:
    raise SystemExit(f"unexpected init version {version}")
print(count)
PY
)
if (( produced != total )); then
    echo "[ERROR] atm_init produced $produced particles, expected $total" >&2
    exit 1
fi

echo "[INFO] Verified particle count: $produced"
echo "[INFO] Wrote $out_file"
