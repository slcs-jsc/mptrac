#! /bin/bash

# Testing components for interoperable use between CLaMS and MPTRAC...

# ------------------------------------------------------------
# (1) Test reading and writing CLaMS position files...
# ------------------------------------------------------------

# Setup...
trac=../../src

# Create directories...
rm -rf data && mkdir -p data

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH

# Run the application...
cat > ./data/conv.ctl <<EOF
NQ = 5
QNT_NAME[0] = theta
QNT_NAME[1] = pv
QNT_NAME[2] = m
QNT_NAME[3] = zeta
QNT_NAME[4] = zeta_d
EOF
$trac/atm_conv ./data/conv.ctl ./data.ref/atm_input.tab 0 ./data/atm_output.nc 4
$trac/atm_conv - ./data/atm_output.nc 4 ./data/atm_output.tab 0
if $trac/atm_conv - ./data.ref/atm_input.tab 0 ./data/invalid_clams.nc 4 \
    > ./data/invalid_clams.log 2>&1; then
  echo 'Accepted CLaMS position output without zeta_d' >&2
  exit 1
fi
grep -q 'requires quantity zeta_d' ./data/invalid_clams.log

# ------------------------------------------------------------
# (2) Test diabatic transport calculations...
# ------------------------------------------------------------

# start and end date
t0=$(${trac}/time2jsec 2016 7 1 0 0 0 0)
t1=$(${trac}/time2jsec 2016 7 1 6 0 0 0)

# Create control parameter file...
cat > ./data/trac.ctl <<EOF
MET_CONVENTION = 1
MET_PRESS_LEVEL_DEF = 5
ATM_TYPE = 3
ATM_TYPE_OUT = 0
ADVECT = 2
MET_CLAMS = 1
ADVECT_VERT_COORD = 1
MET_VERT_COORD = 1
NQ = 7
QNT_NAME[0] = theta
QNT_NAME[1] = pv
QNT_NAME[2] = m
QNT_NAME[3] = zeta
QNT_NAME[4] = zeta_d
QNT_NAME[5] = ps
QNT_NAME[6] = p
METBASE = ../data/erai_vlr
DIRECTION = 1
MET_TROPO = 3
TDEC_TROP = 259200
TDEC_STRAT = 259200
DT_MOD = 180
DT_MET = 21600   
T_START = ${t0}         
T_STOP = ${t1}
CHUNKSZHINT = 163840000
ATM_DT_OUT = 21600
EOF

# Calculate trajectories...
echo "./data" > ./data/dirlist
$trac/trac ./data/dirlist trac.ctl ../data.ref/init/pos_glo_16070100.nc \
	   ATM_BASENAME atm GRID_BASENAME grid

# Compare files...
echo -e "\nCompare results..."
error=0
diff -q -s ./data.ref/atm_output.tab ./data/atm_output.tab || error=1
diff -q -s ./data.ref/atm_2016_07_01_00_00_00.tab ./data/atm_2016_07_01_00_00_00.tab || error=1
diff -q -s ./data.ref/atm_2016_07_01_06_00_00.tab ./data/atm_2016_07_01_06_00_00.tab || error=1

# ------------------------------------------------------------
# (3) Check output isolation across independent directory-list tasks...
# ------------------------------------------------------------

set -euo pipefail

trac=../../src
base=data/output_safety
mkdir -p "$base"
export OMP_NUM_THREADS=1

# Two short synthetic tasks share one executable but have distinct observations.
"$trac/wind" - "$base/met" WIND_T0 0 WIND_NX 3 WIND_NY 3 WIND_NZ 3 > /dev/null
"$trac/wind" - "$base/met" WIND_T0 3600 WIND_NX 3 WIND_NY 3 WIND_NZ 3 > /dev/null
printf '0 1 0 0 2 0 300 300\n0 1 0 0 3 1 300 300\n' > "$base/atm.tab"
printf '0 1 0 0 1\n' > "$base/obs_a.tab"
printf '0 1 0 0 1\n0 1 0 0 2\n' > "$base/obs_b.tab"
printf '0 1\n2 0\n' > "$base/kernel.tab"

for name in a_alone b_alone a_seq b_seq a_rev b_rev; do
  mkdir -p "$base/$name"
  case "$name" in a_*) task=a;; *) task=b;; esac
  if [[ $task == a ]]; then
    kernel="$base/kernel.tab"
  else
    kernel=-
  fi
  cat > "$base/$name/trac.ctl" <<EOF
NQ = 4
QNT_NAME[0] = m
QNT_NAME[1] = ens
QNT_NAME[2] = zeta_d
QNT_NAME[3] = zeta
NENS = 2
MOLMASS = 64
METBASE = $base/met
DT_MET = 3600
DT_MOD = 180
T_STOP = 180
MET_PBL = 0
MET_TROPO = 0
CONV_CAPE = -1
ADVECT = 1
DIFFUSION = 0
ATM_TYPE = 0
ATM_TYPE_OUT = 3
ATM_BASENAME = atm
ATM_DT_OUT = 180
GRID_BASENAME = grid
GRID_DT_OUT = 180
GRID_KERNEL = $kernel
GRID_LON0 = -1
GRID_LON1 = 1
GRID_LAT0 = -1
GRID_LAT1 = 1
GRID_NX = 1
GRID_NY = 1
GRID_Z0 = 0
GRID_Z1 = 2
GRID_NZ = 1
ENS_BASENAME = ens
ENS_DT_OUT = 180
BUDGET_BASENAME = budget
BUDGET_DT_OUT = 180
CSI_BASENAME = csi
CSI_DT_OUT = 180
CSI_OBSFILE = $base/obs_$task.tab
CSI_KERNEL = $kernel
CSI_LON0 = -1
CSI_LON1 = 1
CSI_LAT0 = -1
CSI_LAT1 = 1
CSI_NX = 1
CSI_NY = 1
CSI_Z0 = 0
CSI_Z1 = 2
CSI_NZ = 1
PROF_BASENAME = prof
PROF_OBSFILE = $base/obs_$task.tab
PROF_LON0 = -1
PROF_LON1 = 1
PROF_LAT0 = -1
PROF_LAT1 = 1
PROF_NX = 1
PROF_NY = 1
PROF_Z0 = 0
PROF_Z1 = 2
PROF_NZ = 1
SAMPLE_BASENAME = sample
SAMPLE_OBSFILE = $base/obs_$task.tab
SAMPLE_KERNEL = $kernel
SAMPLE_DX = 50
SAMPLE_DZ = 2
EOF
  cp "$base/atm.tab" "$base/$name/atm.tab"
done

run_list() {
  local list=$1
  shift
  printf '%s\n' "$@" > "$base/$list"
  "$trac/trac" "$base/$list" trac.ctl atm.tab > "$base/$list.log" 2>&1
}
run_list alone_a "$base/a_alone"
run_list alone_b "$base/b_alone"
run_list sequence "$base/a_seq" "$base/b_seq"
run_list reverse "$base/b_rev" "$base/a_rev"

for task in a b; do
  for variant in seq rev; do
    for file in ens_2000_01_01_00_00_00.tab budget.tab csi.tab prof.tab sample.tab \
      grid_2000_01_01_00_00_00.tab grid_2000_01_01_00_03_00.tab \
      traj_fix_3d_00010100_00010100.nc init_fix_00010100.nc; do
      cmp "$base/${task}_alone/$file" "$base/${task}_$variant/$file"
    done
  done
done

# The two-particle task above already covers the first and last ensemble IDs.
for value in 0.5 -1 2 nan inf; do
  printf '0 1 0 0 2 %s 300 300\n' "$value" > "$base/b_alone/atm.tab"
  printf '%s\n' "$base/b_alone" > "$base/invalid_id"
  if "$trac/trac" "$base/invalid_id" trac.ctl atm.tab \
      ATM_BASENAME - BUDGET_BASENAME - CSI_BASENAME - \
      PROF_BASENAME - SAMPLE_BASENAME - > "$base/invalid_id.log" 2>&1; then
    echo "Accepted invalid ensemble ID: $value" >&2
    exit 1
  fi
  grep -q 'Ensemble ID' "$base/invalid_id.log"
done

# Budget and CSI must also reject fractional IDs before integer conversion.
printf '0 1 0 0 2 0.5 300 300\n' > "$base/b_alone/atm.tab"
for output in budget csi; do
  budget=-
  csi=-
  if [[ $output == budget ]]; then budget=budget; else csi=csi; fi
  if "$trac/trac" "$base/invalid_id" trac.ctl atm.tab \
      ATM_BASENAME - GRID_BASENAME - ENS_BASENAME - \
      BUDGET_BASENAME "$budget" CSI_BASENAME "$csi" \
      PROF_BASENAME - SAMPLE_BASENAME - > "$base/invalid_$output.log" 2>&1; then
    echo "Accepted fractional ensemble ID in $output" >&2
    exit 1
  fi
  grep -q 'Ensemble ID' "$base/invalid_$output.log"
done

# A missing CLaMS coordinate must fail before dereferencing a quantity index.
printf '0 1 0 0 2 0 300 300\n' > "$base/b_alone/atm.tab"
sed -i 's/^QNT_NAME\[2\] = zeta_d$/QNT_NAME[2] = theta/' \
  "$base/b_alone/trac.ctl"
if "$trac/trac" "$base/invalid_id" trac.ctl atm.tab \
    ENS_BASENAME - BUDGET_BASENAME - CSI_BASENAME - \
    PROF_BASENAME - SAMPLE_BASENAME - > "$base/missing_zeta.log" 2>&1; then
  echo 'Accepted missing CLaMS zeta_d' >&2
  exit 1
fi
grep -q 'requires quantity zeta_d' "$base/missing_zeta.log"

echo 'Output task isolation and validation passed.'

# Check sample VMR against a uniform 1 ppb tracer in an interior layer and
# a layer clipped by the synthetic meteorology's flat surface.
sample_dir="$base/sample_vmr"
mkdir -p "$sample_dir"
cp "$base/a_seq/trac.ctl" "$sample_dir/trac.ctl"
sample_mass() {
  awk -v z="$1" 'BEGIN {
    pi = atan2(0, -1)
    lower = z - 0.5
    if (lower < 0) lower = 0
    p_lower = 1013.25 * exp(-lower / 7)
    p_upper = 1013.25 * exp(-(z + 0.5) / 7)
    air_column = 100 * (p_lower - p_upper) / 9.80665
    printf "%.17g", 1e-9 * (64 / 28.9644) * air_column * pi * 50 * 50 * 1e6
  }'
}
mass_interior=$(sample_mass 10)
mass_surface=$(sample_mass 0.25)
printf '0 10 0 0 %s 0 300 300\n0 0.25 10 0 %s 1 300 300\n' \
  "$mass_interior" "$mass_surface" > "$sample_dir/atm.tab"
printf '0 10 0 0 1\n0 0.25 10 0 1\n' > "$sample_dir/obs.tab"
printf '%s\n' "$sample_dir" > "$base/sample_list"
"$trac/trac" "$base/sample_list" trac.ctl atm.tab \
  ATM_BASENAME - GRID_BASENAME - ENS_BASENAME - BUDGET_BASENAME - \
  CSI_BASENAME - PROF_BASENAME - SAMPLE_OBSFILE "$sample_dir/obs.tab" \
  SAMPLE_KERNEL - SAMPLE_DZ 0.5 > "$sample_dir/run.log" 2>&1
awk -v m1="$mass_interior" -v m2="$mass_surface" '
  NF == 10 && $1 !~ /^#/ {
    n++
    mass = (n == 1 ? m1 : m2)
    cd = mass / (atan2(0, -1) * 50 * 50 * 1e6)
    if ($7 != 1 || ($8 - cd) / cd > 1e-5 || (cd - $8) / cd > 1e-5 ||
        ($9 - 1e-9) / 1e-9 > 1e-5 || (1e-9 - $9) / 1e-9 > 1e-5) exit 1
  }
  END { if (n != 2) exit 1 }
' "$sample_dir/sample.tab"

# A staggered, unreleased first particle must not define snapshot time.
clock_dir="$base/clams_clock"
mkdir -p "$clock_dir"
cp "$base/a_seq/trac.ctl" "$clock_dir/trac.ctl"
printf '360 1 0 0 2 0 300 300\n0 1 0 0 3 1 300 300\n' > "$clock_dir/atm.tab"
printf '%s\n' "$clock_dir" > "$base/clock_list"
"$trac/trac" "$base/clock_list" trac.ctl atm.tab \
  GRID_BASENAME - ENS_BASENAME - BUDGET_BASENAME - CSI_BASENAME - \
  PROF_BASENAME - SAMPLE_BASENAME - > "$clock_dir/run.log" 2>&1
if command -v ncdump > /dev/null; then
  ncdump -v time "$clock_dir/traj_fix_3d_00010100_00010100.nc" |
    awk '/^ time =/ { gsub(/[,;]/, ""); if ($3 != 0 || $4 != 180) exit 1; found = 1 }
         END { if (!found) exit 1 }'
  ncdump -v time "$clock_dir/init_fix_00010100.nc" |
    awk '/^ time =/ { gsub(/[,;]/, ""); if ($3 != 180) exit 1; found = 1 }
         END { if (!found) exit 1 }'
fi

# The shortened final step must assign each observation once, in either
# trajectory direction, and CSI must flush the last partial interval.
check_windows() {
  local name=$1 start=$2 stop=$3 first=$4 second=$5 direction=$6
  local dir="$base/window_$name"
  mkdir -p "$dir"
  cp "$base/a_seq/trac.ctl" "$dir/trac.ctl"
  printf '%s 1 0 0 5 0 300 300\n' "$start" > "$dir/atm.tab"
  printf '%s 1 0 0 1\n%s 1 0 0 1\n' "$first" "$second" > "$dir/obs.tab"
  printf '%s\n' "$dir" > "$base/window_list"
  "$trac/trac" "$base/window_list" trac.ctl atm.tab \
    T_STOP "$stop" DT_MOD 60 DIRECTION "$direction" NENS 0 \
    ATM_BASENAME - GRID_BASENAME - ENS_BASENAME - BUDGET_BASENAME - \
    CSI_BASENAME csi CSI_DT_OUT 60 CSI_OBSFILE "$dir/obs.tab" \
    PROF_BASENAME prof PROF_OBSFILE "$dir/obs.tab" \
    SAMPLE_OBSFILE "$dir/obs.tab" CSI_KERNEL - SAMPLE_KERNEL - \
    > "$dir/run.log" 2>&1
  awk -v a="$first" -v b="$second" '
    NF == 10 && $1 !~ /^#/ { n++; if ($1 == a) na++; if ($1 == b) nb++ }
    END { if (n != 2 || na != 1 || nb != 1) exit 1 }
  ' "$dir/sample.tab"
  awk -v stop="$stop" '
    NF == 20 && $1 !~ /^#/ { n++; if ($1 == stop) final++ }
    END { if (n != 2 || final != 1) exit 1 }
  ' "$dir/csi.tab"
  awk '$1 !~ /^#/ && NF == 11 { n += $11 } END { if (n != 2) exit 1 }' \
    "$dir/prof.tab"
}
check_windows forward 0 100 75 85 1
check_windows backward 120 20 35 45 -1

echo 'Output units, snapshot time, and partial-step windows passed.'

# Check central moments using large-offset values and unequal kernel weights.
stats_dir="$base/welford"
mkdir -p "$stats_dir"
cat > "$stats_dir/trac.ctl" <<EOF
NQ = 3
QNT_NAME[0] = m
QNT_NAME[1] = ens
QNT_NAME[2] = marker
QNT_UNIT[2] = s
QNT_FORMAT[2] = %.12g
NENS = 1
MOLMASS = 64
METBASE = $base/met
DT_MET = 3600
DT_MOD = 180
T_START = 0
T_STOP = 180
MET_PBL = 0
MET_TROPO = 0
CONV_CAPE = -1
ADVECT = 1
DIFFUSION = 0
ATM_TYPE = 0
ATM_BASENAME = -
GRID_BASENAME = grid
GRID_DT_OUT = 180
GRID_STDDEV = 1
GRID_LON0 = -1
GRID_LON1 = 1
GRID_LAT0 = -1
GRID_LAT1 = 1
GRID_NX = 1
GRID_NY = 1
GRID_Z0 = 0
GRID_Z1 = 2
GRID_NZ = 1
ENS_BASENAME = ens
ENS_DT_OUT = 180
EOF
printf '%s\n' "$stats_dir" > "$base/stats_list"
printf '0 1 0 0 2 0 999999998\n0 1 0 0 2 0 999999999\n0 1 0 0 2 0 1000000000\n0 1 0 0 2 0 1000000001\n0 1 0 0 2 0 1000000002\n' > "$stats_dir/atm.tab"
"$trac/trac" "$base/stats_list" trac.ctl atm.tab > "$stats_dir/run.log" 2>&1
awk '!/^#/ && NF { if ($9 != 5 || $12 != 1e9 ||
    ($15 - sqrt(2))^2 > 1e-18) exit 1; found = 1 }
    END { if (!found) exit 1 }' "$stats_dir/grid_2000_01_01_00_00_00.tab"
awk '!/^#/ && NF { if ($11 != 5 || $7 != 1e9 ||
    ($10 - sqrt(2))^2 > 1e-18) exit 1; found = 1 }
    END { if (!found) exit 1 }' "$stats_dir/ens_2000_01_01_00_00_00.tab"

printf '0.5 1\n1.5 2\n' > "$stats_dir/kernel.tab"
printf 'GRID_KERNEL = %s\n' "$stats_dir/kernel.tab" >> "$stats_dir/trac.ctl"
printf '0 0.5 0 0 2 0 999999999\n0 1.5 0 0 4 0 1000000002\n' > "$stats_dir/atm.tab"
"$trac/trac" "$base/stats_list" trac.ctl atm.tab > "$stats_dir/run.log" 2>&1
awk '!/^#/ && NF { cd = 5 / (1e6 * $5);
    if ($9 != 2 || $12 != 1000000001 ||
        ($15 - sqrt(2))^2 > 1e-18 ||
        ($7 - cd)^2 / (cd * cd) > 1e-10) exit 1; found = 1 }
    END { if (!found) exit 1 }' "$stats_dir/grid_2000_01_01_00_00_00.tab"
awk '!/^#/ && NF { if ($11 != 2 || $7 != 1000000000.5 ||
    ($10 - 1.5)^2 > 1e-18) exit 1; found = 1 }
    END { if (!found) exit 1 }' "$stats_dir/ens_2000_01_01_00_00_00.tab"

# Two observations with a 1 kg/m2 spread have unit population deviation.
cat >> "$stats_dir/trac.ctl" <<EOF
CSI_BASENAME = csi
CSI_DT_OUT = 180
CSI_OBSFILE = $stats_dir/obs.tab
CSI_KERNEL = -
CSI_LON0 = -1
CSI_LON1 = 1
CSI_LAT0 = -1
CSI_LAT1 = 1
CSI_NX = 1
CSI_NY = 1
CSI_Z0 = 0
CSI_Z1 = 2
CSI_NZ = 1
EOF
printf '0 1 0 0 999999999\n0 1 0 0 1000000001\n' > "$stats_dir/obs.tab"
printf '0 1 0 0 0 0 1000000000\n' > "$stats_dir/atm.tab"
"$trac/trac" "$base/stats_list" trac.ctl atm.tab > "$stats_dir/run.log" 2>&1
awk '!/^#/ && NF { if ($19 != -5e17 || $20 != 1) exit 1; found = 1 }
    END { if (!found) exit 1 }' "$stats_dir/csi.tab"

# A zero-weight particle still counts, but cannot change weighted moments.
printf '0.5 0\n1.5 1\n' > "$stats_dir/zero_kernel.tab"
printf '0 0.5 0 0 2 0 1000000000\n0 1.5 0 0 4 0 1000000000\n' \
  > "$stats_dir/atm.tab"
"$trac/trac" "$base/stats_list" trac.ctl atm.tab \
  CSI_BASENAME - GRID_NX 2 GRID_KERNEL "$stats_dir/zero_kernel.tab" \
  > "$stats_dir/run.log" 2>&1
awk '!/^#/ && NF { n++; if (n == 1) {
      if ($9 != 0 || $12 != "nan") exit 1
    } else if ($9 != 2 || $12 != 1e9 || $15 != 0) exit 1 }
    END { if (n != 2) exit 1 }' "$stats_dir/grid_2000_01_01_00_00_00.tab"
awk '!/^#/ && NF { if ($11 != 2 || $7 != 1e9 || $10 != 0) exit 1;
    found = 1 } END { if (!found) exit 1 }' \
  "$stats_dir/ens_2000_01_01_00_00_00.tab"
echo 'Large-offset, weighted, and observation statistics passed.'
exit $error
