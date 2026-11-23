#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4
export LANG=C
export LC_ALL=C

# Setup...
trac=../../src

# Create directories...
rm -rf data && mkdir -p data

# Set timestep and timerange of simulation...
t0=$($trac/time2jsec 2011 6 5 0 0 0 0)
t1=$($trac/time2jsec 2011 6 5 0 1 0 0)

# Create control parameter file...
cat > data/trac.ctl <<EOF
NQ = 4
QNT_NAME[0] = t
QNT_NAME[1] = u
QNT_NAME[2] = v
QNT_NAME[3] = w
METBASE = ../data/ei
ATM_DT_OUT = 10.0
DT_MOD = 10.0
DIFFUSION = 1
DT_MET = 86400.0
T_STOP = $t1
EOF

# Set initial air parcel positions...
$trac/atm_init data/trac.ctl data/atm_init.tab \
	       INIT_T0 "$t0" INIT_T1 "$t0" \
	       INIT_Z0 10.0 INIT_Z1 10.0 \
	       INIT_LON0 -72.117 INIT_LON1 -72.117 \
	       INIT_LAT0 -40.59 INIT_LAT1 -40.59

# Split air parcels...
$trac/atm_split data/trac.ctl data/atm_init.tab data/atm_split.tab \
		SPLIT_N 10000 SPLIT_M 1e9 SPLIT_DX 30.0 SPLIT_DZ 5.0

# Calculate trajectories on pressure levels...
echo "data" > data/dirlist
$trac/trac data/dirlist trac.ctl atm_split.tab  ATM_BASENAME atm_pl

# Convert to grid data...
$trac/atm2grid data/trac.ctl data/atm_pl_2011_06_05_00_00_00.tab \
		GRID_BASENAME data/grid GRID_NX 72 GRID_NY 36 MOLMASS 64.066

$trac/atm_stat data/trac.ctl data/stat_mean.tab mean \
  data/atm_pl_2011_06_05_00_00_00.tab

$trac/atm_dist data/trac.ctl data/dist_mean.tab mean \
 data/atm_pl_2011_06_05_00_00_00.tab \
 data/atm_pl_2011_06_05_00_00_10.tab


# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done

# Enabled periodic outputs require a finite, positive interval. Disabled
# basenames may retain unused interval values.
for pair in \
  ATM_BASENAME:ATM_DT_OUT BUDGET_BASENAME:BUDGET_DT_OUT \
  DEPO_BASENAME:DEPO_DT_OUT CSI_BASENAME:CSI_DT_OUT \
  ENS_BASENAME:ENS_DT_OUT GRID_BASENAME:GRID_DT_OUT \
  VTK_BASENAME:VTK_DT_OUT; do
    basename=${pair%:*}
    interval=${pair#*:}
    for invalid in 0 inf; do
      if $trac/atm_init data/trac.ctl data/interval.tab \
          "$basename" probe "$interval" "$invalid" \
          > data/interval.log 2>&1; then
        echo "Accepted $interval=$invalid with $basename enabled" >&2
        error=1
      fi
      grep -q "$interval must be finite and positive" data/interval.log || error=1
    done
done
for invalid in -1 nan; do
    if $trac/atm_init data/trac.ctl data/interval.tab \
        ATM_BASENAME probe ATM_DT_OUT "$invalid" \
        > data/interval.log 2>&1; then
      echo "Accepted ATM_DT_OUT=$invalid" >&2
      error=1
    fi
    grep -q 'ATM_DT_OUT must be finite and positive' data/interval.log || error=1
done
$trac/atm_init data/trac.ctl data/interval.tab \
    ATM_BASENAME - ATM_DT_OUT 0 CSI_BASENAME - CSI_DT_OUT 0 \
    MET_DT_OUT 0 > data/interval.log 2>&1 || error=1
if $trac/atm_init data/trac.ctl data/interval.tab \
    MET_DT_OUT nan > data/interval.log 2>&1; then
    echo 'Accepted nonfinite MET_DT_OUT with meteo quantities' >&2
    error=1
fi
grep -q 'MET_DT_OUT must be finite' data/interval.log || error=1
exit $error
