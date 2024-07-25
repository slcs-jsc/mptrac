#! /bin/bash

# ToDo: clean up

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export LANG=C
export LC_ALL=C

# Setup...
trac_fortran=../../src

# Create directories...
rm -rf data && mkdir -p data

# Set timestep and timerange of simulation...
t0=$($trac_fortran/time2jsec 2011 6 5 0 0 0 0)
t1=$($trac_fortran/time2jsec 2011 6 8 0 0 0 0)

# Create control parameter file...
cat > data/trac.ctl <<EOF
NQ = 13
QNT_NAME[0] = t
QNT_NAME[1] = u
QNT_NAME[2] = v
QNT_NAME[3] = w
QNT_NAME[4] = zg
QNT_NAME[5] = pv
QNT_NAME[6] = ps
QNT_NAME[7] = pt
QNT_NAME[8] = m
QNT_NAME[9] = stat
QNT_NAME[10] = ens
QNT_NAME[11] = Cccl3f
QNT_NAME[12] = Cx
METBASE = ../data/ei
MET_DT_OUT = 86400.0
SPECIES = SO2
BOUND_LAT0 = -90
BOUND_LAT1 = 90
BOUND_P0 = 1e10
BOUND_P1 = -1e10
BOUND_DPS = 100.0
BOUND_MASS = 0.0
CONV_CAPE = 0.0
H2O2_CHEM_REACTION = 1
TRACER_CHEM = 1
TDEC_TROP = 259200.0
TDEC_STRAT = 259200.0
DRY_DEPO_VDEP = 0.15
DRY_DEPO_DP = 300
MIXING_TROP = 1e-3
MIXING_STRAT = 1e-6
DT_MET = 86400.0
T_STOP = $t1
ATM_BASENAME = atm
ATM_TYPE_OUT = 0
CSI_OBSMIN = 1e-5
CSI_MODMIN = 1e-5
GRID_LON0 = -90
GRID_LON1 = 60
GRID_LAT0 = -60
GRID_LAT1 = -15
GRID_NX = 300
GRID_NY = 90
SAMPLE_DZ = 100
STAT_LON = -22
STAT_LAT = -40
EOF


## Set initial air parcel positions...
#$trac/atm_init data/trac.ctl data/atm_init.tab \
#	       INIT_T0 "$t0" INIT_T1 "$t0" \
#	       INIT_Z0 10.0 INIT_Z1 10.0 \
#	       INIT_LON0 -72.117 INIT_LON1 -72.117 \
#	       INIT_LAT0 -40.59 INIT_LAT1 -40.59
#
## Split air parcels...
#$trac/atm_split data/trac.ctl data/atm_init.tab data/atm_split.tab \
#		SPLIT_N 10000 SPLIT_M 1e9 SPLIT_DX 30.0 SPLIT_DZ 5.0

# Calculate trajectories...
echo "data" > data/dirlist
$trac_fortran/trac_fortran
#valgrind --tool=memcheck --leak-check=full --track-origins=yes -s $trac_fortran/trac_fortran
#valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all --track-origins=yes --max-stackframe=17081687808 --main-stacksize=90388608 -s $trac_fortran/trac_fortran
#valgrind --tool=memcheck --leak-check=full --track-origins=yes --max-stackframe=17081687808 --main-stacksize=190388608 -s $trac_fortran/trac_fortran

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab ) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
