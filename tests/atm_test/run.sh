#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4

# Setup...
trac=../../src

# Create directories...
rm -rf data && mkdir -p data

# Init...
qnt="NQ 3 QNT_NAME[0] aoa QNT_NAME[1] m QNT_NAME[2] vmr"
$trac/atm_init - data/atm_2000_01_01_00_00.tab \
	       INIT_Z0 30 INIT_Z1 30 INIT_UZ 60 \
	       INIT_ULON 360 INIT_ULAT 180 \
	       INIT_EVENLY 1 INIT_REP 10000 INIT_MASS 1e9 INIT_VMR 1e-9 $qnt

$trac/atm_init - data/atm2_2000_01_01_00_00.tab \
	       INIT_Z0 30 INIT_Z1 30 INIT_UZ 60 \
	       INIT_ULON 360 INIT_ULAT 180 \
	       INIT_EVENLY 1 INIT_REP 10000 INIT_MASS 2e9 INIT_VMR 2e-9 $qnt

$trac/atm_init - data/atm_2000_01_02_00_00.tab \
	       INIT_T0 86400 INIT_T1 86400 \
	       INIT_Z0 30 INIT_Z1 30 INIT_UZ 60 \
	       INIT_ULON 360 INIT_ULAT 180 \
	       INIT_EVENLY 1 INIT_REP 10000 INIT_MASS 3e9 INIT_VMR 3e-9 $qnt

export GSL_RNG_SEED=123
$trac/atm_init - data/atm2_2000_01_02_00_00.tab \
	       INIT_T0 86400 INIT_T1 86400 \
	       INIT_Z0 30 INIT_Z1 30 INIT_UZ 60 \
	       INIT_ULON 360 INIT_ULAT 180 \
	       INIT_EVENLY 1 INIT_REP 10000 INIT_MASS 4e9 INIT_VMR 4e-9 $qnt

# Convert...
$trac/atm_conv - data/atm_2000_01_01_00_00.tab 0 data/atm.bin 1 $qnt
$trac/atm_conv - data/atm.bin 1 data/atm_ascii_from_bin.tab 0 $qnt
diff -s data/atm_ascii_from_bin.tab data/atm_2000_01_01_00_00.tab || exit 1

$trac/atm_conv - data/atm_2000_01_01_00_00.tab 0 data/atm.nc 2 $qnt
$trac/atm_conv - data/atm.nc 2 data/atm_ascii_from_nc.tab 0 $qnt
diff -s data/atm_ascii_from_nc.tab data/atm_2000_01_01_00_00.tab || exit 1

# Select...
$trac/atm_select - data/atm_select.tab data/atm_2000_01_01_00_00.tab \
		 SELECT_IP0 1000 SELECT_IP1 9000 SELECT_STRIDE 2 \
		 SELECT_Z0 10 SELECT_Z1 50 \
		 SELECT_LON0 -20 SELECT_LON1 60 \
		 SELECT_LAT0 -30 SELECT_LAT1 50 \
		 SELECT_RLON 20 SELECT_RLAT 10 \
		 SELECT_R0 2000 SELECT_R1 5000 $qnt

# Statistics...
for param in mean stddev min max skew kurt median absdev mad ; do
    $trac/atm_stat - data/stat_$param.tab $param \
		   data/atm_2000_01_01_00_00.tab \
		   STAT_Z0 5 STAT_Z1 55 \
		   STAT_LAT0 -75 STAT_LAT1 75 \
		   STAT_LON0 -150 STAT_LON1 150 $qnt
done

# Distances...
#for param in mean stddev min max skew kurt median absdev mad ; do
# TODO: disabled tests of stddev, skew, and kurt because of numerical differences with different compilers
for param in mean min max median absdev mad ; do
    $trac/atm_dist - data/dist_$param.tab $param \
		   data/atm_2000_01_01_00_00.tab \
		   data/atm2_2000_01_01_00_00.tab \
		   data/atm_2000_01_02_00_00.tab \
		   data/atm2_2000_01_02_00_00.tab \
		   DIST_Z0 5 DIST_Z1 55 \
		   DIST_LAT0 -75 DIST_LAT1 75 \
		   DIST_LON0 -150 DIST_LON1 150 $qnt
done

# Convert to grid data...
$trac/atm2grid - data/atm_2000_01_01_00_00.tab \
	       GRID_BASENAME data/grid GRID_NX 72 GRID_NY 36 $qnt MOLMASS 64.066
$trac/atm2grid - data/atm_2000_01_01_00_00.tab \
	       GRID_BASENAME data/grid GRID_NX 72 GRID_NY 36 $qnt MOLMASS 64.066 GRID_TYPE 1

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
