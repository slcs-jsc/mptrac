#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4

# Setup...
trac=../../src

# Create directories...
rm -rf data && mkdir -p data

# Create tropopause data file...
$trac/tropo - data/tropo_2011_06_05.nc ../data/ei_2011_06_05_00.nc ../data/ei_2011_06_06_00.nc TROPO_DLON 6 TROPO_DLAT 3

# Sample tropopause data...
t0=$($trac/time2jsec 2011 6 5 12 0 0 0)
$trac/atm_init - data/atm.tab \
	       INIT_T0 "$t0" INIT_T1 "$t0" \
	       INIT_LON0 -180 INIT_LON1 180 INIT_DLON 6 \
	       INIT_LAT0 -90 INIT_LAT1 90 INIT_DLAT 3
for var in clp dyn wmo_1st wmo_2nd ; do
    $trac/tropo_sample - data/sample_$var.tab data/tropo_2011_06_05.nc $var data/atm.tab
done

# Calculate zonal means...
for var in clp dyn wmo_1st wmo_2nd ; do
    $trac/tropo_zm - data/zm_$var.tab $var data/tropo_2011_06_05.nc
done

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*) ; do
    diff -q -s data/$(basename "$f") "$f"
    [ $? -ne 0 ] && error=1
done
exit $error
