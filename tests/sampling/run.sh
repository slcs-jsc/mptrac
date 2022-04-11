#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4

# Setup...
trac=../../src

# Create directories...
rm -rf data && mkdir -p data

# Sample map, profile, and zonal mean...
$trac/met_map - data/map.tab ../meteo/ei_2011_06_05_00.nc
$trac/met_prof - data/prof.tab ../meteo/ei_2011_06_05_00.nc
$trac/met_zm - data/zm.tab ../meteo/ei_2011_06_05_00.nc

# Set timerange for sampling...
tm=$($trac/time2jsec 2011 6 5 12 0 0 0)

# Create init file for sampling...
$trac/atm_init - data/atm_init.tab \
	       INIT_T0 $tm INIT_T1 $tm INIT_UT 86400 \
	       INIT_Z0 30 INIT_Z1 30 INIT_UZ 60 \
	       INIT_ULON 360 INIT_ULAT 180 \
	       INIT_EVENLY 1 INIT_REP 10000

# Sample meteo data...
$trac/met_sample - data/sample.tab data/atm_init.tab \
		 METBASE ../meteo/ei DT_MET 86400.0

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/$(basename "$f") "$f"
    [ $? -ne 0 ] && error=1
done
exit $error
