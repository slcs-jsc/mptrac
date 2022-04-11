#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4

# Setup...
trac=../../src

# Create directories...
rm -rf data && mkdir -p data

# Sample map, profile, zonal mean, and spectral data...
$trac/met_map - data/map.tab ../data/ei_2011_06_05_00.nc \
	      MAP_LON0 -180 MAP_LON1 180 MAP_DLON 6 \
	      MAP_LAT0 -90 MAP_LAT1 90 MAP_DLAT 4 \
	      MAP_THETA 380 TROPO 3
$trac/met_prof - data/prof.tab ../data/ei_2011_06_05_00.nc \
	       PROF_Z0 0 PROF_Z1 60 PROF_DZ 0.1 \
	       PROF_LON0 20 PROF_LON1 30 PROF_DLON 1 \
	       PROF_LAT0 40 PROF_LAT1 50 PROD_DLAT 1 \
	       TROPO 4
$trac/met_zm - data/zm.tab ../data/ei_2011_06_05_00.nc \
	     ZM_Z0 0 ZM_Z1 60 ZM_DZ 1 \
	     ZM_LON0 -180 ZM_LON1 180 ZM_DLON 1 \
	     ZM_LAT0 -90 ZM_LAT1 90 ZM_DLAT 3 \
	     TROPO 5
$trac/met_spec - data/spec.tab ../data/ei_2011_06_05_00.nc \
	       SPEC_WAVEMAX 3

# Create init file for sampling...
tm=$($trac/time2jsec 2011 6 5 12 0 0 0)
$trac/atm_init - data/atm_init.tab \
	       INIT_T0 "$tm" INIT_T1 "$tm" INIT_UT 86400 \
	       INIT_Z0 30 INIT_Z1 30 INIT_UZ 60 \
	       INIT_ULON 360 INIT_ULAT 180 \
	       INIT_EVENLY 1 INIT_REP 3000

# Sample meteo data...
$trac/met_sample - data/sample.tab data/atm_init.tab \
		 METBASE ../data/ei DT_MET 86400.0

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/$(basename "$f") "$f"
    [ $? -ne 0 ] && error=1
done
exit $error
