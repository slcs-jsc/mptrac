#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4

# Setup...
trac=../../src

# Create directories...
rm -rf data && mkdir -p data

# Loop over meteo files...
for metfile in ../data/ei_2011_06_05_00.nc ../data/era5ml_2011_06_05_00.nc ; do

    # Set variables...
    base=$(basename $metfile .nc)
    metparam=$([ "$base" = "era5ml_2011_06_05_00" ] && echo "MET_VERT_COORD 1 MET_PRESS_LEVEL_DEF 6")
    
    # Sampling tools...
    $trac/met_lapse - data/lapse_${base}.tab $metfile $metparam
    
    $trac/met_map - data/map_${base}.tab $metfile $metparam \
		  MAP_LON0 -180 MAP_LON1 180 MAP_DLON 6 \
		  MAP_LAT0 -90 MAP_LAT1 90 MAP_DLAT 4 \
		  MAP_THETA 380 TROPO 3 \
		  OH_CHEM_BETA 0.6
    
    $trac/met_prof - data/prof_${base}.tab $metfile $metparam \
		   PROF_Z0 0 PROF_Z1 60 PROF_DZ 0.1 \
		   PROF_LON0 20 PROF_LON1 30 PROF_DLON 1 \
		   PROF_LAT0 40 PROF_LAT1 50 PROD_DLAT 1 \
		   TROPO 4
    
    $trac/met_zm - data/zm_${base}.tab $metfile $metparam \
		 ZM_Z0 0 ZM_Z1 60 ZM_DZ 1 \
		 ZM_LON0 -180 ZM_LON1 180 ZM_DLON 1 \
		 ZM_LAT0 -90 ZM_LAT1 90 ZM_DLAT 3 \
		 TROPO 5
    
    $trac/met_spec - data/spec_${base}.tab $metfile $metparam \
		   SPEC_WAVEMAX 3
    
    $trac/met_subgrid - data/subgrid_${base}.tab \
		      ../data/ei_2011_06_05_00.nc \
		      ../data/ei_2011_06_06_00.nc

done

# Create init file for sampling...
tm=$($trac/time2jsec 2011 6 5 12 0 0 0)
$trac/atm_init - data/atm_init.tab \
	       INIT_T0 "$tm" INIT_T1 "$tm" INIT_UT 86400 \
	       INIT_Z0 30 INIT_Z1 30 INIT_UZ 60 \
	       INIT_ULON 360 INIT_ULAT 180 \
	       INIT_EVENLY 1 INIT_REP 3000

# Sample meteo data...
$trac/met_sample - data/sample_ei_2011_06_05_00.tab data/atm_init.tab \
		 METBASE ../data/ei DT_MET 86400.0

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
