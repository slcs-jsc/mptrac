#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../../libs/build/lib:../../../cmultiscale/cmultiscale/builddir/
export OMP_NUM_THREADS=4

# Setup...
trac=../../src

# Create directories...
rm -rf data && mkdir -p data

# Downsampling and conversion from netCDF to cmultiscale...
$trac/met_conv - ../data/ei_2011_06_05_00.nc 0 data/ei_2011_06_05_00.cms 5 MET_DP 15

# Extract map data...
$trac/met_map - data/map.tab data/ei_2011_06_05_00.cms MET_TYPE 5 MAP_DLON 3 MAP_DLAT 2

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
