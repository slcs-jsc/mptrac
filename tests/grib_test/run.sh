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

$trac/met_zm - ./data/testzm.tab\
		./gribdata/gb_2011060600_small_XX.grb\
		MET_TYPE 6 \
		MET_PRESS_LEVEL_DEF 6 

$trac/met_map - ./data/testmap.tab\
		./gribdata/gb_2011060600_small_XX.grb\
		MET_TYPE 6 \
		MET_PRESS_LEVEL_DEF 6\
		MAP_LON0 5.1\
		MAP_LON1 15\ 
		MAP_LAT0 45\
		MAP_LAT1 54.9\


# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/$(basename "$f") "$f"
    [ $? -ne 0 ] && error=1
done
exit $error
