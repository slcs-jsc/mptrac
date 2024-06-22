#! /bin/bash

# Setup...
trac=../../src
metdir=/p/data1/slmet/datapub/hoffmann1/mptrac/projects/benchmarking/meteo

# Create directories...
rm -rf $metdir && mkdir -p $metdir/erai_nc $metdir/erai_pck $metdir/era5_nc $metdir/era5_pck || exit

# Copy netCDF files...
cp /p/data1/slmet/model_data/mptrac/input/erai/2017/ei_2017_01_01_*.nc \
   /p/data1/slmet/model_data/mptrac/input/erai/2017/ei_2017_01_02_00.nc \
   $metdir/erai_nc/
cp /p/data1/slmet/model_data/mptrac/input/era5/2017/01/era5_2017_01_01_*.nc \
   /p/data1/slmet/model_data/mptrac/input/era5/2017/01/era5_2017_01_02_00.nc \
   $metdir/era5_nc/

# Packing...
for f in $(ls $metdir/erai_nc/*nc) ; do
    $trac/met_conv - "$f" 0 $metdir/erai_pck/$(basename "$f" .nc).pck 2
done
for f in $(ls $metdir/era5_nc/*nc) ; do
    $trac/met_conv - "$f" 0 $metdir/era5_pck/$(basename "$f" .nc).pck 2
done
