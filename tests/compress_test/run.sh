#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4

# Setup...
trac=../../src

# Create directories...
rm -rf data && mkdir -p data

# Downsampling and conversion from netCDF to binary...
$trac/met_conv - ../data/ei_2011_06_05_00.nc 0 data/ei_2011_06_05_00.bin 1 \
	       MET_DX 6 MET_DY 3 MET_DP 2 MET_SX 6 MET_SY 3 MET_SP 2

# Convert from binary to compressed data...
$trac/met_conv - data/ei_2011_06_05_00.bin 1 data/ei_2011_06_05_00.nc 0
$trac/met_conv - data/ei_2011_06_05_00.bin 1 data/ei_2011_06_05_00.pck 2 \
	       MET_COMP_LOGFILE data/ei_2011_06_05_00.pck.log
$trac/met_conv - data/ei_2011_06_05_00.bin 1 data/ei_2011_06_05_00.pckzstd 2 \
	       MET_PCK_ZSTD 1 \
	       MET_COMP_LOGFILE data/ei_2011_06_05_00.pckzstd.log
$trac/met_conv - data/ei_2011_06_05_00.bin 1 data/ei_2011_06_05_00.zfp 3 \
	       MET_COMP_LOGFILE data/ei_2011_06_05_00.zfp.log
$trac/met_conv - data/ei_2011_06_05_00.bin 1 data/ei_2011_06_05_00.zstd 4 \
	       MET_COMP_LOGFILE data/ei_2011_06_05_00.zstd.log
$trac/met_conv - data/ei_2011_06_05_00.bin 1 data/ei_2011_06_05_00.sz3 7 \
	       MET_COMP_LOGFILE data/ei_2011_06_05_00.sz3.log
$trac/met_conv - data/ei_2011_06_05_00.bin 1 data/ei_2011_06_05_00.zfp_scaled 3 \
	       MET_LOSSY_SCALE 1 \
	       MET_COMP_LOGFILE data/ei_2011_06_05_00.zfp_scaled.log
$trac/met_conv - data/ei_2011_06_05_00.bin 1 data/ei_2011_06_05_00.sz3_scaled 7 \
	       MET_LOSSY_SCALE 1 \
	       MET_COMP_LOGFILE data/ei_2011_06_05_00.sz3_scaled.log
$trac/met_conv - data/ei_2011_06_05_00.bin 1 data/ei_2011_06_05_00.lz4 8 \
	       MET_COMP_LOGFILE data/ei_2011_06_05_00.lz4.log

# Convert from compressed to binary data...
$trac/met_conv - data/ei_2011_06_05_00.nc 0 data/nc2bin_2011_06_05_00.bin 1
$trac/met_conv - data/ei_2011_06_05_00.pck 2 data/pck2bin_2011_06_05_00.bin 1
$trac/met_conv - data/ei_2011_06_05_00.pckzstd 2 data/pckzstd2bin_2011_06_05_00.bin 1 \
	       MET_PCK_ZSTD 1
$trac/met_conv - data/ei_2011_06_05_00.zfp 3 data/zfp2bin_2011_06_05_00.bin 1
$trac/met_conv - data/ei_2011_06_05_00.zstd 4 data/zstd2bin_2011_06_05_00.bin 1
$trac/met_conv - data/ei_2011_06_05_00.sz3 7 data/sz32bin_2011_06_05_00.bin 1
$trac/met_conv - data/ei_2011_06_05_00.zfp_scaled 3 data/zfp_scaled2bin_2011_06_05_00.bin 1 \
	       MET_LOSSY_SCALE 1
$trac/met_conv - data/ei_2011_06_05_00.sz3_scaled 7 data/sz3_scaled2bin_2011_06_05_00.bin 1 \
	       MET_LOSSY_SCALE 1
$trac/met_conv - data/ei_2011_06_05_00.lz4 8 data/lz42bin_2011_06_05_00.bin 1

# Compare files...
echo -e "\nCompare results..."
error=0
for f in data.ref/*.bin ; do
    [ -s data/"$(basename "$f")" ] || continue
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
