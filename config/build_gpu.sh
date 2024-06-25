#!/bin/bash

# ==============================================================
# Here is the exmample of a build script for WRF. 
# It should be noted that the configuration files are tested 
# which are based on the the WRF configuration executor.  
# ==============================================================

strHPC=`uname -a`
if_found=false
for strIn in $strHPC ; do
    ## Check which machine
    if [[ $strIn  == *"juwels"* ]] ; then
	# Compiling
        echo "Working on compiling the MPTRAC with GPU"
	export LD_LIBRARY_PATH=$PWD/libs/build/lib:$LD_LIBRARY_PATH
	ml NVHPC ParaStationMPI netCDF  
	cd src
        make clean 
	make COMPILER=nvc GPU=1 STATIC=0 INFO=1 ZSTD=1 ZFP=1 DEFINES="-DNP=20000" || exit
        echo "Done compiling the MPTRAC with GPU"
        if_found=true 
    elif [[ $strIn  == *"jureca"* ]] ; then
	printf "JURECA is not supported.\n" 
	printf "Currently supporting JUWELS-BOOST only\n"
        if_found=true 
    fi 
done

if [[ ! $if_found  ]]; then
    printf "None of the system matches the configuration \n"
    printf "please check that if your machine is registered on the buildbot master\n"
fi 
