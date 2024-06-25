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
    if [[ $strIn == *"juwels"* ]] ; then
	printf "We are on JUWELS ... \n"
        printf "compiling something.. \n "
        ### build script starts here .... 
	module load Intel ParaStationMPI netCDF GSL
        cd ./libs
	sh build.sh -ctsgf
        make check
        ### build script ends here .... 

        if_found=.true.

    elif [[ $strIn == *"jureca"* ]] ; then
	printf "JURECA is not supported.\n" 
    fi 
done

if [[ ! $if_found  ]]; then
    printf "None of the system matches the configuration \n"
    printf "please check that if your machine is registered on the buildbot master\n"
fi 
