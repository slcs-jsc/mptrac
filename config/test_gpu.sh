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
        echo "Working on JUWELS-BOOSTER"
	ml purge
	ml NVHPC ParaStationMPI netCDF
	export LD_LIBRARY_PATH=$PWD/../libs/build/lib:$LD_LIBRARY_PATH
        # Testing
	cd src
	echo "#! /bin/bash" > check.sh
        echo "make gpu_test" >> check.sh
        sbatch --wait --deadline=now+1hour --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=48 --time=00:05:00 --account=slmet --partition=booster --gres=gpu:4 ./check.sh
        [ $? -ne 0 ] && result="FAILED" || result="OK"
        cat slurm*

        # Another Inform
	echo "job result: $result"
        echo "job finished: $(date)"
        if_found=true
    elif [[ $strIn == *"jureca"* ]] ; then
	printf "JURECA is not supported.\n" 	
    fi  
done

if [[ ! $if_found  ]]; then
    printf "None of the system matches the configuration \n"
    printf "please check that if your machine is registered on the buildbot master\n"
fi 
