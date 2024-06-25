#!/bin/bash

strFolderIn=$1
str_in=`uname -a`
if_found=.false.
for i in $str_in ; do 
    if [[ $i  == *"juwels"* ]] ; then
        printf "we are on juwels\n"
        printf "`uname -a`\n"
	printf "PWD: `pwd`\n"
	printf "Loading Modules\n"
        module load Intel ParaStationMPI netCDF GSL
	printf "Compiling MPTRAC using modules on JUWELS ... \n "
	cd ./src
        make STATIC=0
        make check	
        if_found=.true. 

    elif [[ $i == *"jureca"* ]] ; then
        printf "we are on jureca\n"
        printf "`uname -a`\n"
	printf "PWD: `pwd`\n"
	printf "Loading Modules\n"
        module load Intel ParaStationMPI netCDF GSL
	printf "Compiling MPTRAC using modules on JURECA ... \n "
	cd ./src
        make STATIC=0
        make check	
        if_found=.true. 
    else
	printf "None of the machine matches\n" 
    fi 
done

if [[ ! $if_found  ]]; then 
    printf "None of the system matches the configuration \n"
    printf "please check that if your machine is registered on the buildbot master\n"	
fi
