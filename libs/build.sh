#! /bin/bash
# ---------------------------------------------------------------------- # 
#
# Description
#     This script is based on the original build script by Lars Hoffmann
# in order to build the external libraries used by mptrac. The user's
# manual of this script can be found by executing with -h
#
# ---------------------------------------------------------------------- # 
#
# DevLog:
#   03.08.2023: Starting the modification
#
# ---------------------------------------------------------------------- #

###### Presets ######

ifBuildAll=false
ifClean=false
ifBuildGSL=false
ifBuildTHRUST=false
ifBuildZLIB=false
ifBuildZFP=false
ifBuildZSTD=false
ifBuildHDF5=false
ifBuildNETCDF=false
ifSkip=false

numProcs=`nproc --all`

strBuildDir=${PWD}/build
strLibsDir=${PWD}

strFileGSL="gsl-2.7"
strFileHDF5="hdf5-1.12.1"
strFileNETCDF="netcdf-c-4.8.1"
strFileTHRUST="thrustsort-1.2"
strFileZFP="zfp-0.5.5"
strFileZLIB="zlib-1.2.12"
strFileZSTD="zstd-1.5.2"

###### Reminder           ######

if [[ $# -lt 1 ]]; then
    printf "$0: please use -h to open the manual\n"
    ifSkip=true
fi

###### Checking the flags ######

while getopts acgtzfshnp: flag
do
    case "${flag}" in
        a) ifBuildAll=true    ; echo "build all libraries  " ;;	
	c) ifClean=true       ; echo "clean build directory" ;;
	g) ifBuildGSL=true    ; echo "GSL is selected      " ;;
	t) ifBuildTHRUST=true ; echo "THRUST is selected   " ;;
	z) ifBuildZLIB=true   ; echo "ZLIB is selected     " ;;
	f) ifBuildZFP=true    ; echo "ZFP is selected      " ;;
	s) ifBuildZSTD=true   ; echo "ZSTD is selected     " ;;
	d) ifBuildHDF5=true   ; echo "HDF5 is selected     " ;;
	n) ifBuildNETCDF=true ; echo "NETCDF is selected   " ;;
	p) numProcs=${OPTARG}  
           numTotalProcs=`nproc --all`
           if [ $numProcs -lt $numTotalProcs ]; then
                numProcs=$numProcs
           else
		numProcs=$(($numTotalProcs*2/3))
           fi
	   echo "NEW PROCESSES as ${numProcs}" ;;
	h) printf -- "\n\n\nHelping Page: \n"
	   printf -- "--------------------------------------------------------------\n"	
	   printf -- "    Please specify which libs ( or all, -a ) you want to build\n"
	   printf -- "    you can specify within a dash, e.g.: -gtzfs to build      \n"
	   printf -- "    without NETCDF and HDF5                                   \n" 
	   printf -- "--------------------------------------------------------------\n"	
           printf -- "-a         : build all the libs, \$mptrac/libs/build will be  \n"
	   printf -- "             cleaned      \n"
	   printf -- "-c         : clean the \$mptrac/libs/build                    \n"
	   printf -- "-g         : build GSL    \n" 
	   printf -- "-t         : build THRUST     \n" 
	   printf -- "-z         : build ZLIB   \n" 
	   printf -- "-f         : build ZFP    \n" 
	   printf -- "-s         : build ZSTD   \n" 
	   printf -- "-d         : build HDF5   (prerequisite: ZLIB)  \n" 
	   printf -- "-n         : build NETCDF (prerequisite: HDF5)  \n" 
	   printf -- "-p [cores] : how many [cores] to be used by make -j [cores]   \n" 
	   printf -- "           : default is maximum number of cores               \n" 
	   printf -- "--------------------------------------------------------------\n"	
	   printf -- " Example:                                                     \n"	
	   printf -- "     ./build.sh -c         clean the build directory          \n"	
	   printf -- "     ./build.sh -a         build all the libs                 \n"	
	   printf -- "     ./build.sh -gtzfs     build without HDF5 and NetCDF      \n"	
	   printf -- "     ./build.sh -g -p 8    build the GSL with 8 processors    \n"	
	   printf -- "--------------------------------------------------------------\n"	
	   ifSkip=true
	   ;;	 
    esac
done

######

if [ $ifBuildAll = true ] || [ $ifClean = true ]; then
    printf "Clean the folder    \n"
    # Setup...
    rm -rf build
    target=$(mkdir -p build && cd build && pwd)
    printf "Making the folder with subfolders (src, bin, include, lib, man)   \n"
    mkdir -p $target/src $target/bin $target/include $target/lib $target/man/man1 
else
    if ! [ -d build ]; then
        printf "Making the folder with subfolders    \n"
        target=$(mkdir -p build && cd build && pwd)
        mkdir -p $target/src $target/bin $target/include $target/lib $target/man/man1 
    fi
fi

# GSL...
if [ $ifBuildAll = true ] || [ $ifBuildGSL = true ] ; then
    cd $strLibsDir	
    printf "Starting to compile GSL\n"
    strTarget=$strFileGSL
    cp $strTarget.tar.bz2 $strBuildDir/src && cd $strBuildDir/src && tar xvjf $strTarget.tar.bz2
    cd $strBuildDir/src/$strTarget \
        && ./configure --prefix=$strBuildDir \
        && make -j $numProcs && make check && make install && make clean \
    	|| exit
    echo -e "\n***** gsl-config *****\n"
    $strBuildDir/bin/gsl-config --libs --cflags --version
fi

# Thrust sort...
if [ $ifBuildAll = true ] || [ $ifBuildTHRUST = true ] ; then 
    cd $strLibsDir	
    printf "Starting to compile thrustsort\n"
    strTarget=$strFileTHRUST 
    cp $strTarget.tar.bz2 $strBuildDir/src && cd $strBuildDir/src && tar xvjf $strTarget.tar.bz2
    cd $strBuildDir/src/$strTarget \
        && cp -a libthrustsort_gpu.a libthrustsort_cpu.a $strBuildDir/lib/ \
    	|| exit
fi

# zlib...
if [ $ifBuildAll = true ] || [ $ifBuildZLIB = true ] ; then 
    cd $strLibsDir	
    printf "Starting to compile ZLIB\n"
    strTarget=$strFileZLIB
    cp $strTarget.tar.bz2 $strBuildDir/src && cd $strBuildDir/src && tar xvjf $strTarget.tar.bz2
    cd $strBuildDir/src/$strTarget \
        && ./configure --prefix=$strBuildDir \
        && make -j $numProcs && make check && make install && make clean \
    	|| exit
fi

# zfp...
if [ $ifBuildAll = true ] || [ $ifBuildZFP = true ] ; then 
    cd $strLibsDir	
    printf "Starting to compile ZFP\n"
    strTarget=$strFileZFP
    cp $strTarget.tar.bz2 $strBuildDir/src && cd $strBuildDir/src && tar xvjf $strTarget.tar.bz2
    cd $strBuildDir/src/$strTarget \
        && make && make test \
        && cp -a bin/* $strBuildDir/bin/ \
        && cp -a include/* $strBuildDir/include/ \
        && cp -a lib/* $strBuildDir/lib/ \
    	|| exit
fi

# zstd...
if [ $ifBuildAll = true ] || [ $ifBuildZSTD = true ] ; then 
    cd $strLibsDir	
    printf "starting to compile ZSTD \n"
    strTarget=$strFileZSTD
    cp $strTarget.tar.bz2 $strBuildDir/src && cd $strBuildDir/src && tar xvjf $strTarget.tar.bz2
    cd $strBuildDir/src/$strTarget \
        && make -j $numProcs && make check \
        && cp -a lib/libzstd* $strBuildDir/lib/ \
        && cp -a lib/*.h $strBuildDir/include/ \
        && ln -sf $strBuildDir/$strTarget/lib/libzstd.so.1.5.2 $strBuildDir/lib/libzstd.so \
        && make clean \
    	|| exit
fi

# HDF5...
if [ $ifBuildAll = true ] || [ $ifBuildHDF5 = true ] ; then 
    cd $strLibsDir	
    printf "starting to compile HDF5 with zlib\n"
    strTarget=$strFileHDF5
    cp $strTarget.tar.bz2 $strBuildDir/src && cd $strBuildDir/src && tar xvjf $strTarget.tar.bz2
    cd $strBuildDir/src/$strTarget \
        && ./configure --prefix=$strBuildDir --with-zlib=$strBuildDir --enable-hl \
        && make -j $numProcs && make check && make install && make clean \
    	|| exit
fi

# netCDF...
if [ $ifBuildAll = true ] || [ $ifBuildNETCDF = true ] ; then 
    cd $strLibsDir	
    printf "starting to compile netCDF \n"
    strTarget=$strFileNETCDF
    cp $strTarget.tar.bz2 $strBuildDir/src && cd $strBuildDir/src && tar xvjf $strTarget.tar.bz2
    cd $strBuildDir/src/$strTarget \
        && CPPFLAGS=-I$strBuildDir/include LDFLAGS=-L$strBuildDir/lib ./configure --prefix=$strBuildDir --disable-dap --disable-nczarr \
        && make -j $numProcs && make check
    make install && make clean || exit
    echo -e "\n***** nc-config *****"
    $strBuildDir/bin/nc-config --all
fi

    # # Chemistry using KPP...
    # dir=KPP
    # export KPP_HOME=$target/src/$dir
    # export KPP_FLEX_LIB_DIR=$target/src/$dir
    # cd $target/src/$dir/mptrac-chem \
    # 		&& ../bin/kpp chem.kpp && make lib && cp libkpp.a $target/lib && cp *.h $target/include

if [ $ifSkip = false ] || [ $ifBuildAll = true ] ; then 
    printf "All selected compilations are done.\n"
fi
