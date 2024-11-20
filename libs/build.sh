#! /bin/bash

###### Presets ######

ifBuildAll=false
ifBuildMandatory=false
ifClean=false
ifBuildGSL=false
ifBuildTHRUST=false
ifBuildCURL=false
ifBuildZLIB=false
ifBuildSZIP=false
ifBuildZFP=false
ifBuildZSTD=false
ifBuildHDF5=false
ifBuildNETCDF=false
ifBuildKPP=false
ifSkip=false

numProcs=`nproc --all`

strBuildDir=${PWD}/build
strLibsDir=${PWD}

strFileGSL="gsl-2.7.1"
strFileHDF5="hdf5-1.14.4-3"
strFileNETCDF="netcdf-c-4.9.2"
strFileTHRUST="thrustsort-1.2"
strFileZFP="zfp-1.0.1"
strFileCURL="curl-8.10.1"
strFileZLIB="zlib-1.3.1"
strFileSZIP="szip-2.1.1"
strFileZSTD="zstd-1.5.5"
strFileKPP="KPP"

export LD_LIBRARY_PATH=$strBuildDir/lib:$LD_LIBRARY_PATH

###### Reminder           ######

if [[ $# -lt 1 ]]; then
    printf "$0: please use -h to open the manual\n"
    ifSkip=true
fi

###### Checking the flags ######

while getopts amcgtuzifshnkp: flag
do
    case "${flag}" in
        a) ifBuildAll=true        ; echo "build all libraries  " ;;	
        m) ifBuildMandatory=true  ; echo "build mandatory libraries  " ;;	
	c) ifClean=true           ; echo "clean build directory" ;;
	g) ifBuildGSL=true    ; echo "GSL is selected      " ;;
	t) ifBuildTHRUST=true ; echo "THRUST is selected   " ;;
	u) ifBuildCURL=true   ; echo "CURL is selected     " ;;
	z) ifBuildZLIB=true   ; echo "ZLIB is selected     " ;;
	i) ifBuildSZIP=true   ; echo "SZIP is selected     " ;;
	f) ifBuildZFP=true    ; echo "ZFP is selected      " ;;
	s) ifBuildZSTD=true   ; echo "ZSTD is selected     " ;;
	d) ifBuildHDF5=true   ; echo "HDF5 is selected     " ;;
	n) ifBuildNETCDF=true ; echo "NETCDF is selected   " ;;
	k) ifBuildKPP=true    ; echo "KPP is selected      " ;;
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
	   printf -- "    Please specify which libs (or all, -a) you want to build  \n"
	   printf -- "    you can specify within a dash, e.g.: -gtzfs to build      \n"
	   printf -- "    without NETCDF and HDF5                                   \n" 
	   printf -- "--------------------------------------------------------------\n"	
           printf -- "-a         : build all the libs, \$mptrac/libs/build will be  \n"
	   printf -- "             cleaned      \n"
	   printf -- "-m         : build mandatory libs: GSL, netCDF, HDF5, CURL, ZLIB and SZIP\n"
	   printf -- "-c         : clean the \$mptrac/libs/build                    \n"
	   printf -- "-g         : build GSL    \n" 
	   printf -- "-t         : build THRUST     \n"
	   printf -- "-u         : build CURL   \n" 
	   printf -- "-z         : build ZLIB   \n" 
	   printf -- "-f         : build ZFP    \n" 
	   printf -- "-s         : build ZSTD   \n" 
	   printf -- "-d         : build HDF5   (prerequisites: CURL, ZLIB, SZIP)  \n" 
	   printf -- "-n         : build NETCDF (prerequisite: HDF5)  \n" 
	   printf -- "-k         : build KPP    \n" 
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

# Clean directories..
if [ $ifBuildAll = true ] || [ $ifClean = true ]; then
    printf "Clean the folder    \n"
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

# Mandatory ... 
if [ $ifBuildMandatory = true ]; then
    ifBuildGSL=true
    ifBuildCURL=true
    ifBuildZLIB=true
    ifBuildSZIP=true
    ifBuildHDF5=true
    ifBuildNETCDF=true
fi

# GSL...
if [ $ifBuildAll = true ] || [ $ifBuildGSL = true ] ; then
    cd $strLibsDir	
    printf "Starting to compile GSL...\n"
    strTarget=$strFileGSL
    cp $strTarget.tar.bz2 $strBuildDir/src && cd $strBuildDir/src && tar xvjf $strTarget.tar.bz2
    cd $strBuildDir/src/$strTarget \
        && ./configure --prefix=$strBuildDir \
        && make -j $numProcs && make check && make install && make clean \
    	    || exit
fi

# Thrust sort...
if [ $ifBuildAll = true ] || [ $ifBuildTHRUST = true ] ; then 
    cd $strLibsDir	
    printf "Starting to compile thrustsort...\n"
    strTarget=$strFileTHRUST 
    cp $strTarget.tar.bz2 $strBuildDir/src && cd $strBuildDir/src && tar xvjf $strTarget.tar.bz2
    cd $strBuildDir/src/$strTarget \
        && cp -a libthrustsort_gpu.a libthrustsort_cpu.a $strBuildDir/lib/ \
    	    || exit
fi

# curl...
if [ $ifBuildAll = true ] || [ $ifBuildCURL = true ] ; then 
    cd $strLibsDir	
    printf "Starting to compile CURL...\n"
    export CFLAGS="-fPIC"
    strTarget=$strFileCURL
    cp $strTarget.tar.bz2 $strBuildDir/src && cd $strBuildDir/src && tar xvjf $strTarget.tar.bz2
    cd $strBuildDir/src/$strTarget \
        && ./configure --prefix=$strBuildDir --without-ssl --without-libpsl \
        && make -j $numProcs && make check
    make install && make clean \
    	    || exit
fi

# zlib...
if [ $ifBuildAll = true ] || [ $ifBuildZLIB = true ] ; then 
    cd $strLibsDir	
    printf "Starting to compile ZLIB...\n"
    export CFLAGS="-fPIC"
    strTarget=$strFileZLIB
    cp $strTarget.tar.bz2 $strBuildDir/src && cd $strBuildDir/src && tar xvjf $strTarget.tar.bz2
    cd $strBuildDir/src/$strTarget \
        && ./configure --prefix=$strBuildDir \
        && make -j $numProcs && make check && make install && make clean \
    	    || exit
fi

# szip...
if [ $ifBuildAll = true ] || [ $ifBuildSZIP = true ] ; then 
    cd $strLibsDir	
    printf "Starting to compile SZIP...\n"
    export CFLAGS="-fPIC"
    strTarget=$strFileSZIP
    cp $strTarget.tar.bz2 $strBuildDir/src && cd $strBuildDir/src && tar xvjf $strTarget.tar.bz2
    cd $strBuildDir/src/$strTarget \
        && ./configure --prefix=$strBuildDir \
        && make -j $numProcs && make check && make install && make clean \
    	    || exit
fi

# zfp...
if [ $ifBuildAll = true ] || [ $ifBuildZFP = true ] ; then 
    cd $strLibsDir	
    printf "Starting to compile ZFP...\n"
    strTarget=$strFileZFP
    cp $strTarget.tar.bz2 $strBuildDir/src && cd $strBuildDir/src && tar xvjf $strTarget.tar.bz2
    cd $strBuildDir/src/$strTarget \
        && make BUILD_SHARED_LIBS=0 && make test \
        && cp -a bin/* $strBuildDir/bin/ \
        && cp -a include/* $strBuildDir/include/ \
        && cp -a lib/* $strBuildDir/lib/ \
	&& make BUILD_SHARED_LIBS=1 \
	&& cp -a lib/* $strBuildDir/lib/ \
    	    || exit
fi

# zstd...
if [ $ifBuildAll = true ] || [ $ifBuildZSTD = true ] ; then 
    cd $strLibsDir	
    printf "starting to compile ZSTD...\n"
    strTarget=$strFileZSTD
    cp $strTarget.tar.bz2 $strBuildDir/src && cd $strBuildDir/src && tar xvjf $strTarget.tar.bz2
    cd $strBuildDir/src/$strTarget \
        && make -j $numProcs && make check \
        && cp -a lib/libzstd* $strBuildDir/lib/ \
        && cp -a lib/*.h $strBuildDir/include/ \
        && ln -sf $strBuildDir/$strTarget/lib/libzstd.so.1.5.5 $strBuildDir/lib/libzstd.so \
        && make clean \
    	    || exit
fi

# HDF5...
if [ $ifBuildAll = true ] || [ $ifBuildHDF5 = true ] ; then 
    cd $strLibsDir	
    printf "starting to compile HDF5...\n"
    strTarget=$strFileHDF5
    cp $strTarget.tar.bz2 $strBuildDir/src && cd $strBuildDir/src && tar xvjf $strTarget.tar.bz2
    cd $strBuildDir/src/$strTarget \
        && ./configure --prefix=$strBuildDir --with-zlib=$strBuildDir --with-szlib=$strBuildDir --enable-hl --disable-fortran \
        && make -j $numProcs && make check && make install && make clean \
	    || exit
fi

# netCDF...
if [ $ifBuildAll = true ] || [ $ifBuildNETCDF = true ] ; then 
    cd $strLibsDir	
    printf "starting to compile netCDF...\n"
    strTarget=$strFileNETCDF
    cp $strTarget.tar.bz2 $strBuildDir/src && cd $strBuildDir/src && tar xvjf $strTarget.tar.bz2
    cd $strBuildDir/src/$strTarget \
        && CPPFLAGS=-I$strBuildDir/include LDFLAGS=-L$strBuildDir/lib ./configure --prefix=$strBuildDir --with-plugin-dir=$strBuildDir/share/netcdf-plugins --disable-dap --disable-byterange --disable-nczarr --disable-libxml2 \
        && make -j $numProcs && make check
    make install && make clean \
	    || exit
fi

# KPP chemistry...
if [ $ifBuildAll = true ] || [ $ifBuildKPP = true ] ; then 
    cd $strLibsDir	
    printf "starting to compile KPP...\n"
    strTarget=$strFileKPP
    cp $strTarget.tar.bz2 $strBuildDir/src && cd $strBuildDir/src && tar xvjf $strTarget.tar.bz2
    export KPP_HOME=$strBuildDir/src/$strTarget
    export KPP_FLEX_LIB_DIR=$strBuildDir/src/$strTarget
    cd $KPP_HOME/src && make clean && make && cd ${strLibsDir}/../projects/kpp && ./build_KPP.sh chem 0 \
	    || exit		
fi

# Finish...
if [ $ifSkip = false ] || [ $ifBuildAll = true ] ; then 
    printf "All selected compilations are done.\n"
    echo -e "\n***** gcc *****\n"
    gcc --version
    echo -e "\n***** gsl-config *****\n"
    $strBuildDir/bin/gsl-config --libs --cflags --version
    echo -e "\n***** nc-config *****"
    $strBuildDir/bin/nc-config --all    
fi
