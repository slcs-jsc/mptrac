#! /bin/bash

###### Presets ######

# Flags...
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
ifBuildECCODES=false
ifBuildSZ3=false
ifBuildMGARD=false

# System...
numProcs=$(nproc --all)

# Directories...
BuildDir=${PWD}/build
LibsDir=${PWD}

export LD_LIBRARY_PATH=$BuildDir/lib:$LD_LIBRARY_PATH

###### Functions          ######

print_help() {
    cat << EOF

MPTRAC Library Build Script
---------------------------

This script builds the required libraries and dependencies for the
MPTRAC Lagrangian particle dispersion model.

Usage:
  ./build.sh [options]

Options:
  -a           Clean build directory and build ALL libraries
  -m           Build MANDATORY libraries (GSL, netCDF, HDF5)
  -c           Clean the build directory only
  -g           Build GSL
  -t           Build THRUST
  -u           Build CURL
  -z           Build ZLIB
  -i           Build SZIP
  -f           Build ZFP
  -s           Build ZSTD
  -j           Build SZ3
  -l           Build MGARD
  -d           Build HDF5 (requires CURL, ZLIB, SZIP)
  -n           Build NETCDF (requires HDF5)
  -k           Build KPP chemistry
  -e           Build ECCODES
  -p [cores]   Number of CPU cores to use for compilation with make -j [cores]
               Default: all available cores
  -h           Show this help message and exit

Examples:
  ./build.sh -c
      Clean the build directory.

  ./build.sh -a
      Build all libraries.

  ./build.sh -gtzfs
      Build GSL, THRUST, ZLIB, ZFP, and SZIP (without HDF5 and NETCDF).

  ./build.sh -g -p 8
      Build GSL using 8 processor cores.

EOF
    exit 0
}

###### Check arguments    ######

[[ $# -lt 1 ]] && print_help
while getopts amcgtuzijlfshnkep: flag
do
    case "${flag}" in
        a) ifBuildAll=true       ; echo "build all libraries      " ;;
        m) ifBuildMandatory=true ; echo "build mandatory libraries" ;;
	c) ifClean=true          ; echo "clean build directory    " ;;
	g) ifBuildGSL=true       ; echo "GSL is selected          " ;;
	t) ifBuildTHRUST=true    ; echo "THRUST is selected       " ;;
	u) ifBuildCURL=true      ; echo "CURL is selected         " ;;
	z) ifBuildZLIB=true      ; echo "ZLIB is selected         " ;;
	i) ifBuildSZIP=true      ; echo "SZIP is selected         " ;;
	f) ifBuildZFP=true       ; echo "ZFP is selected          " ;;
	s) ifBuildZSTD=true      ; echo "ZSTD is selected         " ;;
	j) ifBuildSZ3=true       ; echo "SZ3 is selected          " ;;
	l) ifBuildMGARD=true     ; echo "MGARD is selected        " ;;
	d) ifBuildHDF5=true      ; echo "HDF5 is selected         " ;;
	n) ifBuildNETCDF=true    ; echo "NETCDF is selected       " ;;
	k) ifBuildKPP=true       ; echo "KPP is selected          " ;;
	e) ifBuildECCODES=true   ; echo "ECCODES is selected      " ;;
	p) numProcs=${OPTARG} ;;
	h) print_help ;;
    esac
done

# Clean build directory...
if [ "$ifBuildAll" = true ] || [ "$ifClean" = true ] ; then
    printf "Cleaning build folder...\n"
    rm -rf build
fi

# Create build directory...
printf "Creating build folder with subfolders...\n"
mkdir -p build/src build/bin build/include build/lib build/man/man1 || exit

# Select mandatory libraries... 
if [ $ifBuildMandatory = true ] ; then
    ifBuildGSL=true
    ifBuildCURL=true
    ifBuildZLIB=true
    ifBuildSZIP=true
    ifBuildHDF5=true
    ifBuildNETCDF=true
fi

# Build GSL...
if [ $ifBuildAll = true ] || [ $ifBuildGSL = true ] ; then
    Target="gsl-2.7.1"
    printf "Starting to compile $Target...\n"
    cd $LibsDir	&& cp $Target.tar.bz2 $BuildDir/src && cd $BuildDir/src && tar xvjf $Target.tar.bz2 || exit
    cd $BuildDir/src/$Target \
        && ./configure --prefix=$BuildDir \
        && make -j $numProcs && make check && make install && make clean \
    	    || exit
fi

# Build Thrust sort...
if [ $ifBuildAll = true ] || [ $ifBuildTHRUST = true ] ; then 
    Target="thrustsort-1.2"
    printf "Starting to compile $Target...\n"
    cd $LibsDir && cp $Target.tar.bz2 $BuildDir/src && cd $BuildDir/src && tar xvjf $Target.tar.bz2 || exit
    cd $BuildDir/src/$Target \
        && cp -a libthrustsort_gpu.a libthrustsort_cpu.a $BuildDir/lib/ \
    	    || exit
fi

# Build curl...
if [ $ifBuildAll = true ] || [ $ifBuildCURL = true ] ; then 
    Target="curl-8.10.1"
    printf "Starting to compile $Target...\n"
    cd $LibsDir && cp $Target.tar.bz2 $BuildDir/src && cd $BuildDir/src && tar xvjf $Target.tar.bz2 || exit
    cd $BuildDir/src/$Target \
        && CFLAGS="-fPIC" ./configure --prefix=$BuildDir --without-ssl --without-libpsl \
        && make -j $numProcs && make check || true && make install && make clean \
    	    || exit
fi

# Build zlib...
if [ $ifBuildAll = true ] || [ $ifBuildZLIB = true ] ; then 
    Target="zlib-1.3.1"
    printf "Starting to compile $Target...\n"
    cd $LibsDir && cp $Target.tar.bz2 $BuildDir/src && cd $BuildDir/src && tar xvjf $Target.tar.bz2 || exit
    cd $BuildDir/src/$Target \
        && CFLAGS="-fPIC" ./configure --prefix=$BuildDir \
        && make -j $numProcs && make check && make install && make clean \
    	    || exit
fi

# Build szip...
if [ $ifBuildAll = true ] || [ $ifBuildSZIP = true ] ; then 
    Target="szip-2.1.1"
    printf "Starting to compile $Target...\n"
    cd $LibsDir && cp $Target.tar.bz2 $BuildDir/src && cd $BuildDir/src && tar xvjf $Target.tar.bz2 || exit
    cd $BuildDir/src/$Target \
        && CFLAGS="-fPIC" ./configure --prefix=$BuildDir \
        && make -j $numProcs && make check && make install && make clean \
    	    || exit
fi

# Build zfp...
if [ $ifBuildAll = true ] || [ $ifBuildZFP = true ] ; then 
    Target="zfp-1.0.1"
    printf "Starting to compile $Target...\n"
    cd $LibsDir && cp $Target.tar.bz2 $BuildDir/src && cd $BuildDir/src && tar xvjf $Target.tar.bz2 || exit
    cd $BuildDir/src/$Target \
        && make BUILD_SHARED_LIBS=0 && make test \
        && cp -a bin/* $BuildDir/bin/ \
        && cp -a include/* $BuildDir/include/ \
        && cp -a lib/* $BuildDir/lib/ \
	&& make BUILD_SHARED_LIBS=1 \
	&& cp -a lib/* $BuildDir/lib/ \
    	    || exit
fi

# Build zstd...
if [ $ifBuildAll = true ] || [ $ifBuildZSTD = true ] ; then 
    Target="zstd-1.5.5"
    printf "Starting to compile $Target...\n"
    cd $LibsDir && cp $Target.tar.bz2 $BuildDir/src && cd $BuildDir/src && tar xvjf $Target.tar.bz2 || exit
    cd $BuildDir/src/$Target \
        && make CFLAGS=-fPIC -j $numProcs && make check \
        && cp -a lib/libzstd* $BuildDir/lib/ \
        && cp -a lib/*.h $BuildDir/include/ \
        && make clean \
    	    || exit
fi

# Build HDF5...
if [ $ifBuildAll = true ] || [ $ifBuildHDF5 = true ] ; then 
    Target="hdf5-1.14.4-3"
    printf "Starting to compile $Target...\n"
    cd $LibsDir && cp $Target.tar.bz2 $BuildDir/src && cd $BuildDir/src && tar xvjf $Target.tar.bz2 || exit
    cd $BuildDir/src/$Target \
        && ./configure --prefix=$BuildDir --with-zlib=$BuildDir --with-szlib=$BuildDir --enable-hl --disable-fortran \
        && make -j $numProcs && make check && make install && make clean \
	    || exit
fi

# Build netCDF...
if [ $ifBuildAll = true ] || [ $ifBuildNETCDF = true ] ; then 
    Target="netcdf-c-4.9.2"
    printf "Starting to compile $Target...\n"
    cd $LibsDir && cp $Target.tar.bz2 $BuildDir/src && cd $BuildDir/src && tar xvjf $Target.tar.bz2 || exit
    cd $BuildDir/src/$Target \
        && CPPFLAGS=-I$BuildDir/include LDFLAGS=-L$BuildDir/lib ./configure --prefix=$BuildDir --with-plugin-dir=$BuildDir/share/netcdf-plugins --disable-dap --disable-byterange --disable-nczarr --disable-libxml2 \
        && make -j $numProcs && make check || true && make install && make clean \
	    || exit
fi

# Build KPP...
if [ $ifBuildAll = true ] || [ $ifBuildKPP = true ] ; then 
    Target="KPP"
    printf "Starting to compile $Target...\n"
    cd $LibsDir && cp $Target.tar.bz2 $BuildDir/src && cd $BuildDir/src && tar xvjf $Target.tar.bz2 || exit
    export KPP_HOME=$BuildDir/src/$Target
    export KPP_FLEX_LIB_DIR=$BuildDir/src/$Target
    cd $KPP_HOME/src && make clean && make && cd ${LibsDir}/../projects/kpp && ./build_KPP.sh chem 0 \
	    || exit		
fi

# Build ecCodes...
if [ $ifBuildAll = true ] || [ $ifBuildECCODES = true ] ; then
    Target="eccodes-2.38.3"
    printf "Starting to compile $Target...\n"
    cd $LibsDir && cp $Target.tar.bz2 $BuildDir/src && cd $BuildDir/src && tar xvjf $Target.tar.bz2 || exit
    mkdir -p $BuildDir/src/$Target/build && cd $BuildDir/src/$Target/build || exit
    cmake -DCMAKE_INSTALL_PREFIX=$BuildDir -DBUILD_SHARED_LIBS=BOTH -DENABLE_AEC=OFF -DENABLE_JPEG=OFF -DENABLE_FORTRAN=OFF -DENABLE_NETCDF=OFF .. && make -j $numProcs && ctest -j $numProcs && make install \
            || exit
fi

# Build SZ3...
if [ $ifBuildAll = true ] || [ $ifBuildSZ3 = true ] ; then
    Target="SZ3-3.2.1"
    printf "Starting to compile $Target...\n"
    cd $LibsDir && cp $Target.tar.bz2 $BuildDir/src && cd $BuildDir/src && tar xvjf $Target.tar.bz2 || exit
    mkdir -p $BuildDir/src/$Target/build_shared && cd $BuildDir/src/$Target/build_shared || exit
    cmake -DCMAKE_INSTALL_PREFIX=$BuildDir -DBUILD_SHARED_LIBS=ON .. && make -j $numProcs && make install || exit
    mkdir -p $BuildDir/src/$Target/build_static && cd $BuildDir/src/$Target/build_static || exit
    cmake -DCMAKE_INSTALL_PREFIX=$BuildDir -DBUILD_SHARED_LIBS=OFF .. && make -j $numProcs && make install || exit
fi

# Build MGARD...
if [ $ifBuildAll = true ] || [ $ifBuildMGARD = true ] ; then
    Target="MGARD-1.5.2"
    printf "Starting to compile $Target...\n"
    cd $LibsDir && cp $Target.tar.bz2 $BuildDir/src && cd $BuildDir/src && tar xvjf $Target.tar.bz2 || exit
    for type in shared static; do
        builddir=$BuildDir/src/$Target/build_$type
        cmake -S $BuildDir/src/$Target -B $builddir -DCMAKE_INSTALL_PREFIX=$BuildDir -DMGARD_ENABLE_OPENMP=ON \
              -DBUILD_SHARED_LIBS=$([ "$type" = shared ] && echo ON || echo OFF) || exit
        cmake --build $builddir -j $numProcs && cmake --install $builddir || exit
    done
fi

# Finish...
printf "All selected compilations are done.\n"
echo -e "\n***** gcc *****\n"
gcc --version
echo -e "\n***** gsl-config *****\n"
$BuildDir/bin/gsl-config --libs --cflags --version
echo -e "\n***** nc-config *****"
$BuildDir/bin/nc-config --all
