#! /bin/bash

# --------------------------------------------------
# Presets...
# --------------------------------------------------

# Build flags...
declare -A build
for lib in all mandatory clean parallel gsl thrust zlib szip zfp zstd hdf5 netcdf kpp eccodes sz3 ; do
    build[$lib]=false
done

# Directories...
build_dir="${PWD}/build"
libs_dir="${PWD}"

# System...
export LD_LIBRARY_PATH="$build_dir"/lib:$LD_LIBRARY_PATH
nprocs=$(nproc --all)

# Ensure PIC code...
export CFLAGS="-fPIC"
export CXXFLAGS="-fPIC"

# --------------------------------------------------
# Functions...
# --------------------------------------------------

# Print help message and exit...
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
  -m           Build MANDATORY libraries (GSL, HDF5, NETCDF)
  -c           Clean the build directory only
  -g           Build GSL
  -t           Build THRUST
  -z           Build ZLIB
  -i           Build SZIP
  -f           Build ZFP
  -s           Build ZSTD
  -j           Build SZ3
  -d           Build HDF5 (requires ZLIB and SZIP)
  -n           Build NETCDF (requires HDF5)
  -l           Enable parallel build of HDF5 and NETCDF
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

# Unpack tar file...
unpack() {
    local target=$1
    printf "Starting to compile $target...\n"
    cp "$libs_dir/$target.tar.bz2" "$build_dir/src" \
	&& cd "$build_dir/src" \
	&& tar xvjf "$target.tar.bz2" \
	&& cd "$target" \
	    || exit
}

# --------------------------------------------------
# Compile...
# --------------------------------------------------

[[ $# -lt 1 ]] && print_help
while getopts amcgtzijfshnlkep: flag ; do
    case "${flag}" in
        a) build[all]=true       ; echo "build all libraries      " ;;
        m) build[mandatory]=true ; echo "build mandatory libraries" ;;
	c) build[clean]=true     ; echo "clean build directory    " ;;
	g) build[gsl]=true       ; echo "GSL is selected          " ;;
	t) build[thrust]=true    ; echo "THRUST is selected       " ;;
	z) build[zlib]=true      ; echo "ZLIB is selected         " ;;
	i) build[szip]=true      ; echo "SZIP is selected         " ;;
	f) build[zfp]=true       ; echo "ZFP is selected          " ;;
	s) build[zstd]=true      ; echo "ZSTD is selected         " ;;
	j) build[sz3]=true       ; echo "SZ3 is selected          " ;;
	d) build[hdf5]=true      ; echo "HDF5 is selected         " ;;
	n) build[netcdf]=true    ; echo "NETCDF is selected       " ;;
	l) build[parallel]=true  ; echo "Enable parallel NETCDF   " ;;
	k) build[kpp]=true       ; echo "KPP is selected          " ;;
	e) build[eccodes]=true   ; echo "ECCODES is selected      " ;;
	p) nprocs=${OPTARG} ;;
	h) print_help ;;
    esac
done

# Clean build directory...
if ${build[all]} || ${build[clean]} ; then
    printf "Cleaning build folder...\n"
    rm -rf build
fi

# Create build directory...
printf "Creating build folder with subfolders...\n"
mkdir -p build/src build/bin build/include build/lib build/man/man1 || exit

# Select mandatory libraries... 
if ${build[mandatory]} ; then
    build[gsl]=true
    build[zlib]=true
    build[szip]=true
    build[hdf5]=true
    build[netcdf]=true
fi

# Check availability of MPI...
if ${build[parallel]} ; then
    command -v mpicc >/dev/null 2>&1 || { echo "Error: mpicc not found, cannot do parallel build!" ; exit 1 ; }
fi

# Build GSL...
if ${build[all]} || ${build[gsl]} ; then
    target="gsl-2.7.1"
    unpack "$target"
    ./configure --prefix="$build_dir" && make -j $nprocs && make check && make install && make clean || exit
fi

# Build Thrust sort...
if ${build[all]} || ${build[thrust]} ; then 
    target="thrustsort-1.3"
    unpack "$target"
    cp -a libthrustsort_gpu.a libthrustsort_cpu.a "$build_dir"/lib/ || exit
    cp -a thrustsort.h "$build_dir"/include/ || exit
fi

# Build zlib...
if ${build[all]} || ${build[zlib]} ; then 
    target="zlib-1.3.1"
    unpack "$target"
    ./configure --prefix="$build_dir" && make -j $nprocs && make check && make install && make clean || exit
fi

# Build szip...
if ${build[all]} || ${build[szip]} ; then 
    target="szip-2.1.1"
    unpack "$target"
    ./configure --prefix="$build_dir" && make -j $nprocs && make check && make install && make clean || exit
fi

# Build zfp...
if ${build[all]} || ${build[zfp]} ; then 
    target="zfp-1.0.1"
    unpack "$target"
    make BUILD_SHARED_LIBS=0 && make test \
        && cp -a bin/* "$build_dir"/bin/ \
        && cp -a include/* "$build_dir"/include/ \
        && cp -a lib/* "$build_dir"/lib/ || exit
    make BUILD_SHARED_LIBS=1 \
	&& cp -a lib/* "$build_dir"/lib/ || exit
fi

# Build zstd...
if ${build[all]} || ${build[zstd]} ; then 
    target="zstd-1.5.5"
    unpack "$target"
    make -j $nprocs && make check \
        && cp -a lib/libzstd* "$build_dir"/lib/ \
        && cp -a lib/*.h "$build_dir"/include/ \
        && make clean || exit
fi

# Build HDF5...
if ${build[all]} || ${build[hdf5]} ; then 
    target="hdf5-1.14.4-3"
    unpack "$target"
    cfg_opt="--prefix=$build_dir --with-zlib=$build_dir --with-szlib=$build_dir --enable-hl --disable-fortran --enable-shared --enable-static --disable-dap --disable-logging --without-curl"
    if ${build[parallel]} ; then
        CC=mpicc CXX=mpicxx ./configure $cfg_opt --enable-parallel || exit
    else
	./configure $cfg_opt || exit
    fi
    make -j $nprocs && make check && make install && make clean || exit
fi

# Build netCDF...
if ${build[all]} || ${build[netcdf]} ; then 
    target="netcdf-c-4.9.2"
    unpack "$target"
    cfg_opt="--prefix=$build_dir --with-hdf5=$build_dir --with-plugin-dir=$build_dir/share/netcdf-plugins --disable-dap --disable-byterange --disable-nczarr --disable-libxml2 --with-zlib=$build_dir"
    if ${build[parallel]}; then
        CPPFLAGS="-I$build_dir/include" LDFLAGS="-L$build_dir/lib" CC=mpicc ./configure $cfg_opt --enable-netcdf-4 --enable-parallel4 || exit
    else
        CPPFLAGS="-I$build_dir/include" LDFLAGS="-L$build_dir/lib" ./configure $cfg_opt || exit
    fi
    make -j $nprocs && make check || true && make install && make clean || exit
fi

# Build KPP...
if ${build[all]} || ${build[kpp]} ; then 
    target="KPP"
    unpack "$target"
    cd - 
    export KPP_HOME="$build_dir"/src/"$target"
    export KPP_FLEX_LIB_DIR="$build_dir"/src/"$target"
    cd $KPP_HOME/src && make clean && make && cd ${libs_dir}/../projects/kpp && ./build_KPP.sh chem 0 || exit		
fi

# Build ecCodes...
if ${build[all]} || ${build[eccodes]} ; then
    target="eccodes-2.38.3"
    unpack "$target"
    cd -
    mkdir -p "$build_dir"/src/"$target"/build && cd "$build_dir"/src/"$target"/build || exit
    cmake -DCMAKE_INSTALL_PREFIX="$build_dir" -DBUILD_SHARED_LIBS=BOTH -DENABLE_AEC=OFF -DENABLE_JPEG=OFF -DENABLE_FORTRAN=OFF -DENABLE_NETCDF=OFF .. \
	&& make -j $nprocs && ctest -j $nprocs && make install || exit
fi

# Build SZ3...
if ${build[all]} || ${build[sz3]} ; then
    target="SZ3-3.2.1"
    unpack "$target"
    cd -
    mkdir -p "$build_dir"/src/"$target"/build_shared && cd "$build_dir"/src/"$target"/build_shared || exit
    cmake -DCMAKE_INSTALL_PREFIX="$build_dir" -DBUILD_SHARED_LIBS=ON .. && make -j $nprocs && make install || exit
    mkdir -p "$build_dir"/src/"$target"/build_static && cd "$build_dir"/src/"$target"/build_static || exit
    cmake -DCMAKE_INSTALL_PREFIX="$build_dir" -DBUILD_SHARED_LIBS=OFF .. && make -j $nprocs && make install || exit
fi

# Finish...
printf "All selected compilations are done.\n"
echo -e "\n***** gcc *****\n"
gcc --version
echo -e "\n***** gsl-config *****\n"
"$build_dir"/bin/gsl-config --libs --cflags --version
echo -e "\n***** nc-config *****"
"$build_dir"/bin/nc-config --all
