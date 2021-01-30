#! /bin/bash

# Check arguments...
if [ $# -ne 1 ] ; then
    echo "usage: $0 <nc3|nc4>"
    exit
fi

# Setup...
target=$(rm -rf build && mkdir -p build && cd build && pwd)
threads=$(cat /proc/cpuinfo | grep processor | wc -l)

# Prepare directories...
mkdir -p $target/src $target/bin $target/lib $target/man/man1 \
    && cp *tar.gz $target/src \
    && cd $target/src \
    && for f in $(ls *tar.gz) ; do tar xvzf $f ; done \
	|| exit

# GSL...
dir=gsl-2.6
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make -j$threads && make check && make install && make clean \
	|| exit

# Check netCDF flag...
if [ $1 = "nc4" ] ; then
    
    # zlib...
    dir=zlib-1.2.11
    cd $target/src/$dir \
	&& ./configure --prefix=$target \
	&& make -j$threads && make check && make install && make clean \
	    || exit
    
    # HDF5...
    dir=hdf5-1.12.0
    cd $target/src/$dir \
	&& ./configure --prefix=$target --with-zlib=$target --enable-hl \
	&& make -j$threads ; make -j$threads && make check && make install && make clean \
	    || exit
fi

# netCDF...
dir=netcdf-c-4.7.4
cd $target/src/$dir \
    && CPPFLAGS=-I$target/include LDFLAGS=-L$target/lib ./configure --prefix=$target --disable-dap $([ $1 != "nc4" ] && echo "--disable-netcdf-4") \
    && make -j$threads && make check && make install && make clean \
	|| exit
