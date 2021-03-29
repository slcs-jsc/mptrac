#! /bin/bash

# Setup...
target=$(rm -rf build && mkdir -p build && cd build && pwd)
threads=$(cat /proc/cpuinfo | grep processor | wc -l)

# Download libraries...
[ -s gsl-2.6.tar.gz ] || wget https://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz
[ -s netcdf-c-4.7.4.tar.gz ] || wget https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-c-4.7.4.tar.gz

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

# netCDF...
dir=netcdf-c-4.7.4
cd $target/src/$dir \
    && ./configure --prefix=$target --disable-dap --disable-netcdf-4 \
    && make -j$threads && make check && make install && make clean \
	|| exit
