#! /bin/bash

# Setup...
rm -rf build
target=$(mkdir -p build && cd build && pwd)

# Prepare directories...
mkdir -p $target/src $target/bin $target/include $target/lib $target/man/man1 \
    && cp *tar.bz2 $target/src \
    && cd $target/src \
    && for f in $(ls *tar.bz2) ; do tar xvjf $f ; done \
	|| exit

# GSL...
dir=gsl-2.7
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make -j && make check && make install && make clean \
	|| exit

# Thrust sort...
dir=thrustsort-1.2
cd $target/src/$dir \
    && cp -a libthrustsort_gpu.a libthrustsort_cpu.a $target/lib/ \
	|| exit

# zlib...
dir=zlib-1.2.12
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make -j && make check && make install && make clean \
	|| exit

# zfp...
dir=zfp-0.5.5
cd $target/src/$dir \
    && make && make test \
    && cp -a bin/* $target/bin/ \
    && cp -a include/* $target/include/ \
    && cp -a lib/* $target/lib/ \
	|| exit

# zstd...
dir=zstd-1.5.2
cd $target/src/$dir \
    && make -j && make check \
    && cp -a lib/libzstd* $target/lib/ \
    && cp -a lib/*.h $target/include/ \
    && ln -sf $target/lib/libzstd.so.1.5.2 $target/lib/libzstd.so \
    && make clean \
	|| exit

# HDF5...
dir=hdf5-1.12.1
cd $target/src/$dir \
    && ./configure --prefix=$target --with-zlib=$target --enable-hl \
    && make -j ; make -j && make check && make install && make clean \
	|| exit

# netCDF...
dir=netcdf-c-4.8.1
cd $target/src/$dir \
    && CPPFLAGS=-I$target/include LDFLAGS=-L$target/lib ./configure --prefix=$target --disable-dap --disable-nczarr \
    && make -j && make check && make install && make clean \
	|| exit

# # Chemistry using KPP...
# dir=KPP
# export KPP_HOME=$target/src/$dir
# export KPP_FLEX_LIB_DIR=$target/src/$dir
# cd $target/src/$dir/mptrac-chem \
# 		&& ../bin/kpp chem.kpp && make lib && cp libkpp.a $target/lib && cp *.h $target/include

# Summary...
echo -e "\n***** gsl-config *****\n"
$target/bin/gsl-config --libs --cflags --version
echo -e "\n***** nc-config *****"
$target/bin/nc-config --all
