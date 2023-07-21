#! /bin/bash

# Check arguments...
if [ $# -ne 1 ] ; then
    echo "usage: $0 <kppdir>"
    exit
fi

# Set directory...
kppdir=$1

# Copy Makefile and header file...
cp ./Makefile $kppdir/
cp ./kpp_chem.h $kppdir/

# Set target directory...
target=$(cd ../libs/build && pwd)
export KPP_HOME=$target/src/KPP
export KPP_FLEX_LIB_DIR=$target/src/KPP

# Build KPP code...
cd $kppdir  
$KPP_HOME/bin/kpp chem.kpp \
    && make clean \
    && make lib INFO=0 \
    && cp libkpp.a $target/lib \
    && cp *.h $target/include \
    && rm *.o
