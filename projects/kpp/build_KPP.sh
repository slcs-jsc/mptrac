#! /bin/bash

# Check arguments...
if [ $# -ne 1 ] ; then
    echo "usage: $0 <kppdir>"
    exit
fi

# Set directory...
kppdir=$1
REBUILD=1



# Set target directory...
target=$(cd ../../libs/build && pwd)
export KPP_HOME=$target/src/KPP
export KPP_FLEX_LIB_DIR=$target/src/KPP

# Build KPP code... 
cd $kppdir
[ $REBUILD -eq 1 ] && cp ../Makefile ./ \
&& cp ../kpp_chem.h ./ \
&& $KPP_HOME/bin/kpp chem.kpp 

make clean \
&& make lib INFO=0 \
&& cp libkpp.a $target/lib \
&& cp *.h $target/include \
&& rm *.o
