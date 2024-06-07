#! /bin/bash

# Check arguments...
if [ $# -ne 2 ] ; then
    echo "usage: $0 <kppdir> <1:GPU on; 0:GPU off>"
    exit
fi

# Set directory...
kppdir=$1
GPU=$2
REBUILD=1

# Set target directory...
target=$(cd ../../libs/build && pwd)
export KPP_HOME=$target/src/KPP
export KPP_FLEX_LIB_DIR=$target/src/KPP

# Build KPP code... 
cd $kppdir
[ $REBUILD -eq 1 ] && cp ../Makefile ./Makefile \
&& cp ../kpp_chem.h ./kpp_chem.h \
&& $KPP_HOME/bin/kpp chem.kpp 

# [ $kppdir = 'chem2/' ] && \
# cp ../chem_Integrator2.c ./chem_Integrator.c

make clean \
&& make lib INFO=0 GPU=$GPU\
&& cp libkpp.a $target/lib \
&& cp ./*.h $target/include \
&& rm ./*.o
