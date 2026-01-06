#! /bin/bash

# Check arguments...
if [ $# -ne 2 ] ; then
    echo "Usage: $0 <kppdir> <1:GPU on; 0:GPU off>"
    echo "Example: ./build.sh chem 0"
    echo "Note: now the GPU code of kpp_chem module is still under testing and not finished yet."
    exit
fi

# Set directory...
kppdir=$1
GPU=$2
REBUILD=1

# Set target directory...
target=$(cd ../../libs/build && pwd)
export KPP_HOME="$target"/src/KPP
export KPP_FLEX_LIB_DIR="$target"/src/KPP

# Build KPP code... 
cd "$kppdir" || exit
[ $REBUILD -eq 1 ] && cp ../Makefile ./Makefile \
    && cp ../kpp_chem.h ./kpp_chem.h \
    && "$KPP_HOME"/bin/kpp chem.kpp 

make clean \
    && make lib INFO=0 GPU=$GPU\
    && cp libkpp.a "$target"/lib \
    && cp ./*.h "$target"/include \
    && rm ./*.o
