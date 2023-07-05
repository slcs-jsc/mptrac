#! /bin/bash

if [ $# -ne 1 ] ; then
    echo "usage: $0 <kppdir>"
    exit
fi

kppdir=$1

target=$(cd ../libs/build && pwd)
export KPP_HOME=$target/src/KPP
export KPP_FLEX_LIB_DIR=$target/src/KPP

cd $kppdir  
$KPP_HOME/bin/kpp chem.kpp && make clean && make lib INFO=0 && cp libkpp.a $target/lib && cp *.h $target/include && rm *.o