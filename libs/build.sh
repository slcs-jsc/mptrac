#! /bin/bash
# ---------------------------------------------------------------------- # 
#
# Description
#     This script is based on the original build script by Lars Hoffmann
# in order to build the external libraries used by mptrac. The user's
# manual of this script can be found by executing with -h
#
# ---------------------------------------------------------------------- # 
#
# DevLog:
#   03.08.2023: Starting the modification
#
# ---------------------------------------------------------------------- #

###### Presets ######

ifBuildAll=false
ifBuildGSL=false
ifBuildTS=false
ifBuildZLIB=false
ifBuildZFP=false
ifBuildZSTD=false
ifBuildHDF5=false
ifBuildNETCDF=false
ifSkip=false
numProcs=4

###### Reminder           ######

if [[ $# -lt 1 ]]; then
    printf "please using -h to open the manual\n"
    ifSkip=true
fi

###### Checking the flags ######

while getopts agtzfshnp: flag
do
    case "${flag}" in
        a) ifBuildAll=true    ; echo "All the libs will be built"    ;;	
	g) ifBuildGSL=true    ; echo "GSL"     ;;
	t) ifBuildTS=true     ; echo "TS "     ;;
	z) ifBuildZLIB=true   ; echo "ZLIB"    ;;
	f) ifBuildZFP=true    ; echo "ZFP"     ;;
	s) ifBuildZSTD=true   ; echo "ZSTD"    ;;
	d) ifBuildHDF5=true   ; echo "HDF5"    ;;
	n) ifBuildNETCDF=true ; echo "NETCDF"  ;;
	p) numProcs=${OPTARG}  
           numTotalProcs=`nproc --all`
           if [ $numProcs -lt $numTotalProcs ]; then
                numProcs=$numProcs
           else
		numProcs=$(($numTotalProcs*2/3))
           fi
	   echo "NEW PROCESSES as ${numProcs}" ;;
	h) printf -- "Helping Page: \n"
	   printf -- "    Please specify which libs ( or all, -a ) you want to build\n"
	   printf -- "    you can specify within a dash, e.g.: -gtzfs to build      \n"
	   printf -- "    without NETCDF and HDF5                                   \n" 
	   printf -- "---------------------------------------------------------------"	
           printf -- "-a         : build all the libs     \n"
	   printf -- "-g         : build GSL    \n" 
	   printf -- "-t         : build TS     \n" 
	   printf -- "-z         : build ZLIB   \n" 
	   printf -- "-f         : build ZFP    \n" 
	   printf -- "-s         : build ZSTD   \n" 
	   printf -- "-d         : build HDF5   (prerequisite: ZLIB)  \n" 
	   printf -- "-n         : build NETCDF (prerequisite: HDF5)  \n" 
	   printf -- "-p [cores] : how many [cores] to be used by make -j [cores]  \n" 
	   printf -- "           : default is 4                                    \n" 
	   ifSkip=true
	   ;;	 
    esac
done

###### 

#printf "processes: %d \n"  $numProcs
#printf "building the following Libs:" 

######


if [ $ifSkip = false ]; then
    printf "preparing the files \n"
    # Setup...
    rm -rf build
    target=$(mkdir -p build && cd build && pwd)
    
    printf "Using %d CPUs for make ... \n" $numProcs
    
    # Prepare directories...
    mkdir -p $target/src $target/bin $target/include $target/lib $target/man/man1 \
        && cp *tar.bz2 $target/src \
        && cd $target/src \
        && for f in $(ls *tar.bz2) ; do tar xvjf $f ; done \
    	|| exit
fi 
if [ $ifBuildAll = true ] || [ $ifBuildGSL = true ] ; then 
    printf "\n"
    # GSL...
    dir=gsl-2.7
    cd $target/src/$dir \
        && ./configure --prefix=$target \
        && make -j $numProcs && make check && make install && make clean \
    	|| exit
fi

if [ $ifBuildAll = true ] || [ $ifBuildTS = true ] ; then 
    printf "\n"
    # Thrust sort...
    dir=thrustsort-1.2
    cd $target/src/$dir \
        && cp -a libthrustsort_gpu.a libthrustsort_cpu.a $target/lib/ \
    	|| exit
fi

if [ $ifBuildAll = true ] || [ $ifBuildZLIB = true ] ; then 
    printf "\n"
    # zlib...
    dir=zlib-1.2.12
    cd $target/src/$dir \
        && ./configure --prefix=$target \
        && make -j $numProcs && make check && make install && make clean \
    	|| exit
fi

if [ $ifBuildAll = true ] || [ $ifBuildZFP = true ] ; then 
    printf "\n"
    # zfp...
    dir=zfp-0.5.5
    cd $target/src/$dir \
        && make && make test \
        && cp -a bin/* $target/bin/ \
        && cp -a include/* $target/include/ \
        && cp -a lib/* $target/lib/ \
    	|| exit
fi

if [ $ifBuildAll = true ] || [ $ifBuildZSTD = true ] ; then 
    printf "\n"
    # zstd...
    dir=zstd-1.5.2
    cd $target/src/$dir \
        && make -j $numProcs && make check \
        && cp -a lib/libzstd* $target/lib/ \
        && cp -a lib/*.h $target/include/ \
        && ln -sf $target/lib/libzstd.so.1.5.2 $target/lib/libzstd.so \
        && make clean \
    	|| exit
fi

if [ $ifBuildAll = true ] || [ $ifBuildHDF5 = true ] ; then 
    printf "\n"
    # HDF5...
    dir=hdf5-1.12.1
    cd $target/src/$dir \
        && ./configure --prefix=$target --with-zlib=$target --enable-hl \
        && make -j $numProcs && make check && make install && make clean \
    	|| exit
fi

if [ $ifBuildAll = true ] || [ $ifBuildNETCDF = true ] ; then 
    printf "starting to compile netCDF \n"
    # netCDF...
    dir=netcdf-c-4.8.1
    cd $target/src/$dir \
        && CPPFLAGS=-I$target/include LDFLAGS=-L$target/lib ./configure --prefix=$target --disable-dap --disable-nczarr \
        && make -j $numProcs && make check && make install && make clean \
    	|| exit
fi

    # # Chemistry using KPP...
    # dir=KPP
    # export KPP_HOME=$target/src/$dir
    # export KPP_FLEX_LIB_DIR=$target/src/$dir
    # cd $target/src/$dir/mptrac-chem \
    # 		&& ../bin/kpp chem.kpp && make lib && cp libkpp.a $target/lib && cp *.h $target/include
    

if [ $ifSkip = false ] ; then 
    printf "Summary of the compilation \n"
    # Summary...
    echo -e "\n***** gsl-config *****\n"
    $target/bin/gsl-config --libs --cflags --version
    echo -e "\n***** nc-config *****"
    $target/bin/nc-config --all
    printf "All selected compilations are done"
fi
