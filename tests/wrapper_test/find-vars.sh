#!/bin/bash

########## Compare Structures

# remove variable files if exist
rm -f data/var_*

# Preprocess sources with the same user-defined macros used for the build.
fortran_src=data/mptrac_fortran.pp.f90
c_defines=data/mptrac_defines_c.asc
${FC:-gfortran} ${FPPFLAGS:--cpp -E -P} ../../src/mptrac_fortran.f90 > "$fortran_src"
${CC:-gcc} -E -dM ${DEFINES} ../../src/mptrac.h > "$c_defines"

# extract structures
awk '/TYPE, bind\(c\) :: ctl_t/{f=1;next} /END TYPE ctl_t/{f=0} f' < "$fortran_src" > data/test_ctl_f.asc
awk '/TYPE, bind\(c\) :: met_t/{f=1;next} /END TYPE met_t/{f=0} f' < "$fortran_src" > data/test_met_f.asc
awk '$0=="  TYPE, bind(c) :: clim_t"{f=1;next} /END TYPE clim_t/{f=0} f' < "$fortran_src" > data/test_clim_f.asc
awk '/TYPE, bind\(c\) :: clim_zm_t/{f=1;next} /END TYPE clim_zm_t/{f=0} f' < "$fortran_src" > data/test_clim_zm_f.asc
awk '/TYPE, bind\(c\) :: clim_ts_t/{f=1;next} /END TYPE clim_ts_t/{f=0} f' < "$fortran_src" > data/test_clim_ts_f.asc
awk '/TYPE, bind\(c\) :: clim_photo_t/{f=1;next} /END TYPE clim_photo_t/{f=0} f' < "$fortran_src" > data/test_clim_photo_f.asc
awk '/TYPE, bind\(c\) :: atm_t/{f=1;next} /END TYPE atm_t/{f=0} f' < "$fortran_src" > data/test_atm_f.asc
awk '/TYPE, bind\(c\) :: cache_t/{f=1;next} /END TYPE cache_t/{f=0} f' < "$fortran_src" > data/test_cache_f.asc
awk '/TYPE, bind\(c\) :: dd_t/{f=1;next} /END TYPE dd_t/{f=0} f' < "$fortran_src" > data/test_dd_f.asc

awk '/brief Control parameters/{f=1;next} /\} ctl_t;/{f=0} f' < ../../src/mptrac.h > data/test_ctl_c.asc
awk '/Meteo data structure/{f=1;next} /\} met_t;/{f=0} f' < ../../src/mptrac.h > data/test_met_c.asc
awk '/climatological data containing/{f=1;next} /\} clim_t;/{f=0} f' < ../../src/mptrac.h > data/test_clim_c.asc
awk '/Climatological data in the form of zonal means./{f=1;next} /\} clim_zm_t;/{f=0} f' < ../../src/mptrac.h > data/test_clim_zm_c.asc
awk '/Climatological data in the form of time series./{f=1;next} /\} clim_ts_t;/{f=0} f' < ../../src/mptrac.h > data/test_clim_ts_c.asc
awk '/Climatological data in the form of photolysis rates./{f=1;next} /\} clim_photo_t;/{f=0} f' < ../../src/mptrac.h > data/test_clim_photo_c.asc
awk '/Air parcel data./{f=1;next} /\} atm_t;/{f=0} f' < ../../src/mptrac.h > data/test_atm_c.asc
awk '/Cache data structure./{f=1;next} /\} cache_t;/{f=0} f' < ../../src/mptrac.h > data/test_cache_c.asc
awk '/Domain decomposition data structure./{f=1;skip=0;next} /\} dd_t;/{f=0} f {if ($0 ~ /#ifdef DD/) {skip=1; next}; if ($0 ~ /#endif/) {skip=0; next}; if (!skip) print}' < ../../src/mptrac.h > data/test_dd_c.asc

# compare number of structure entries and variables inside 
grep ";" data/test_ctl_c.asc | wc -l > x
wc -l < data/test_ctl_f.asc > y
for varc in $(grep ";" data/test_ctl_c.asc| awk '{print $2}') ; do
    echo "$varc" | awk -F'[;[]' '{print $1}' >> data/var_ctl_c.asc
done

for varf in $(grep "::" data/test_ctl_f.asc|awk '{print $NF}') ; do
    echo "$varf" >> data/var_ctl_f.asc
done
if cmp -s x y && cmp -s data/var_ctl_c.asc data/var_ctl_f.asc; then
    echo "Number of control parameters are equal:"
    cat x
    rm x y
else
    echo "Number of control parameters are not equal. C has: "
    cat x
    echo "and Fortran has: "
    cat y
    echo "Difference between var_ctl_c.asc and var_ctl_f.asc"
    diff data/var_ctl_c.asc data/var_ctl_f.asc
    rm x y
    exit_code=99
    exit $exit_code
fi

grep ";" data/test_met_c.asc | wc -l > x
wc -l < data/test_met_f.asc > y
for varc in $(grep ";" data/test_met_c.asc| awk '{print $2}') ; do
    echo "$varc" |awk -F'[;[]' '{print $1}' >> data/var_met_c.asc
done

for varf in $(grep "::" data/test_met_f.asc|awk '{print $NF}') ; do
    echo "$varf" >> data/var_met_f.asc
done

if cmp -s x y && cmp -s data/var_met_c.asc data/var_met_f.asc; then
    echo "Number of meteo data structures are equal:"
    cat x
    rm x y
else
    echo "Number of meteo data structures are not equal. C has: "
    cat x
    echo "and Fortran has: "
    cat y
    echo "Difference between var_met_c.asc and var_met_f.asc"
    diff data/var_met_c.asc data/var_met_f.asc
    rm x y
    exit_code=99
    exit $exit_code
fi

grep ";" data/test_clim_c.asc | wc -l > x
wc -l < data/test_clim_f.asc > y
for varc in $(grep ";" data/test_clim_c.asc| awk '{print $2}') ; do
    echo "$varc" | awk -F'[;[]' '{print $1}' >> data/var_clim_c.asc
done

for varf in $(grep "::" data/test_clim_f.asc|awk '{print $NF}') ; do
    echo "$varf" >> data/var_clim_f.asc
done
if cmp -s x y && cmp -s data/var_clim_c.asc data/var_clim_f.asc; then
    echo "Number of climatological data structures are equal:"
    cat x
    rm x y
else
    echo "Number of climatological data structures are not equal. C has: "
    cat x
    echo "and Fortran has: "
    cat y
    echo "Difference between var_clim_c.asc and var_clim_f.asc"
    diff data/var_clim_c.asc data/var_clim_f.asc
    rm x y
    exit_code=99
    exit $exit_code
fi

grep ";" data/test_clim_zm_c.asc | wc -l > x
wc -l < data/test_clim_zm_f.asc > y
for varc in $(grep ";" data/test_clim_zm_c.asc| awk '{print $2}') ; do
    echo "$varc" | awk -F'[;[]' '{print $1}' >> data/var_clim_zm_c.asc
done

for varf in $(grep "::" data/test_clim_zm_f.asc|awk '{print $NF}') ; do
    echo "$varf" >> data/var_clim_zm_f.asc
done
if cmp -s x y && cmp -s data/var_clim_zm_c.asc data/var_clim_zm_f.asc; then
    echo "Number of climatological data in the form of zonal means are equal:"
    cat x
    rm x y
else
    echo "Number of climatological datain the form of zonal means are not equal. C has: "
    cat x
    echo "and Fortran has: "
    cat y
    echo "Difference between var_clim_zm_c.asc and var_clim_zm_f.asc"
    diff data/var_clim_zm_c.asc data/var_clim_zm_f.asc
    rm x y
    exit_code=99
    exit $exit_code
fi

grep ";" data/test_clim_ts_c.asc | wc -l > x
wc -l < data/test_clim_ts_f.asc > y
for varc in $(grep ";" data/test_clim_ts_c.asc| awk '{print $2}') ; do
    echo "$varc" | awk -F'[;[]' '{print $1}' >> data/var_clim_ts_c.asc
done

for varf in $(grep "::" data/test_clim_ts_f.asc|awk '{print $NF}') ; do
    echo "$varf" >> data/var_clim_ts_f.asc
done

if cmp -s x y && cmp -s data/var_clim_ts_c.asc data/var_clim_ts_f.asc; then
    echo "Number of climatological data in the form of time series are equal:"
    cat x
    rm x y
else
    echo "Number of climatological datain the form of time series are not equal. C has: "
    cat x
    echo "and Fortran has: "
    cat y
    echo "Difference between var_clim_ts_c.asc and var_clim_ts_f.asc"
    diff data/var_clim_ts_c.asc data/var_clim_ts_f.asc
    rm x y    
    exit_code=99
    exit $exit_code
fi

grep ";" data/test_clim_photo_c.asc | wc -l > x
wc -l < data/test_clim_photo_f.asc > y
for varc in $(grep ";" data/test_clim_photo_c.asc| awk '{print $2}') ; do
    echo "$varc" | awk -F'[;[]' '{print $1}' >> data/var_clim_photo_c.asc
done

for varf in $(grep "::" data/test_clim_photo_f.asc|awk '{print $NF}') ; do
    echo "$varf" >> data/var_clim_photo_f.asc
done
if cmp -s x y && cmp -s data/var_clim_photo_c.asc data/var_clim_photo_f.asc ; then
    echo "Number of climatological data in the form of photolysis rates are equal:"
    cat x
    rm x y
else
    echo "Number of climatological data in the form of photolysis rates are not equal. C has: "
    cat x
    echo "and Fortran has: "
    cat y
    echo "Difference between var_clim_photo_c.asc and var_clim_photo_f.asc"
    diff data/var_clim_photo_c.asc data/var_clim_photo_f.asc
    rm x y
    exit_code=99
    exit $exit_code
fi

grep ";" data/test_atm_c.asc | wc -l > x
wc -l < data/test_atm_f.asc > y
for varc in $(grep ";" data/test_atm_c.asc| awk '{print $2}') ; do
    echo "$varc" |awk -F'[;[]' '{print $1}' >> data/var_atm_c.asc
done

for varf in $(grep "::" data/test_atm_f.asc|awk '{print $NF}') ; do
    echo "$varf" >> data/var_atm_f.asc
done
if cmp -s x y && cmp -s data/var_atm_c.asc data/var_atm_f.asc ; then
    echo "Number of air parcel data are equal:"
    cat x
    rm x y
else
    echo "Number of air parcel data are not equal. C has: "
    cat x
    echo "and Fortran has: "
    cat y
    echo "Difference between var_atm_c.asc and var_atm_f.asc"
    diff data/var_atm_c.asc data/var_atm_f.asc
    rm x y
    exit_code=99
    exit $exit_code
fi

grep ";" data/test_cache_c.asc | wc -l > x
wc -l < data/test_cache_f.asc > y
for varc in $(grep ";" data/test_cache_c.asc| awk '{print $2}') ; do
    echo "$varc" | awk -F'[;[]' '{print $1}' >> data/var_cache_c.asc
done

for varf in $(grep "::" data/test_cache_f.asc|awk '{print $NF}') ; do
    echo "$varf" >> data/var_cache_f.asc
done
if cmp -s x y && cmp -s data/var_cache_c.asc data/var_cache_f.asc ; then
    echo "Number of cache data are equal:"
    cat x
    rm x y
else
    echo "Number of cache data are not equal. C has: "
    cat x
    echo "and Fortran has: "
    cat y
    echo "Difference between var_cache_c.asc and var_cache_f.asc"
    diff data/var_cache_c.asc data/var_cache_f.asc
    rm x y
    exit_code=99
    exit $exit_code
fi

grep ";" data/test_dd_c.asc | wc -l > x
wc -l < data/test_dd_f.asc > y
for varc in $(grep ";" data/test_dd_c.asc | awk '{print $2}') ; do
    echo "$varc" | awk -F'[;[]' '{print $1}' >> data/var_dd_c.asc
done

for varf in $(grep "::" data/test_dd_f.asc|awk '{print $NF}') ; do
    echo "$varf" >> data/var_dd_f.asc
done
if cmp -s x y && cmp -s data/var_dd_c.asc data/var_dd_f.asc ; then
    echo "Number of domain decomposition data are equal:"
    cat x
    rm x y
else
    echo "Number of domain decomposition data are not equal. C has: "
    cat x
    echo "and Fortran has: "
    cat y
    echo "Difference between var_dd_c.asc and var_dd_f.asc"
    diff data/var_dd_c.asc data/var_dd_f.asc
    rm x y
    exit_code=99
    exit $exit_code
fi


########## Compare Dimensions 

check_dimension() {
    macro=$1
    fvar=$2
    grep "^#define $macro " "$c_defines" | awk '{print $3}' > xx
    grep -i "^[[:space:]]*INTEGER, PARAMETER :: $fvar =" "$fortran_src" | awk '{print $6}' > yy

    if cmp -s xx yy; then
	echo "Dimension $macro (C) and $fvar (Fortran) are equal"
	rm xx yy
    else
	echo "Dimension $macro (C) and $fvar (Fortran) are not equal. $macro is:"
	cat xx
	echo "and $fvar is:"
	cat yy
	rm xx yy
	exit_code=99
	exit $exit_code
    fi
}

check_dimension EX ex
check_dimension EY ey
check_dimension EP ep
check_dimension METVAR metvar
check_dimension NP npp
check_dimension NQ nqq
check_dimension LEN length
check_dimension CY cyy
check_dimension CO3 co3
check_dimension CP cp
check_dimension CSZA csza
check_dimension CT ct
check_dimension CTS cts
check_dimension DD_NNMAX dd_nnmax
