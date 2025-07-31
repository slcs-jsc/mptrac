#!/bin/bash

########## Compare Structures

# remove variable files if exist
rm -f data/var_*

# extract structures
awk '/TYPE, bind\(c\) :: ctl_t/{f=1;next} /END TYPE ctl_t/{f=0} f' <  ../../src/mptrac_fortran.f90 > data/test_ctl_f.asc
awk '/TYPE, bind\(c\) :: met_t/{f=1;next} /END TYPE met_t/{f=0} f' < ../../src/mptrac_fortran.f90 > data/test_met_f.asc
awk '$0=="  TYPE, bind(c) :: clim_t"{f=1;next} /END TYPE clim_t/{f=0} f' < ../../src/mptrac_fortran.f90 > data/test_clim_f.asc
awk '/TYPE, bind\(c\) :: clim_zm_t/{f=1;next} /END TYPE clim_zm_t/{f=0} f' < ../../src/mptrac_fortran.f90 > data/test_clim_zm_f.asc
awk '/TYPE, bind\(c\) :: clim_ts_t/{f=1;next} /END TYPE clim_ts_t/{f=0} f' < ../../src/mptrac_fortran.f90 > data/test_clim_ts_f.asc
awk '/TYPE, bind\(c\) :: clim_photo_t/{f=1;next} /END TYPE clim_photo_t/{f=0} f' < ../../src/mptrac_fortran.f90 > data/test_clim_photo_f.asc
awk '/TYPE, bind\(c\) :: atm_t/{f=1;next} /END TYPE atm_t/{f=0} f' < ../../src/mptrac_fortran.f90 > data/test_atm_f.asc
awk '/TYPE, bind\(c\) :: cache_t/{f=1;next} /END TYPE cache_t/{f=0} f' < ../../src/mptrac_fortran.f90 > data/test_cache_f.asc

awk '/brief Control parameters/{f=1;next} /\} ctl_t;/{f=0} f' < ../../src/mptrac.h > data/test_ctl_c.asc
awk '/Meteo data structure/{f=1;next} /\} met_t;/{f=0} f' < ../../src/mptrac.h > data/test_met_c.asc
awk '/climatological data containing/{f=1;next} /\} clim_t;/{f=0} f' < ../../src/mptrac.h > data/test_clim_c.asc
awk '/Climatological data in the form of zonal means./{f=1;next} /\} clim_zm_t;/{f=0} f' < ../../src/mptrac.h > data/test_clim_zm_c.asc
awk '/Climatological data in the form of time series./{f=1;next} /\} clim_ts_t;/{f=0} f' < ../../src/mptrac.h > data/test_clim_ts_c.asc
awk '/Climatological data in the form of photolysis rates./{f=1;next} /\} clim_photo_t;/{f=0} f' < ../../src/mptrac.h > data/test_clim_photo_c.asc
awk '/Air parcel data./{f=1;next} /\} atm_t;/{f=0} f' < ../../src/mptrac.h > data/test_atm_c.asc
awk '/Cache data structure./{f=1;next} /\} cache_t;/{f=0} f' < ../../src/mptrac.h > data/test_cache_c.asc

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


########## Compare Dimensions 

grep "ex =" ../../src/mptrac_fortran.f90 | awk '{print $6}' > y
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > yy; fi; done < y
rm y

grep "#define EX " ../../src/mptrac.h | awk '{print $3}' > x
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > xx; fi; done < x
rm x

if cmp -s xx yy; then
    echo "Dimension EX (C) and ex (Fortran) are equal"
    rm xx yy
else
    echo "Dimension EX (C) and ex (Fortran) are not equal. EX is:"
    cat xx
    echo "and ex is:"
    cat yy
    rm xx yy
    exit_code=99
    exit $exit_code
fi

grep "ey ="  ../../src/mptrac_fortran.f90 | awk '{print $6}' > y
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > yy; fi; done < y
rm y

grep "#define EY "  ../../src/mptrac.h | awk '{print $3}' > x
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > xx; fi; done < x
rm x

if cmp -s xx yy; then
    echo "Dimension EY (C) and ey (Fortran) are equal"
    rm xx yy
else
    echo "Dimension EY (C) and ey (Fortran) are not equal. EY is:"
    cat xx
    echo "and ey is:"
    cat yy
    rm xx yy
    exit_code=99
    exit $exit_code
fi

grep "ep ="  ../../src/mptrac_fortran.f90 | awk '{print $6}' > y
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > yy; fi; done < y
rm y

grep "#define EP "  ../../src/mptrac.h | awk '{print $3}' > x
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > xx; fi; done < x
rm x

if cmp -s xx yy; then
    echo "Dimension EP (C) and ep (Fortran) are equal"
    rm xx yy
else
    echo "Dimension EP (C) and ep (Fortran) are not equal. EP is:"
    cat xx
    echo "and ep is:"
    cat yy
    rm xx yy
    exit_code=99
    exit $exit_code
fi

grep "npp ="  ../../src/mptrac_fortran.f90 | awk '{print $6}' > y
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > yy; fi; done < y
rm y

grep "#define NP "  ../../src/mptrac.h | awk '{print $3}' > x
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > xx; fi; done < x
rm x

if cmp -s xx yy; then
    echo "Dimension NP (C) and npp (Fortran) are equal"
    rm xx yy
else
    echo "Dimension NP (C) and npp (Fortran) are not equal. NP is:"
    cat xx
    echo "and npp is:"
    cat yy
    rm xx yy
    exit_code=99
    exit $exit_code
fi

grep "nqq =" ../../src/mptrac_fortran.f90 | awk '{print $6}' > y
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > yy; fi; done < y
rm y

grep "#define NQ "  ../../src/mptrac.h | awk '{print $3}' > x
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > xx; fi; done < x
rm x

if cmp -s xx yy; then
    echo "Dimension NQ (C) and nqq (Fortran) are equal"
    rm xx yy
else
    echo "Dimension NQ (C) and nqq (Fortran) are not equal. NQ is:"
    cat xx
    echo "and nqq is:"
    cat yy
    rm xx yy
    exit_code=99
    exit $exit_code
fi

grep "length =" ../../src/mptrac_fortran.f90 | awk '{print $6}' > y
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > yy; fi; done < y
rm y

grep "#define LEN " ../../src/mptrac.h | awk '{print $3}' > x
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > xx; fi; done < x
rm x

if cmp -s xx yy; then
    echo "Dimension LEN (C) and length (Fortran) are equal"
    rm xx yy
else
    echo "Dimension LEN (C) and length (Fortran) are not equal. LEN is:"
    cat xx
    echo "and length is:"
    cat yy
    rm xx yy
    exit_code=99
    exit $exit_code
fi

grep "cyy =" ../../src/mptrac_fortran.f90 | awk '{print $6}' > y
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > yy; fi; done < y
rm y

grep "#define CY " ../../src/mptrac.h | awk '{print $3}' > x
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > xx; fi; done < x
rm x

if cmp -s xx yy; then
    echo "Dimension CY (C) and cyy (Fortran) are equal"
    rm xx yy
else
    echo "Dimension CY (C) and cyy (Fortran) are not equal. CY is:"
    cat xx
    echo "and cyy is:"
    cat yy
    rm xx yy
    exit_code=99
    exit $exit_code
fi

grep "co3 =" ../../src/mptrac_fortran.f90 | awk '{print $6}' > y
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > yy; fi; done < y
rm y

grep "#define CO3 " ../../src/mptrac.h | awk '{print $3}' > x
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > xx; fi; done < x
rm x

if cmp -s xx yy; then
    echo "Dimension CO3 (C) and co3 (Fortran) are equal"
    rm xx yy
else
    echo "Dimension CO3 (C) and co3 (Fortran) are not equal. CO3 is:"
    cat xx
    echo "and co3 is:"
    cat yy
    rm xx yy
    exit_code=99
    exit $exit_code
fi

grep "cp =" ../../src/mptrac_fortran.f90 | awk '{print $6}' > y
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > yy; fi; done < y
rm y

grep "#define CP " ../../src/mptrac.h | awk '{print $3}' > x
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > xx; fi; done < x
rm x

if cmp -s xx yy; then
    echo "Dimension CP (C) and cp (Fortran) are equal"
    rm xx yy
else
    echo "Dimension CP (C) and cp (Fortran) are not equal. CP is:"
    cat xx
    echo "and cp is:"
    cat yy
    rm xx yy
    exit_code=99
    exit $exit_code
fi

grep "csza =" ../../src/mptrac_fortran.f90 | awk '{print $6}' > y
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > yy; fi; done < y
rm y

grep "#define CSZA " ../../src/mptrac.h | awk '{print $3}' > x
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > xx; fi; done < x
rm x

if cmp -s xx yy; then
    echo "Dimension CSZA (C) and csza (Fortran) are equal"
    rm xx yy
else
    echo "Dimension CSZA (C) and csza (Fortran) are not equal. CSZA is:"
    cat xx
    echo "and csza is:"
    cat yy
    rm xx yy
    exit_code=99
    exit $exit_code
fi

echo "Difference ct, CT"
grep "ct =" ../../src/mptrac_fortran.f90 | awk '{print $6}' > y
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > yy; fi; done < y
rm y

grep "#define CT " ../../src/mptrac.h | awk '{print $3}' > x
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > xx; fi; done < x
rm x

if cmp -s xx yy; then
    echo "Dimension CT (C) and ct (Fortran) are equal"
    rm xx yy
else
    echo "Dimension CT (C) and ct (Fortran) are not equal. CT is:"
    cat xx
    echo "and ct is:"
    cat yy
    rm xx yy
    exit_code=99
    exit $exit_code
fi

grep "cts =" ../../src/mptrac_fortran.f90 | awk '{print $6}' > y
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > yy; fi; done < y
rm y

grep "#define CTS " ../../src/mptrac.h | awk '{print $3}' > x
while read line ; do  if [[ $line = *[[:"digit":]]* ]]; then echo $line > xx; fi; done < x
rm x

if cmp -s xx yy; then
    echo "Dimension CTS (C) and cts (Fortran) are equal"
    rm xx yy
else
    echo "Dimension CTS (C) and cts (Fortran) are not equal. CTS is:"
    cat xx
    echo "and cts is:"
    cat yy
    rm xx yy
    exit_code=99
    exit $exit_code
fi

