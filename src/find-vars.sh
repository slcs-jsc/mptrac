#!/bin/bash

#ToDo expand script:
# add check for constants
# add check for dimensions

# remove variable files if exist
rm -f var_*

# extract structures
cat mptrac_fortran.f90| awk '/TYPE, bind\(c\) :: ctl_t/{f=1;next} /END TYPE ctl_t/{f=0} f' > test_ctl_f.asc
cat mptrac_fortran.f90| awk '/TYPE, bind\(c\) :: met_t/{f=1;next} /END TYPE met_t/{f=0} f' > test_met_f.asc
cat mptrac_fortran.f90| awk '$0=="  TYPE, bind(c) :: clim_t"{f=1;next} /END TYPE clim_t/{f=0} f'  > test_clim_f.asc
cat mptrac_fortran.f90| awk '/TYPE, bind\(c\) :: clim_zm_t/{f=1;next} /END TYPE clim_zm_t/{f=0} f' > test_clim_zm_f.asc
cat mptrac_fortran.f90| awk '/TYPE, bind\(c\) :: clim_ts_t/{f=1;next} /END TYPE clim_ts_t/{f=0} f' > test_clim_ts_f.asc
cat mptrac_fortran.f90| awk '/TYPE, bind\(c\) :: clim_photo_t/{f=1;next} /END TYPE clim_photo_t/{f=0} f' > test_clim_photo_f.asc

cat mptrac.h| awk '/Control parameters/{f=1;next} /\} ctl_t;/{f=0} f' > test_ctl_c.asc
cat mptrac.h| awk '/Meteo data structure/{f=1;next} /\} met_t;/{f=0} f' > test_met_c.asc
cat mptrac.h| awk '/climatological data containing/{f=1;next} /\} clim_t;/{f=0} f' > test_clim_c.asc
cat mptrac.h| awk '/Climatological data in the form of zonal means./{f=1;next} /\} clim_zm_t;/{f=0} f' > test_clim_zm_c.asc
cat mptrac.h| awk '/Climatological data in the form of time series./{f=1;next} /\} clim_ts_t;/{f=0} f' > test_clim_ts_c.asc
cat mptrac.h| awk '/Climatological data in the form of photolysis rates./{f=1;next} /\} clim_photo_t;/{f=0} f' > test_clim_photo_c.asc

echo "Compare lines of structure files"
echo "Control paramters C"
cat test_ctl_c.asc| grep ";"| wc -l
echo "Control paramters Fortran"
cat test_ctl_f.asc|wc -l
echo "Meteo data structure C"
cat test_met_c.asc| grep ";"| wc -l
echo "Meteo data structure Fortran"
cat test_met_f.asc|wc -l
echo "Climatological data structure C"
cat test_clim_c.asc| grep ";"| wc -l
echo "Climatological data structure Fortran"
cat test_clim_f.asc|wc -l
echo "Climatological data in the form of zonal means C"
cat test_clim_zm_c.asc| grep ";"| wc -l
echo "Climatological data in the form of zonal means Fortran"
cat test_clim_zm_f.asc|wc -l
echo "Climatological data in the form of time series C"
cat test_clim_ts_c.asc| grep ";"| wc -l
echo "Climatological data in the form of time series Fortran"
cat test_clim_ts_f.asc|wc -l
echo "Climatological data in the form of photolysis rates C"
cat test_clim_photo_c.asc| grep ";"| wc -l
echo "Climatological data in the form of photolysis rates Fortran"
cat test_clim_photo_f.asc|wc -l

echo "Diff between var_ctl_c.asc and var_ctl_f.asc"
for vari in $(grep ";" test_ctl_c.asc| awk '{print $2}')
do
    echo ${vari}|awk -F'[;[]' '{print $1}' >> var_ctl_c.asc
    #mystring=$(echo ${vari}|awk -F'[;[]' '{print $1}')
done

for varf in $(grep "::" test_ctl_f.asc|awk '{print $NF}')
do
    echo $varf >> var_ctl_f.asc
done

diff var_ctl_c.asc var_ctl_f.asc

# ########################## met

echo "Diff between var_met_c.asc and var_met_f.asc"
for vari in $(grep ";" test_met_c.asc| awk '{print $2}')
do
    echo ${vari}|awk -F'[;[]' '{print $1}' >> var_met_c.asc
    #mystring=$(echo ${vari}|awk -F'[;[]' '{print $1}')
done

for varf in $(grep "::" test_met_f.asc|awk '{print $NF}')
do
    echo $varf >> var_met_f.asc
done

diff var_met_c.asc var_met_f.asc

# ########################## clim

echo "Diff between var_clim_c.asc and var_clim_f.asc"
for vari in $(grep ";" test_clim_c.asc| awk '{print $2}')
do
    echo ${vari}|awk -F'[;[]' '{print $1}' >> var_clim_c.asc
    #mystring=$(echo ${vari}|awk -F'[;[]' '{print $1}')
done

for varf in $(grep "::" test_clim_f.asc|awk '{print $NF}')
do
    echo $varf >> var_clim_f.asc
done

diff var_clim_c.asc var_clim_f.asc

# ########################## clim_zm

echo "Diff between var_clim_zm_c.asc and var_clim_zm_f.asc"
for vari in $(grep ";" test_clim_zm_c.asc| awk '{print $2}')
do
    echo ${vari}|awk -F'[;[]' '{print $1}' >> var_clim_zm_c.asc
    #mystring=$(echo ${vari}|awk -F'[;[]' '{print $1}')
done

for varf in $(grep "::" test_clim_zm_f.asc|awk '{print $NF}')
do
    echo $varf >> var_clim_zm_f.asc
done

diff var_clim_zm_c.asc var_clim_zm_f.asc

# ########################## clim_ts

echo "Diff between var_clim_ts_c.asc and var_clim_ts_f.asc"
for vari in $(grep ";" test_clim_ts_c.asc| awk '{print $2}')
do
    echo ${vari}|awk -F'[;[]' '{print $1}' >> var_clim_ts_c.asc
    #mystring=$(echo ${vari}|awk -F'[;[]' '{print $1}')
done

for varf in $(grep "::" test_clim_ts_f.asc|awk '{print $NF}')
do
    echo $varf >> var_clim_ts_f.asc
done

diff var_clim_ts_c.asc var_clim_ts_f.asc

# ########################## clim_photo

echo "Diff between var_clim_photo_c.asc and var_clim_photo_f.asc"
for vari in $(grep ";" test_clim_photo_c.asc| awk '{print $2}')
do
    echo ${vari}|awk -F'[;[]' '{print $1}' >> var_clim_photo_c.asc
    #mystring=$(echo ${vari}|awk -F'[;[]' '{print $1}')
done

for varf in $(grep "::" test_clim_photo_f.asc|awk '{print $NF}')
do
    echo $varf >> var_clim_photo_f.asc
done

diff var_clim_photo_c.asc var_clim_photo_f.asc
