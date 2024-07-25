#!/bin/bash

#ToDo expand script:
# add all other data types, clim_t, clim_zm_t, clim_zt_t, clim_photo_t
# add check for constants
# add check for dimensions

#extract structures
cat mptrac_fortran.f90| awk '/TYPE, bind\(c\) :: ctl_t/{f=1;next} /END TYPE ctl_t/{f=0} f' > test_ctl_f.asc
cat mptrac_fortran.f90| awk '/TYPE, bind\(c\) :: met_t/{f=1;next} /END TYPE met_t/{f=0} f' > test_met_f.asc

cat mptrac.h| awk '/Control parameters/{f=1;next} /\} ctl_t;/{f=0} f' > test_ctl_c.asc
cat mptrac.h| awk '/Meteo data structure/{f=1;next} /\} met_t;/{f=0} f' > test_met_c.asc


cat test_ctl_c.asc| grep ";"| wc -l
cat test_ctl_f.asc|wc -l
cat test_met_c.asc| grep ";"| wc -l
cat test_met_f.asc|wc -l 

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
