#! /bin/bash

# directory of the trajectory code...
trac=../../src

# directory of the reanalysis data...
metbase=../data/clams/erai

# directory for the output...
output=./data

# directory of the reference and init...
input=./data.ref/

# start and end date
start_date=16070100
end_date=16070106

year=20${start_date:0:2}
month=${start_date:2:2}
day=${start_date:4:2}
hour=${start_date:6:6}

end_year=20${end_date:0:2}
end_month=${end_date:2:2}
end_day=${end_date:4:2}
end_hour=${end_date:6:6}

t0=$(${trac}/time2jsec ${year} ${month} ${day} ${hour} 0 0 0)
t1=$(${trac}/time2jsec ${end_year} ${end_month} ${end_day} ${end_hour} 0 0 0)

# Create control parameter file...
cat > ${output}/trac.ctl <<EOF
MET_CONVENTION = 1
CPL_ZETA_PRESS_MODULES = 0
INTERPOLATION_TYPE = 0
PRESS_LEVEL_DEF = 2
ATM_TYPE = 3
ATM_TYPE_OUT = 0
ADVECT = 2
CLAMS_MET_DATA = 1
VERT_COORD_AP = 1
VERT_COORD_MET = 1
NQ = 7
QNT_NAME[0] = theta
QNT_NAME[1] = pv
QNT_NAME[2] = m
QNT_NAME[3] = zeta
QNT_NAME[4] = zeta_d
QNT_NAME[5] = ps
QNT_NAME[6] = p
METBASE = ${metbase}
DIRECTION = 1
MET_TROPO = 3
TDEC_TROP = 259200
TDEC_STRAT = 259200
DT_MOD = 180
DT_MET = 21600   
T_START = ${t0}         
T_STOP = ${t1}
CHUNKSZHINT = 163840000
ATM_DT_OUT = 21600
EOF

cp ${input}/init/pos_glo_16070100.nc ${output}/pos_glo_16070100.nc

echo "Run 6 hours forward calculations with ERA-Interim."
echo ${output} > ${output}/dirlist
$trac/trac ${output}/dirlist trac.ctl pos_glo_16070100.nc\
     ATM_BASENAME atm GRID_BASENAME grid \
     TURB_MESOX 0 TURB_MESOZ 0 \
     TURB_DX_TROP 0 TURB_DZ_STRAT 0 \
>> ${output}/log_trac_diff.txt

error=0

diff -q -s ${input}/diabatic/atm_2016_07_01_00_00.tab ${output}/atm_2016_07_01_00_00.tab
diff -q -s ${input}/diabatic/atm_2016_07_01_06_00.tab ${output}/atm_2016_07_01_06_00.tab

if [ $? -ne 0 ]; then
	error=1
	cat ${output}/log_trac_diff.txt
	exit ${error}
fi

rm ${output}/pos_glo_16070100.nc
rm ${output}/trac.ctl
rm ${output}/log_trac_diff.txt
rm ${output}/dirlist
#rm ${output}/*.tab

exit ${error}
