#! /bin/bash

# Testing components for better interoperable use between CLaMS and MPTRAC...

# (1) Test reading and writing CLaMS position files...

# Create directories...
rm -rf data && mkdir -p data

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH

# Path to log files...
log=./data

# Pathway to the application...
atm_conv=../../src/atm_conv

# Path to the atmosphere data that containts air parcel positions...
atm_in=./data.ref/atm_input.tab

# Path to the output atmosphere data...
atm_out_nc=./data/atm_output.nc

# Run the application...
echo "Test reading atm_type 0 data and write it to atm_type 4 data."
echo ${atm_conv} - ${atm_in} 0 ${atm_out_nc} 4
${atm_conv} - ${atm_in} 0 ${atm_out_nc} 4 >> ${log}/log_0to4.txt

# Path to the atmosphere data that containts air parcel positions...
atm_in=${atm_out_nc}

# Path to the output atmosphere data...
atm_out_tab=./data/atm_output.tab

# Path to the reference...
atm_ref_tab=./data.ref/atm_output.tab

# Run the application...
echo "Reading atmospheric data of atm_type 4 and write it as atm_type 0."
echo ${atm_conv} - ${atm_in} 4 ${atm_out_tab} 0
${atm_conv} - ${atm_in} 4 ${atm_out_tab} 0 > ${log}/log_4to0.txt

# Compare files...
error=0
diff -q -s ${atm_out_tab} ${atm_ref_tab}

if [ $? -ne 0 ]; then
    error=1
    cat ${log}/log_0to4.txt
    cat ${log}/log_4to0.txt
    exit $error
fi

rm ${atm_out_tab} ${atm_out_nc} ${log}/*.txt

# (2) Test if the diabatic transport calculations have been changed...

# directory of the trajectory code...
trac=../../src

# directory of the reanalysis data...
metbase=../data/clams/erai_vlr

# directory for the output...
output=./data

# directory of the reference and init...
input=./data.ref/

# start and end date
t0=$(${trac}/time2jsec 2016 7 1 0 0 0 0)
t1=$(${trac}/time2jsec 2016 7 1 6 0 0 0)

# Create control parameter file...
cat > ${output}/trac.ctl <<EOF
MET_CONVENTION = 1
MET_PRESS_LEVEL_DEF = 2
ATM_TYPE = 3
ATM_TYPE_OUT = 0
ADVECT = 2
MET_CLAMS = 1
ADVECT_VERT_COORD = 1
MET_VERT_COORD = 1
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
    | tee ${output}/log_trac_diff.txt

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
rm ${output}/*.tab

exit ${error}
