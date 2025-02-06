#! /bin/bash

# Testing components for interoperable use between CLaMS and MPTRAC...

# ------------------------------------------------------------
# (1) Test reading and writing CLaMS position files...
# ------------------------------------------------------------

# Setup...
trac=../../src

# Create directories...
rm -rf data && mkdir -p data

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH

# Run the application...
$trac/atm_conv - ./data.ref/atm_input.tab 0 ./data/atm_output.nc 4
$trac/atm_conv - ./data/atm_output.nc 4 ./data/atm_output.tab 0

# ------------------------------------------------------------
# (2) Test diabatic transport calculations...
# ------------------------------------------------------------

# start and end date
t0=$(${trac}/time2jsec 2016 7 1 0 0 0 0)
t1=$(${trac}/time2jsec 2016 7 1 6 0 0 0)

# Create control parameter file...
cat > ./data/trac.ctl <<EOF
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
METBASE = ../data/clams/erai_vlr
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

# Calculate trajectories...
echo "./data" > ./data/dirlist
$trac/trac ./data/dirlist trac.ctl ../data.ref/init/pos_glo_16070100.nc \
	   ATM_BASENAME atm GRID_BASENAME grid \
	   TURB_MESOX 0 TURB_MESOZ 0 \
	   TURB_DX_TROP 0 TURB_DZ_STRAT 0

# Compare files...
diff -q -s ./data.ref/atm_output.tab ./data/atm_output.tab || error=1
diff -q -s ./data.ref/diabatic/atm_2016_07_01_00_00.tab ./data/atm_2016_07_01_00_00.tab || error=1
diff -q -s ./data.ref/diabatic/atm_2016_07_01_06_00.tab ./data/atm_2016_07_01_06_00.tab || error=1
exit ${error}
