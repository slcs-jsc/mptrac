#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4

# Setup...
trac=../../src

# Create directories...
rm -rf data && mkdir -p data

# Set time range of simulation...
t0=$($trac/time2jsec 2011 6 5 0 0 0 0)
t1=$($trac/time2jsec 2011 6 8 0 0 0 0)

# Create control parameter file...
cat > data/trac.ctl <<EOF
# Set quantities...
NQ = 1
QNT_NAME[0] = zg

# Meteorological input data...
METBASE = ../../tests/data/ei
DT_MET = 86400.0

# Sampling frequency of meteo data along the trajectories...
MET_DT_OUT = 86400
EOF

# Set initial air parcel positions...
$trac/atm_init data/trac.ctl data/atm_init.tab \
	       INIT_T0 "$t0" INIT_T1 "$t0" \
	       INIT_Z0 10.0 INIT_Z1 10.0 \
	       INIT_LON0 -72.117 INIT_LON1 -72.117 \
	       INIT_LAT0 -40.59 INIT_LAT1 -40.59

# Split air parcels...
$trac/atm_split data/trac.ctl data/atm_init.tab data/atm_split.tab \
		SPLIT_N 100000 SPLIT_M 1e9 SPLIT_DX 30.0 SPLIT_DZ 1.0

# Calculate trajectories...
echo "data" > data/dirlist
$trac/trac data/dirlist trac.ctl atm_split.tab T_STOP $t1 \
	   VTK_BASENAME atm_plane VTK_SCALE 5
$trac/trac data/dirlist trac.ctl atm_split.tab T_STOP $t1 \
	   VTK_BASENAME atm_sphere VTK_SCALE 100 VTK_SPHERE 1
