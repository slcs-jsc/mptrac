#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4

# Setup...
trac=../../src

# Create directories...
rm -rf data && mkdir -p data

# Set timestep and timerange of simulation...
t0=$($trac/time2jsec 2011 6 5 0 0 0 0)
t1=$($trac/time2jsec 2011 6 8 0 0 0 0)

# Create control parameter file...
cat > data/trac.ctl <<EOF
NQ = 10
QNT_NAME[0] = t
QNT_NAME[1] = u
QNT_NAME[2] = v
QNT_NAME[3] = w
QNT_NAME[4] = z
QNT_NAME[5] = pv
QNT_NAME[6] = ps
QNT_NAME[7] = pt
QNT_NAME[8] = m
QNT_NAME[9] = stat
METBASE = ../data/ei
MET_DT_OUT = 86400.0
SPECIES = SO2
BOUND_MASS = 0.0
BOUND_DPS = 100.0
CONV_CAPE = 0.0
TDEC_TROP = 259200.0
TDEC_STRAT = 259200.0
DT_MET = 86400.0
T_STOP = $t1
GRID_LON0 = -90
GRID_LON1 = 60
GRID_LAT0 = -60
GRID_LAT1 = -15
GRID_NX = 300
GRID_NY = 90
SAMPLE_DZ = 100
STAT_LON = -22
STAT_LAT = -40
EOF

# Create observation file...
echo | awk -v tobs=$($trac/time2jsec 2011 6 7 0 0 0 0) '{
  for(lon=-25; lon<=-15; lon+=0.5)
    for(lat=-50; lat<=-25; lat+=1)
      printf("%.2f %g %g %g %g\n", tobs, 0, lon, lat, 0)
}' > data/obs.tab

# Set initial air parcel positions...
$trac/atm_init data/trac.ctl data/atm_init.tab \
	       INIT_T0 "$t0" INIT_T1 "$t0" \
	       INIT_Z0 10.0 INIT_Z1 10.0 \
	       INIT_LON0 -72.117 INIT_LON1 -72.117 \
	       INIT_LAT0 -40.59 INIT_LAT1 -40.59

# Split air parcels...
$trac/atm_split data/trac.ctl data/atm_init.tab data/atm_split.tab \
		SPLIT_N 10000 SPLIT_M 1e9 SPLIT_DX 30.0 SPLIT_DZ 1.0

# Calculate trajectories...
echo "data" > data/dirlist
$trac/trac data/dirlist trac.ctl atm_split.tab \
	   ATM_BASENAME atm GRID_BASENAME grid STAT_BASENAME station \
           SAMPLE_BASENAME sample SAMPLE_OBSFILE data/obs.tab \
	   PROF_BASENAME prof PROF_OBSFILE data/obs.tab

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/$(basename "$f") "$f"
    [ $? -ne 0 ] && error=1
done
exit $error
