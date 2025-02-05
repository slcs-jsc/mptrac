#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4
export LANG=C
export LC_ALL=C

# Setup...
trac=../../src

# Create directories...
rm -rf data && mkdir -p data

# Set timestep and timerange of simulation...
t0=$($trac/time2jsec 2011 6 5 0 0 0 0)
t1=$($trac/time2jsec 2011 6 8 0 0 0 0)

# Create control parameter file...
cat > data/trac.ctl <<EOF
NQ = 12
QNT_NAME[0] = m
QNT_NAME[1] = Cx 
QNT_NAME[2] = Co1d 
QNT_NAME[3] = Ch 
QNT_NAME[4] = Co3 
QNT_NAME[5] = Co3p 
QNT_NAME[6] = Cho2 
QNT_NAME[7] = Coh 
QNT_NAME[8] = Ch2o2
QNT_NAME[9] = Cho2
QNT_NAME[10] = mloss_kpp
QNT_NAME[11] = mloss_wet
METBASE = ../data/ei
MET_DT_OUT = 86400.0
KPP_CHEM = 1
DT_KPP = 7200
MIXING_DT = 7200
MIXING_TROP = 1
MIXING_STRAT = 1
CHEMGRID_LON0 = -90
CHEMGRID_LON1 = 60
CHEMGRID_LAT0 = -60
CHEMGRID_LAT1 = -15
CHEMGRID_Z0 = 0
CHEMGRID_Z1 = 20
CHEMGRID_NZ = 20
CHEMGRID_NX = 150
CHEMGRID_NY = 45
MOLMASS = 64
CONV_CAPE = 0.0
TDEC_TROP = 259200.0
TDEC_STRAT = 259200.0
DT_MET = 86400.0
T_STOP = $t1
MIXING_TROP = 1e-3
MIXING_STRAT = 1e-6
DIFFUSION = 1
RNG_TYPE = 0
EOF

# Create observation file...
echo | awk -v tobs=$($trac/time2jsec 2011 6 7 0 0 0 0) '{
  for(lon=-25; lon<=-15; lon+=0.5)
    for(lat=-50; lat<=-25; lat+=1) {
      if(lon>=-24 && lon<=-21 && lat>=-36 && lat<=-28)
        obs=0.005
      else
        obs=0
      printf("%.2f %g %g %g %g\n", tobs, 0, lon, lat, obs)
    }
}' > data/obs.tab

# Set initial air parcel positions...
$trac/atm_init data/trac.ctl data/atm_init.tab \
	       INIT_T0 "$t0" INIT_T1 "$t0" \
	       INIT_Z0 10.0 INIT_Z1 10.0 \
	       INIT_LON0 -72.117 INIT_LON1 -72.117 \
	       INIT_LAT0 -40.59 INIT_LAT1 -40.59

# Split air parcels...
$trac/atm_split data/trac.ctl data/atm_init.tab data/atm_split.tab \
		SPLIT_N 100 SPLIT_M 1e9 SPLIT_DX 30.0 SPLIT_DZ 1.0

# Calculate trajectories...
echo "data" > data/dirlist
$trac/trac data/dirlist trac.ctl atm_split.tab \
	   ATM_BASENAME atm 

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/$(basename "$f") "$f"
    [ $? -ne 0 ] && error=1
done
exit $error
