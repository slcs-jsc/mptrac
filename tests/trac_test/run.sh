#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4
export LANG=C
export LC_ALL=C
export GSL_RNG_SEED=123

# Setup...
trac=../../src

# Create directories...
rm -rf data && mkdir -p data

# Set timestep and timerange of simulation...
t0=$($trac/time2jsec 2011 6 5 0 0 0 0)
t1=$($trac/time2jsec 2011 6 8 0 0 0 0)

# Create control parameter file...
cat > data/trac.ctl <<EOF
NQ = 13
QNT_NAME[0] = t
QNT_NAME[1] = u
QNT_NAME[2] = v
QNT_NAME[3] = w
QNT_NAME[4] = zg
QNT_NAME[5] = pv
QNT_NAME[6] = ps
QNT_NAME[7] = pt
QNT_NAME[8] = m
QNT_NAME[9] = stat
QNT_NAME[10] = ens
QNT_NAME[11] = Cccl3f
QNT_NAME[12] = Cx
METBASE = ../data/ei
MET_DT_OUT = 86400.0
SPECIES = SO2
BOUND_LAT0 = -90
BOUND_LAT1 = 90
BOUND_P0 = 1e10
BOUND_P1 = -1e10
BOUND_DPS = 100.0
BOUND_MASS = 0.0
CONV_CAPE = 0.0
H2O2_CHEM_REACTION = 1
TRACER_CHEM = 1
CHEMGRID_NX = 72
CHEMGRID_NY = 36
CHEMGRID_NZ = 30
DIFFUSION = 1
TDEC_TROP = 259200.0
TDEC_STRAT = 259200.0
DRY_DEPO_VDEP = 0.15
DRY_DEPO_DP = 300
MIXING_TROP = 1e-3
MIXING_STRAT = 1e-6
DT_MET = 86400.0
DT_MOD = 300.0
T_STOP = $t1
CSI_OBSMIN = 1e-5
CSI_MODMIN = 1e-5
GRID_LON0 = -90
GRID_LON1 = 60
GRID_LAT0 = -60
GRID_LAT1 = -15
GRID_NX = 300
GRID_NY = 90
GRID_SPARSE = 1
SAMPLE_DZ = 100
STAT_LON = -22
STAT_LAT = -40
VTK_STRIDE = 100
EOF

# Create observation file...
echo | awk -v tobs="$($trac/time2jsec 2011 6 7 0 0 0 0)" '{
  for(lon=-25; lon<=-15; lon+=0.5)
    for(lat=-50; lat<=-25; lat+=1) {
      if(lon>=-24 && lon<=-21 && lat>=-36 && lat<=-28)
        obs=0.005
      else
        obs=0
      printf("%.2f %g %g %g %g\n", tobs, 0, lon, lat, obs)
    }
}' > data/obs.tab

# Initialize globally distributed, well-mixed air parcels...
$trac/atm_init data/trac.ctl data/atm_init.tab \
	       INIT_T0 "$t0" INIT_T1 "$t0" \
	       INIT_Z0 0.0 INIT_Z1 0.0 \
	       INIT_LON0 0.0 INIT_LON1 0.0 INIT_ULON 360.0 \
	       INIT_LAT0 0.0 INIT_LAT1 0.0 INIT_ULAT 180.0 \
	       INIT_COL_MASS 1 INIT_WELL_MIXED 1 INIT_EVENLY 1 \
	       INIT_NP 10000 INIT_MASS 1e9

# Calculate trajectories on model levels...
echo "data" > data/dirlist
$trac/trac data/dirlist trac.ctl atm_init.tab \
	   METBASE ../data/era5ml MET_PRESS_LEVEL_DEF 6 \
	   MET_VERT_COORD 1 ADVECT_VERT_COORD 2 \
	   ATM_BASENAME atm_ml GRID_BASENAME grid_ml \
	   ENS_BASENAME ens_ml STAT_BASENAME station_ml \
	   CSI_BASENAME csi_ml CSI_OBSFILE data/obs.tab \
	   PROF_BASENAME prof_ml PROF_OBSFILE data/obs.tab \
           SAMPLE_BASENAME sample_ml SAMPLE_OBSFILE data/obs.tab \
	   VTK_BASENAME atm_ml

# Calculate trajectories on pressure levels...
$trac/trac data/dirlist trac.ctl atm_init.tab \
	   ATM_BASENAME atm_pl GRID_BASENAME grid_pl \
	   ENS_BASENAME ens_pl STAT_BASENAME station_pl \
	   CSI_BASENAME csi_pl CSI_OBSFILE data/obs.tab \
	   PROF_BASENAME prof_pl PROF_OBSFILE data/obs.tab \
           SAMPLE_BASENAME sample_pl SAMPLE_OBSFILE data/obs.tab \
	   VTK_BASENAME atm_pl

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab data.ref/*.vtk) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
