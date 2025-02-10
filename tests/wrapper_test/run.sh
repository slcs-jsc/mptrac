#! /bin/bash

# Set environment...
ulimit -s unlimited
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4
export LANG=C
export LC_ALL=C

# Setup...
trac=../../src

# Create directories...
rm -rf data && mkdir -p data
mkdir -p data/c
mkdir -p data/fortran

# Test order and missing variables between C and Fortran structure...
./find-vars.sh
[ $? == 99 ] && exit 1

# Set timestep and timerange of simulation...
t0=$($trac/time2jsec 2011 6 5 0 0 0 0)
#t1=$($trac/time2jsec 2011 6 6 0 0 0 0)
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
DIFFUSION = 1
TDEC_TROP = 259200.0
TDEC_STRAT = 259200.0
DRY_DEPO_VDEP = 0.15
DRY_DEPO_DP = 300
MIXING_TROP = 1e-3
MIXING_STRAT = 1e-6
DRY_DEPO_VDEP = 0.15
DRY_DEPO_DP = 300
MIXING_TROP = 1e-3
MIXING_STRAT = 1e-6
DT_MET = 86400.0
T_STOP = $t1
CSI_OBSMIN = 1e-5
CSI_MODMIN = 1e-5
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

# Set initial air parcel positions...
$trac/atm_init data/trac.ctl data/atm_init.tab \
	       INIT_T0 "$t0" INIT_T1 "$t0" \
	       INIT_Z0 10.0 INIT_Z1 10.0 \
	       INIT_LON0 -72.117 INIT_LON1 -72.117 \
	       INIT_LAT0 -40.59 INIT_LAT1 -40.59

# Split air parcels...
$trac/atm_split data/trac.ctl data/atm_init.tab data/atm_split.tab \
		SPLIT_N 10000 SPLIT_M 1e9 SPLIT_DX 30.0 SPLIT_DZ 5.0

# Calculate trajectories with original C code...
echo "C Code..."
echo "data" > data/dirlist
$trac/trac data/dirlist trac.ctl atm_split.tab \
	   ATM_BASENAME atm GRID_BASENAME grid \
	   ENS_BASENAME ens STAT_BASENAME station \
	   CSI_BASENAME csi CSI_OBSFILE data/obs.tab \
	   PROF_BASENAME prof PROF_OBSFILE data/obs.tab \
           SAMPLE_BASENAME sample SAMPLE_OBSFILE data/obs.tab \
	   VTK_BASENAME atm

# trac.c produces the reference data set for the Fortran wrapper
mv data/*00.tab data/*.vtk data/c/

# Calculate trajectories with Fortran Wrapper...
echo "Fortran Wrapper..."
echo "data" > data/dirlist
$trac/trac_fortran data/dirlist trac.ctl atm_split.tab \
	   ATM_BASENAME atm GRID_BASENAME grid \
	   ENS_BASENAME ens STAT_BASENAME station \
	   CSI_BASENAME csi CSI_OBSFILE data/obs.tab \
	   PROF_BASENAME prof PROF_OBSFILE data/obs.tab \
           SAMPLE_BASENAME sample SAMPLE_OBSFILE data/obs.tab \
	   VTK_BASENAME atm

mv data/*00.tab data/*.vtk data/fortran/

# Compare files...
## Files are limited to time, altitude, longitude and latitude, because Fortran wrapper includes only advection 
#echo -e "\nCompare results..."
#error=0
#for f in $(ls data/c/atm*00.tab ) ; do
#    f2=data/fortran/"$(basename "$f")"
#    a=$(awk '{printf $1}' "$f")
#    aa=$(awk '{printf $1}' "$f2")
#    rm -f x y
#    if [ $a == "#" ]; then
#	cut -d' ' -f1-5 $f >> x
#    else
#	cut -d' ' -f1-4 $f >> x
#    fi
#    if [ $aa == "#" ]; then
#	cut -d' ' -f1-5 $f2 >> y
#    else
#	cut -d' ' -f1-4 $f2 >> y
#    fi
#    mv x $f
#    mv y $f2
#    diff -q -s $f $f2 || error=1
#done
#exit $error

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data/c/*.tab data/c/*.vtk) ; do
    diff -q -s data/fortran/"$(basename "$f")" "$f" || error=1
    #    xxdiff data/fortran/"$(basename "$f")" "$f" || error=1
done
exit $error

## Compare files...
#echo -e "\nCompare results..."
#error=0
#for f in $(ls data.ref/*.tab data.ref/*.vtk) ; do
#    diff -q -s data/"$(basename "$f")" "$f" || error=1
#done
#exit $error
