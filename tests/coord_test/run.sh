#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4

# Setup...
trac=../../src
set -e errexit

# Create directories...
rm -rf data && mkdir -p data

# Set timestep and timerange of simulation...
t0=$($trac/time2jsec 2025 5 1 0 0 0 0)
t1=$($trac/time2jsec 2025 5 1 2 0 0 0)

cat > data/trac.ctl <<EOF
NQ = 4
QNT_NAME[0] = t
QNT_NAME[1] = u
QNT_NAME[2] = v
QNT_NAME[3] = w
METBASE = ../data/era5_utm32
TRACER_CHEM = 0
DIFFUSION = 1
DT_MET = 3600.0
T_STOP = $t1
EOF

# Set initial air parcel positions...
$trac/atm_init data/trac.ctl data/atm_init.tab \
	       INIT_T0 "$t0" INIT_T1 "$t0" \
	       INIT_Z0 2.0 INIT_Z1 2.0 \
	       INIT_LON0 691090 INIT_LON1 691090 \
         INIT_LAT0 5336247 INIT_LAT1 5336247
#         INIT_LON0 11.5692782 INIT_LON1 11.5692782 \
#         INIT_LAT0 48.1507476 INIT_LAT1 48.1507476

# Split air parcels...
$trac/atm_split data/trac.ctl data/atm_init.tab data/atm_split.tab \
		SPLIT_N 1000 SPLIT_M 1e9 SPLIT_DX 1.0 SPLIT_DZ 1.0

echo "data" > data/dirlist
$trac/trac data/dirlist trac.ctl atm_split.tab \
    ATM_BASENAME atm MET_CAPE 0 DT_MOD 600 ATM_DT_OUT 600 MET_COORD_TYPE 1

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error