#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4
export LANG=C
export LC_ALL=C

# Setup...
trac=../../src
rm -rf data && mkdir -p data
t0=$($trac/time2jsec 2011 6 5 0 0 0 0)
t1=$($trac/time2jsec 2011 6 5 0 3 0 0)

# Configure a single deposition-only model step...
cat > data/trac.ctl <<EOF
NQ = 6
QNT_NAME[0] = rp
QNT_NAME[1] = rhop
QNT_NAME[2] = m
QNT_NAME[3] = mloss_dry
QNT_NAME[4] = ps
QNT_NAME[5] = t
METBASE = ../data/ei
MET_DT_OUT = 180
ATM_DT_OUT = 180
BUDGET_DT_OUT = 180
DT_MET = 86400
DT_MOD = 180
T_STOP = $t1
DRY_DEPO_DZ = 30
EOF

# Positions are near the surface at this location; one parcel is above the layer.
cat > data/atm_init.tab <<EOF
$t0 0.055 -82.267 12.0748 10 1000 100 0 0 0
$t0 0.065 -82.267 12.0748 0 1000 100 0 0 0
$t0 0.100 -82.267 12.0748 10 1000 100 0 0 0
$t0 0.055 -82.267 12.0748 10 0 100 0 0 0
EOF

echo "data" > data/dirlist
$trac/trac data/dirlist trac.ctl atm_init.tab \
    DRY_DEPO_MODE 1 DRY_DEPO_VDEP 0.01 \
    ATM_BASENAME atm_vdep BUDGET_BASENAME budget_vdep
$trac/trac data/dirlist trac.ctl atm_init.tab \
    DRY_DEPO_MODE 2 DRY_DEPO_VDEP 0 \
    ATM_BASENAME atm_aero BUDGET_BASENAME budget_aero

# Test wet scavenging in and below clouds using ERA5 model-level data...
tw0=$($trac/time2jsec 2011 6 7 0 0 0 0)
tw1=$($trac/time2jsec 2011 6 7 0 3 0 0)
cat > data/wet.ctl <<EOF
NQ = 11
QNT_NAME[0] = m
QNT_NAME[1] = mloss_wet
QNT_NAME[2] = t
QNT_NAME[3] = ps
QNT_NAME[4] = pct
QNT_NAME[5] = pcb
QNT_NAME[6] = cl
QNT_NAME[7] = lwc
QNT_NAME[8] = iwc
QNT_NAME[9] = rwc
QNT_NAME[10] = swc
METBASE = ../data/era5ml
MET_PRESS_LEVEL_DEF = 6
MET_VERT_COORD = 1
MET_DT_OUT = 180
ATM_DT_OUT = 180
BUDGET_DT_OUT = 180
DT_MET = 86400
DT_MOD = 180
T_STOP = $tw1
DRY_DEPO_MODE = 0
RADIO_DEPO = 0
DIFFUSION = 0
CONV_CAPE = 0
TDEC_TROP = 0
TDEC_STRAT = 0
EOF

# In-cloud, below-cloud, and above-cloud-top parcels...
cat > data/wet_init.tab <<EOF
$tw0 5 -96 30 100 0 0 0 0 0 0 0 0 0 0
$tw0 1 42 24 100 0 0 0 0 0 0 0 0 0 0
$tw0 14 -96 30 100 0 0 0 0 0 0 0 0 0 0
EOF

$trac/trac data/dirlist wet.ctl wet_init.tab \
    'WET_DEPO_IC_H[0]' 1000 'WET_DEPO_IC_H[1]' 3000 \
    'WET_DEPO_BC_H[0]' 1000 'WET_DEPO_BC_H[1]' 3000 \
    ATM_BASENAME atm_wet_gas BUDGET_BASENAME budget_wet_gas
$trac/trac data/dirlist wet.ctl wet_init.tab \
    WET_DEPO_IC_A 0.0001 WET_DEPO_IC_B 1 \
    WET_DEPO_BC_A 0.0002 WET_DEPO_BC_B 1 \
    ATM_BASENAME atm_wet_aerosol BUDGET_BASENAME budget_wet_aerosol
$trac/trac data/dirlist wet.ctl wet_init.tab \
    WET_DEPO_IC_A 0.0001 WET_DEPO_IC_B 1 \
    ATM_BASENAME atm_wet_ic_only
$trac/trac data/dirlist wet.ctl wet_init.tab \
    WET_DEPO_BC_A 0.0002 WET_DEPO_BC_B 1 \
    ATM_BASENAME atm_wet_bc_only
$trac/trac data/dirlist wet.ctl wet_init.tab ATM_BASENAME atm_wet_off

# Test radioactive decay, dry deposition, and ground inventories...
cat > data/radio.ctl <<EOF
NQ = 14
QNT_NAME[0] = Arn222
QNT_NAME[1] = Apb210
QNT_NAME[2] = Abe7
QNT_NAME[3] = Acs137
QNT_NAME[4] = Ai131
QNT_NAME[5] = Axe133
QNT_NAME[6] = ps
QNT_NAME[7] = pct
QNT_NAME[8] = cl
QNT_NAME[9] = t
QNT_NAME[10] = lwc
QNT_NAME[11] = iwc
QNT_NAME[12] = rwc
QNT_NAME[13] = swc
METBASE = ../data/ei
MET_DT_OUT = 180
ATM_DT_OUT = 180
DEPO_DT_OUT = 180
DT_MET = 86400
DT_MOD = 180
T_STOP = $t1
DRY_DEPO_MODE = 0
DRY_DEPO_DZ = 30
DIFFUSION = 0
CONV_CAPE = -1
TDEC_TROP = 0
TDEC_STRAT = 0
GRID_LON0 = -83
GRID_LON1 = -82
GRID_NX = 1
GRID_LAT0 = 12
GRID_LAT1 = 13
GRID_NY = 1
EOF
for i in {0..5}; do
    echo "QNT_FORMAT[$i] = %.12g" >> data/radio.ctl
done

# Active, partial-step, unreleased, and Pb-210-only parcels...
tp=$($trac/time2jsec 2011 6 5 0 1 30 0)
tf=$($trac/time2jsec 2011 6 5 0 6 0 0)
cat > data/radio_decay_init.tab <<EOF
$t0 5 -82.267 12.0748 1000 1000 1000 1000 1000 1000 0 0 0 0 0 0 0 0
$tp 5 -82.267 12.0748 1000 1000 1000 1000 1000 1000 0 0 0 0 0 0 0 0
$tf 5 -82.267 12.0748 1000 1000 1000 1000 1000 1000 0 0 0 0 0 0 0 0
$t0 5 -82.267 12.0748 0 1000 1000 1000 1000 1000 0 0 0 0 0 0 0 0
EOF
$trac/trac data/dirlist radio.ctl radio_decay_init.tab \
    RADIO_DECAY 1 RADIO_DEPO 0 ATM_BASENAME atm_radio_decay

# Two parcels deposit into one cell; the third stays above the surface layer.
cat > data/radio_dry_init.tab <<EOF
$t0 0.055 -82.267 12.0748 0 1000 1000 1000 1000 0 0 0 0 0 0 0 0 0
$t0 0.055 -82.267 12.0748 0 500 500 500 500 0 0 0 0 0 0 0 0 0
$t0 0.100 -82.267 12.0748 0 1000 1000 1000 1000 0 0 0 0 0 0 0 0 0
EOF
$trac/trac data/dirlist radio.ctl radio_dry_init.tab \
    RADIO_DECAY 0 RADIO_DEPO 1 ATM_BASENAME atm_radio_dry \
    DEPO_BASENAME depo_radio_dry

# Reuse the verified ERA5 cloud regimes for radioactive wet deposition.
cat > data/radio_wet_init.tab <<EOF
$tw0 5 -96 30 0 1000 1000 1000 1000 0 0 0 0 0 0 0 0 0
$tw0 1 42 24 0 1000 1000 1000 1000 0 0 0 0 0 0 0 0 0
$tw0 14 -96 30 0 1000 1000 1000 1000 0 0 0 0 0 0 0 0 0
EOF
$trac/trac data/dirlist radio.ctl radio_wet_init.tab \
    METBASE ../data/era5ml MET_PRESS_LEVEL_DEF 6 MET_VERT_COORD 1 \
    T_STOP "$tw1" RADIO_DECAY 0 RADIO_DEPO 1 \
    ATM_BASENAME atm_radio_wet DEPO_BASENAME depo_radio_wet \
    GRID_LON0 -100 GRID_LON1 50 GRID_LAT0 20 GRID_LAT1 35

# Two longer steps resolve decay of already deposited I-131 activity.
tc=$($trac/time2jsec 2011 6 5 1 0 0 0)
cat > data/radio_coupled_init.tab <<EOF
$t0 0.055 -82.267 12.0748 0 0 0 0 1000 0 0 0 0 0 0 0 0 0
EOF
$trac/trac data/dirlist radio.ctl radio_coupled_init.tab \
    DT_MOD 1800 MET_DT_OUT 1800 ATM_DT_OUT 1800 DEPO_DT_OUT 1800 \
    T_STOP "$tc" RADIO_DECAY 1 RADIO_DEPO 1 \
    ATM_BASENAME atm_radio_coupled DEPO_BASENAME depo_radio_coupled

# Omit the empty separator line after the single ground-grid cell.
sed -i '${/^$/d;}' data/depo_radio_*.tab

# Compare files...
echo -e "\nCompare results..."
error=0
for f in data.ref/*.tab ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
