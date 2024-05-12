#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4
export LANG=C
export LC_ALL=C

# Setup...
trac=../../src

# Create directories...
rm -rf data plots && mkdir -p data plots

# Set time range of simulation...
t0=$($trac/time2jsec 2011 6 5 0 0 0 0)
t1=$($trac/time2jsec 2011 6 8 0 0 0 0)

# Create control parameter file...
cat > data/trac.ctl <<EOF
# Set quantities...
NQ = 9
QNT_NAME[0] = t
QNT_NAME[1] = u
QNT_NAME[2] = v
QNT_NAME[3] = w
QNT_NAME[4] = zg
QNT_NAME[5] = pv
QNT_NAME[6] = ps
QNT_NAME[7] = pt
QNT_NAME[8] = m

# Meteorological input data...
METBASE = ../../tests/data/ei
DT_MET = 86400.0

# Activate air parcel output...
ATM_BASENAME = atm

# Activate grid output...
GRID_BASENAME = grid
GRID_LON0 = -90
GRID_LON1 = 60
GRID_NX = 300
GRID_LAT0 = -65
GRID_LAT1 = -15
GRID_NY = 100
EOF

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
$trac/trac data/dirlist trac.ctl atm_split.tab T_STOP "$t1"

# Plot air parcel data...
echo
for f in $(ls data/atm_2011*tab) ; do
    echo "Plot $f ..."
    t=$(basename "$f" .tab | awk 'BEGIN{FS="_"}{print $2"-"$3"-"$4", "$5":"$6" UTC"}')
    gnuplot <<EOF
set out "plots/$(basename "$f").png"
set term pngcairo truecolor crop linewidth 2 font "Helvetica,24" size 1440,900
set size ratio 0.75
set grid lw 0.75

set cbla "Altitude [km]"
set cbra [5:15]
set pal def

set xtics geographic
set format x "%D°%E"
set xra [-90:60]
set xtics 30
set mxtics 6

set ytics geographic
set format y "%D°%N"
set yra [-65:-15]
set ytics 10
set mytics 5

set title "MPTRAC | $t"
plot "$f" u 3:4:(1.*\$2) w d lc pal z t "", \
    "../../data/wcl.tab" u 1:2 w l lt -1 t "", \
    "-" u 1:2 w p pt 9 ps 3 lc rgbcolor "red" t ""
-72.117 -40.59
e
EOF
done

# Plot grid data...
for f in $(ls data/grid_2011*tab) ; do
    echo "Plot $f ..."
    t=$(basename "$f" .tab | awk 'BEGIN{FS="_"}{print $2"-"$3"-"$4", "$5":"$6" UTC"}')
    gnuplot <<EOF
set out "plots/$(basename "$f").png"
set term pngcairo truecolor crop linewidth 2 font "Helvetica,24" size 1440,900
set size ratio 0.75
set grid lw 0.75

set pm3d map interp 4,4
set pal def (0 'gray90', 1 'blue', 2 'cyan', 3 'green', 4 'yellow', 5 'orange', 6 'red')
set cbla "Column density [g/m^2]"
set cbra [0:1]

set xtics geographic
set format x "%D°%E"
set xra [-90:60]
set xtics 30
set mxtics 6

set ytics geographic
set format y "%D°%N"
set yra [-65:-15]
set ytics 10
set mytics 5

set title "MPTRAC | $t"
splot "$f" u 3:4:(1e3*\$7) t "", \
    "../../data/wcl.tab" u 1:2:(0) w l lt -1 t "", \
    "-" u 1:2:3 w p pt 9 ps 3 lc rgbcolor "red" t ""
-72.117 -40.59 0
e
EOF
done

# Compare files...
echo
error=0
for f in $(ls data.ref/atm*tab data.ref/grid*tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
