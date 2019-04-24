#! /bin/bash

# Setup...
export OMP_NUM_THREADS=4
trac=../src

# Create directories...
rm -rf data plots && mkdir -p data plots

# Set time range of simulations...
t0=$($trac/time2jsec 2011 6 5 0 0 0 0)
t1=$($trac/time2jsec 2011 6 8 0 0 0 0)

# Create control parameter file...
cat > data/trac.ctl <<EOF
NQ = 9
QNT_NAME[0] = t
QNT_NAME[1] = u
QNT_NAME[2] = v
QNT_NAME[3] = w
QNT_NAME[4] = z
QNT_NAME[5] = pv
QNT_NAME[6] = ps
QNT_NAME[7] = pt
QNT_NAME[8] = m
MET_TROPO = 3
MET_GEOPOT = meteo/ei_geopot_1deg.tab
EOF

# Set initial air parcel positions...
$trac/atm_init data/trac.ctl data/atm_init.tab \
	       INIT_T0 $t0 INIT_T1 $t0 \
	       INIT_Z0 10.0 INIT_Z1 10.0 \
	       INIT_LON0 -72.117 INIT_LON1 -72.117 \
	       INIT_LAT0 -40.59 INIT_LAT1 -40.59 \
    | tee data/log_init.tab

# Split air parcels...
$trac/atm_split data/trac.ctl data/atm_init.tab data/atm_split.tab \
		SPLIT_N 10000 SPLIT_M 1e9 SPLIT_DX 30.0 SPLIT_DZ 1.0 \
    | tee data/log_split.tab

# Calculate trajectories (without diffusion)...
echo "data" > data/dirlist
$trac/trac data/dirlist trac.ctl atm_split.tab meteo/ei \
	   DT_MET 259200 T_STOP $t1 DIRECTION 1 \
	   TURB_MESOX 0 TURB_MESOZ 0 \
	   TURB_DX_TROP 0 TURB_DZ_STRAT 0 \
	   ATM_BASENAME atm_nodiff ATM_DT_OUT 43200 \
	   GRID_BASENAME grid_nodiff GRID_DT_OUT 43200 \
	   GRID_LON0 -90 GRID_LON1 60 \
	   GRID_LAT0 -60 GRID_LAT1 -15 \
	   GRID_NX 300 GRID_NY 90 \
    | tee data/log_trac_nodiff.tab

# Calculate trajectories (with diffusion)...
echo "data" > data/dirlist
$trac/trac data/dirlist trac.ctl atm_split.tab meteo/ei \
	   DT_MET 259200 T_STOP $t1 DIRECTION 1 \
	   ATM_BASENAME atm_diff ATM_DT_OUT 43200 \
	   GRID_BASENAME grid_diff GRID_DT_OUT 43200 \
	   GRID_LON0 -90 GRID_LON1 60 \
	   GRID_LAT0 -60 GRID_LAT1 -15 \
	   GRID_NX 300 GRID_NY 90 \
    | tee data/log_trac_diff.tab

# Calculate deviations...
$trac/atm_dist data/trac.ctl data/dist_nodiff.tab absdev $(for f in $(ls data.org/atm_nodiff_2011_*tab) ; do echo data/$(basename $f) $f ; done) \
    | tee data/log_dist_nodiff.tab
$trac/atm_dist data/trac.ctl data/dist_diff.tab absdev $(for f in $(ls data.org/atm_nodiff_2011_*tab) ; do echo data/$(basename $f | sed s/nodiff/diff/g) $f ; done) \
    | tee data/log_dist_diff.tab

# Use gnuplot to plot air parcels...
for f in $(ls data/atm_*_2011*tab) ; do
    echo "Plot $f ..."
    t=$(basename $f .tab | awk 'BEGIN{FS="_"}{print $3"-"$4"-"$5", "$6":"$7" UTC"}')
    gnuplot <<EOF
set out "plots/$(basename $f).png"
set term png truecolor crop linewidth 2 font "Helvetica" 24 size 1440,900
set size ratio 0.75
set pal def
set cbra [5:15]
set cbla "altitude [km]"
set xla "longitude [deg]"
set yla "latitude [deg]"
set xtics 30
set ytics 10
set mxtics 6
set mytics 5
set xra [-90:60]
set yra [-60:-15]
set grid
set title "MPTRAC | $t"
plot "$f" u 3:4:(1.*\$2) w d lc pal z t "", \
    "wcl.tab" u 1:2 w l lt -1 t "", \
    "-" u 1:2 w p pt 9 ps 3 lc rgbcolor "red" t ""
-72.117 -40.59
e
EOF
done

# Use gnuplot to plot grid data...
for f in $(ls data/grid_*_2011*tab) ; do
    echo "Plot $f ..."
    t=$(basename $f .tab | awk 'BEGIN{FS="_"}{print $3"-"$4"-"$5", "$6":"$7" UTC"}')
    gnuplot <<EOF
set out "plots/$(basename $f).png"
set term png truecolor crop linewidth 2 font "Helvetica" 24 size 1440,900
set size ratio 0.75
set pm3d map
set pal def (0 'gray90', 1 'blue', 2 'cyan', 3 'green', 4 'yellow', 5 'orange', 6 'red')
set cbla "column density [g/m^2]"
set cbra [0:10]
set xla "longitude [deg]"
set yla "latitude [deg]"
set xtics 30
set ytics 10
set mxtics 6
set mytics 5
set xra [-90:60]
set yra [-60:-15]
set grid
set title "MPTRAC | $t"
splot "$f" u 3:4:(1e3*\$8) t "", \
    "wcl.tab" u 1:2:(0) w l lt -1 t "", \
    "-" u 1:2:3 w p pt 9 ps 3 lc rgbcolor "red" t ""
-72.117 -40.59 0
e
EOF
done

# Check deviations...
echo -e "\nDeviations from reference output (without diffusion):"
diff -s data/dist_nodiff.tab data.org/dist_nodiff.tab
echo -e "\nDeviations from reference output (with diffusion):"
diff -s data/dist_diff.tab data.org/dist_diff.tab

# Show timers...
echo -e "\nTimers, memory usage, and problem size:"
grep TIMER data/log_trac_diff.tab
grep MEMORY data/log_trac_diff.tab
grep SIZE data/log_trac_diff.tab
