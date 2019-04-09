#! /bin/bash

# Create directories...
rm -rf data plots && mkdir -p data plots

# Set time range of simulations...
t0=$(../src/time2jsec 2011 6 5 0 0 0 0)
t1=$(../src/time2jsec 2011 6 5 1 0 0 0)

# Add "mass" to air parcel quantities...
qnt="NQ 1 QNT_NAME[0] m"

# Set initial air parcel positions at 10 km height above Puyehue volcano...
../src/atm_init - data/atm_init.tab $qnt \
		INIT_T0 $t0 INIT_T1 $t0 \
		INIT_Z0 10.0 INIT_Z1 10.0 \
		INIT_LON0 -72.117 INIT_LON1 -72.117 \
		INIT_LAT0 -40.59 INIT_LAT1 -40.59

# Split into 10000 air parcels with a total mass of 10^9 kg
# and scatter around initial position...
../src/atm_split - data/atm_init.tab data/atm_split.tab $qnt \
		 SPLIT_N 10000 SPLIT_M 1e9 SPLIT_DX 30.0 SPLIT_DZ 1.0

# Calculate trajectories, write air parcel and grid output every 12 h...
echo "data" > data/dirlist
../src/trac data/dirlist - atm_split.tab meteo/ei DT_MET 259200 $qnt \
	    T_STOP $t1 DIRECTION 1 \
	    ATM_BASENAME atm ATM_DT_OUT 600 \
	    GRID_BASENAME grid GRID_DT_OUT 600 \
	    GRID_LON0 -100 GRID_LON1 40 \
	    GRID_LAT0 -65 GRID_LAT1 -20 \
	    GRID_NX 280 GRID_NY 90 \
            TURB_DX_TROP 0 \
            TURB_DX_STRAT 0 \
            TURB_DZ_TROP 0 \
            TURB_DZ_STRAT 0 \
            TURB_MESOX 0 \
            TURB_MESOZ 0 
#../src/trac data/dirlist - atm_split.tab meteo/ei DT_MET 259200 $qnt \
#	    T_STOP $t1 DIRECTION 1 \
#	    ATM_BASENAME atm ATM_DT_OUT 43200 \
#	    GRID_BASENAME grid GRID_DT_OUT 43200 \
#	    GRID_LON0 -100 GRID_LON1 40 \
#	    GRID_LAT0 -65 GRID_LAT1 -20 \
#	    GRID_NX 280 GRID_NY 90

# Use gnuplot to plot air parcels...
for f in $(ls data/atm_2011*tab) ; do
    echo "Plot $f ..."
    t=$(basename $f .tab | awk 'BEGIN{FS="_"}{print $2"-"$3"-"$4", "$5":"$6" UTC"}')
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
set xra [-100:40]
set yra [-65:-20]
set grid
set title "MPTRAC | $t"
plot "$f" u 3:4:(1.*\$2) w d lc palette z t "", \
    "wcl.tab" u 1:2 w l lt -1 t "", "-" u 1:2 w p pt 9 ps 3 lc rgbcolor "red" t ""
-72.117 -40.59
e
EOF
done

# Use gnuplot to plot grid data...
for f in $(ls data/grid*tab) ; do
    echo "Plot $f ..."
    t=$(basename $f .tab | awk 'BEGIN{FS="_"}{print $2"-"$3"-"$4", "$5":"$6" UTC"}')
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
set xra [-100:40]
set yra [-65:-20]
set grid
set title "MPTRAC | $t"
splot "$f" u 3:4:(1e3*\$8) t "", \
    "wcl.tab" u 1:2:(0) w l lt -1 t "", \
    "-" u 1:2:3 w p pt 9 ps 3 lc rgbcolor "red" t ""
-72.117 -40.59 0
e
EOF
done

# Remove data directory...
#rm -rf data
