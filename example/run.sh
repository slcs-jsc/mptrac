#! /bin/bash

# Slurm configuration for JUWELS-Cluster...
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=00:05:00
#SBATCH --account=slmet
#SBATCH --partition=batch

## Slurm configuration for JUWELS-Booster...
##SBATCH --nodes=1
##SBATCH --ntasks=1
##SBATCH --ntasks-per-node=1
##SBATCH --cpus-per-task=48
##SBATCH --time=00:05:00
##SBATCH --account=slmet
##SBATCH --partition=booster
##SBATCH --gres=gpu:4

# Load modules (uncomment for JUWELS)...
#ml purge
#ml gnuplot
#ml GCC ParaStationMPI     # for MPI runs
#ml NVHPC ParaStationMPI   # for GPU runs

# Check arguments...
if [ $# -ne 4 ] ; then
    echo "usage: $0 <meteo> <dtmet> <np> <nthreads>"
    exit
fi

# Set environment...
export LD_LIBRARY_PATH=../libs/build/lib:$LD_LIBRARY_PATH
[ $4 -gt 0 ] && export OMP_NUM_THREADS=$4

# Set variables...
trac=../src
meteo=$1
dtmet=$2
np=$3

# Download meteo data...
metdir=meteo_${meteo}_${dtmet}s
if [ ! -d $metdir ] ; then
    echo "Downloading meteo data (this may take a while)..."
    mkdir -p $metdir && cd $metdir || exit
    wget --progress=dot:giga https://datapub.fz-juelich.de/slcs/mptrac/data/test_case/${meteo}_${dtmet}s.zip \
	&& unzip ${meteo}_${dtmet}s.zip && rm ${meteo}_${dtmet}s.zip || exit
    cd -
fi
metbase=$(basename $(ls $metdir/*.nc | head -1) | awk 'BEGIN{FS="_"}{print $1}')

# Create directories...
rm -rf data plots && mkdir -p data plots

# Set timestep and timerange of simulation...
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
METBASE = $metdir/$metbase
TDEC_TROP = 259200
TDEC_STRAT = 259200
DT_MET = $dtmet
T_STOP = $t1
GRID_LON0 = -90
GRID_LON1 = 60
GRID_LAT0 = -60
GRID_LAT1 = -15
GRID_NX = 300
GRID_NY = 90
EOF

# Set initial air parcel positions...
$trac/atm_init data/trac.ctl data/atm_init.tab \
	       INIT_T0 $t0 INIT_T1 $t0 \
	       INIT_Z0 10.0 INIT_Z1 10.0 \
	       INIT_LON0 -72.117 INIT_LON1 -72.117 \
	       INIT_LAT0 -40.59 INIT_LAT1 -40.59

# Split air parcels...
$trac/atm_split data/trac.ctl data/atm_init.tab data/atm_split.tab \
		SPLIT_N $np SPLIT_M 1e9 SPLIT_DX 30.0 SPLIT_DZ 1.0

# Calculate trajectories...
echo "data" > data/dirlist
$trac/trac data/dirlist trac.ctl atm_split.tab \
	   ATM_BASENAME atm GRID_BASENAME grid

# Plot air parcel data...
echo
for f in $(ls data/atm_2011*tab) ; do
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
    "$metdir/wcl.tab" u 1:2 w l lt -1 t "", \
    "-" u 1:2 w p pt 9 ps 3 lc rgbcolor "red" t ""
-72.117 -40.59
e
EOF
done

# Plot grid data...
for f in $(ls data/grid_2011*tab) ; do
    echo "Plot $f ..."
    t=$(basename $f .tab | awk 'BEGIN{FS="_"}{print $3"-"$4"-"$5", "$6":"$7" UTC"}')
    gnuplot <<EOF
set out "plots/$(basename $f).png"
set term png truecolor crop linewidth 2 font "Helvetica" 24 size 1440,900
set size ratio 0.75
set pm3d map interp 4,4
set pal def (0 'gray90', 1 'blue', 2 'cyan', 3 'green', 4 'yellow', 5 'orange', 6 'red')
set cbla "column density [g/m^2]"
set cbra [0:1]
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
    "$metdir/wcl.tab" u 1:2:(0) w l lt -1 t "", \
    "-" u 1:2:3 w p pt 9 ps 3 lc rgbcolor "red" t ""
-72.117 -40.59 0
e
EOF
done

# Compare atm and grid file...
echo
error=0
for f in $(ls data.ref/atm*tab data.ref/grid*tab) ; do
    diff -q -s data/$(basename $f) $f
    [ $? -ne 0 ] && error=1
done
exit $error
