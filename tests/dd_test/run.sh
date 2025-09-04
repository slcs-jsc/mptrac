#! /bin/bash
set -e

# Usage: ./run.sh [skip_compile|skip] [hpc]
# Arguments:
#   skip_compile|skip: Skip compilation step
#   hpc: Run in HPC mode (uses ml and srun), otherwise uses mpirun
# Arguments can be provided in any order

#################################################################################################################
# Setup directory structure, OMP, modules, library paths
#################################################################################################################

# Get absolute path to mptrac directory
mptrac_dir=$(cd ../../ && pwd)
work_dir=$mptrac_dir/tests/dd_test
echo "[INFO] MPTRAC directory: $mptrac_dir"
echo "[INFO] Working directory: $work_dir"
cd $work_dir
rm -rf data
mkdir -p data/wind_data

# Parse arguments - check all arguments for keywords
skip_compile=""
hpc_mode=""
for arg in "$@"; do
    case "$arg" in
        "skip"|"skip_compile")
            skip_compile="$arg"
            ;;
        "hpc")
            hpc_mode="hpc"
            ;;
        *)
            echo "[WARNING] Unknown argument: $arg"
            ;;
    esac
done

# Check if running on HPC system
if [[ "$hpc_mode" == "hpc" ]]; then
    echo "[INFO] Running in HPC mode"
    ml purge
    ml NVHPC ParaStationMPI bzip2/.1.0.8
else
    echo "[INFO] Running in non-HPC mode"
fi

# Set environment
export LD_LIBRARY_PATH=$mptrac_dir/libs/build/lib/:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=true
export OMP_PLACES=sockets
#export LANG=C
#export LC_ALL=C

#################################################################################################################
# Define simulation parameters
#################################################################################################################

# Set number of domains in each direction
domains_lon=3                       # number of subdomains in longitudinal direction
domains_lat=3                       # number of subdomains in latitudinal direction

# Set parameters for atmospheric initialization
particles_subdomain_lon_stride=4    # number of particles in longitudinal direction per subdomain
particles_subdomain_lat_stride=4    # number of particles in latitudinal direction per subdomain
grid_lon_start=0                    # starting longitude of the grid for the atmospheric initialization
grid_lon_d=1                        # longitude direction (1 or -1)
grid_lon_size=360                   # size of the grid for the atmospheric initialization in longitudinal direction
grid_lat_start=-90                  # starting latitude of the grid for the atmospheric initialization
grid_lat_d=1                        # latitude direction (1 or -1)
grid_lat_size=180                   # size of the grid for the atmospheric initialization in latitudinal direction

# Set meteo file parameters
metbase=$work_dir/data/wind_data/wind    # path to meteo files
met_press_level_def=-1  # pressure level definition
met_clams=0             # follow CLaMS filename guideline
met_vert_coord=0        # original vertical coordinate type
dt_met=3600             # read time interval for meteo files
dt_mod=600              # model time step
met_dt_out=3600         # output time interval for meteo files
atm_type=3              # original atmospheric data type
met_grid_nx=360         # number of grid points in longitudinal direction of meteo files
met_grid_ny=181         # number of grid points in latitudinal direction of meteo files
met_grid_nz=60          # number of vertical levels of meteo files
dex_glob=364            # allocated size of meteo grid in longitudinal direction
dey_glob=186            # allocated size of meteo grid in latitudinal direction
dez_glob=90             # allocated size of meteo grid in vertical direction
wind_alpha=10           # angle for wind rotation  in degrees
wind_speed=38.5876601   # wind speed in m/s

# Calculate derived parameters
processes=$((domains_lon * domains_lat))
particles_num_subdomain=$((particles_subdomain_lon_stride * particles_subdomain_lat_stride))
npmax=$((domains_lon * domains_lat * particles_num_subdomain * 2))
subgrid_lon_size=$(echo "scale=6; $grid_lon_size / $domains_lon" | bc)
subgrid_lat_size=$(echo "scale=6; $grid_lat_size / $domains_lat" | bc)
subgrid_pos_dlon=$(echo "scale=6; $subgrid_lon_size / $particles_subdomain_lon_stride" | bc)
subgrid_pos_dlat=$(echo "scale=6; $subgrid_lat_size / $particles_subdomain_lat_stride" | bc)
subgrid_pos_lon_offset=$(echo "scale=6; $subgrid_pos_dlon / 2" | bc)
subgrid_pos_lat_offset=$(echo "scale=6; $subgrid_pos_dlat / 2" | bc)
dex=$((met_grid_nx/domains_lon + 10))
dey=$((met_grid_ny/domains_lat + 10))
grid_nx=$((met_grid_nx * grid_lon_size / 360))
grid_ny=$((met_grid_ny * grid_lat_size / 180))

#print parameters
echo "[INFO] processes: $processes"
echo "[INFO] particles_num_subdomain: $particles_num_subdomain"
echo "[INFO] npmax: $npmax"
echo "[INFO] subgrid_lon_size: $subgrid_lon_size"
echo "[INFO] subgrid_lat_size: $subgrid_lat_size"
echo "[INFO] subgrid_pos_dlon: $subgrid_pos_dlon"
echo "[INFO] subgrid_pos_dlat: $subgrid_pos_dlat"
echo "[INFO] subgrid_pos_lon_offset: $subgrid_pos_lon_offset"
echo "[INFO] subgrid_pos_lat_offset: $subgrid_pos_lat_offset"
echo "[INFO] dex: $dex"
echo "[INFO] dey: $dey"


#################################################################################################################
# Compile MPTRAC with domain decomposition
#################################################################################################################

defs="-DNP=${npmax} -DEX_GLOB=${dex_glob} -DEY_GLOB=${dey_glob} -DEP=${dez_glob} -DEX=${dex} -DEY=${dey} -DNQ=5"

echo "[INFO] defs: $defs"
echo "[INFO] LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

# Check if compilation should be skipped
if [[ "$skip_compile" == "skip" || "$skip_compile" == "skip_compile" ]]; then
    echo "[INFO] Skipping compilation (argument: $skip_compile)"
else
    echo "[INFO] Compile MPTRAC with domain decomposition"
    cd $mptrac_dir/src
    make clean
    make -j8 DD=1 MPI=1 STATIC=0 GPU=0 THRUST=1 COMPILER=nvc DEFINES="$defs" || exit 1
    cd $work_dir
fi

#################################################################################################################
# Create atmospheric initialization files & control files for each subdomain
#################################################################################################################

# Set timerange of simulation
t0=$($mptrac_dir/src/time2jsec 2022 6 2 0 0 0 0)
t1=$($mptrac_dir/src/time2jsec 2022 6 2 6 0 0 0)

echo "[INFO] write dirlist.tab and config.ctl to data directory"
> data/dirlist.tab  # Clear the file first
for i in $(seq 0 $((processes-1))); do
    echo "data.$i" >> data/dirlist.tab
done

cat > data/config.ctl <<EOF
MET_PRESS_LEVEL_DEF = ${met_press_level_def}
MET_CLAMS = ${met_clams}
MET_VERT_COORD = ${met_vert_coord}
MET_DT_OUT = ${met_dt_out}
ATM_TYPE = ${atm_type}
ATM_TYPE_OUT = 4
ADVECT = 2
ADVECT_VERT_COORD = 1
TURB_DX_TROP = 0
TURB_DX_STRAT = 0
TURB_DZ_TROP = 0
TURB_DZ_STRAT = 0.0
TURB_MESOX = 0.0
TURB_MESOZ = 0.0
DIRECTION = 1
METBASE = ${metbase}
MET_TROPO = 3
TDEC_TROP = 259200
TDEC_STRAT = 259200
DT_MOD = ${dt_mod}
DT_MET = ${dt_met}
T_START = ${t0}
T_STOP = ${t1}
ATM_DT_OUT = ${met_dt_out}
GRID_DT_OUT = 25920000
GRID_LON0 = ${grid_lon_start}
GRID_LON1 = $((grid_lon_start + grid_lon_size))
GRID_LAT0 = ${grid_lat_start}
GRID_LAT1 = $((grid_lat_start + grid_lat_size))
GRID_NX = ${grid_nx}
GRID_NY = ${grid_ny}
ENSEMBLE = 1
NQ = 5
QNT_NAME[0] = idx
QNT_NAME[1] = zeta
QNT_NAME[2] = m
QNT_NAME[3] = subdomain
QNT_NAME[4] = destination
DD_HALOS_SIZE = 2
DD_SUBDOMAINS_MERIDIONAL = ${domains_lat}
DD_SUBDOMAINS_ZONAL = ${domains_lon}
NAP_MAX = ${npmax}
EOF

echo "[INFO] Create initialization files for each subdomain"
for i in $(seq 0 $((processes-1))); do
    mkdir -p data/data.$i

    ilat=$(($i % domains_lat))
    ilon=$(($i / domains_lon))

    # Calculate floating point coordinates
    init_lon0=$(echo "scale=6; $grid_lon_start + $grid_lon_d * ($ilon * $subgrid_lon_size + $subgrid_pos_lon_offset)" | bc)
    init_lon1=$(echo "scale=6; $grid_lon_start + $grid_lon_d * ($ilon * $subgrid_lon_size + $subgrid_lon_size - $subgrid_pos_lon_offset)" | bc)
    init_lat0=$(echo "scale=6; $grid_lat_start + $grid_lat_d * ($ilat * $subgrid_lat_size + $subgrid_pos_lat_offset)" | bc)
    init_lat1=$(echo "scale=6; $grid_lat_start + $grid_lat_d * ($ilat * $subgrid_lat_size + $subgrid_lat_size - $subgrid_pos_lat_offset)" | bc)
    init_dlon=$(echo "scale=6; $subgrid_pos_dlon" | bc)
    init_dlat=$(echo "scale=6; $subgrid_pos_dlat" | bc)

    # take the smaller value for lat0 or lon0
    if (( $(echo "$init_lon0 > $init_lon1" | bc -l) )); then
        tmp=$init_lon0
        init_lon0=$init_lon1
        init_lon1=$tmp
    fi
    if (( $(echo "$init_lat0 > $init_lat1" | bc -l) )); then
        tmp=$init_lat0
        init_lat0=$init_lat1
        init_lat1=$tmp
    fi

    # print initialization box
    echo "[INFO] Domain $i (ilon=$ilon, ilat=$ilat):"
    echo "[INFO]   Lon0: $init_lon0"
    echo "[INFO]   Lon1: $init_lon1"
    echo "[INFO]   Lat0: $init_lat0"
    echo "[INFO]   Lat1: $init_lat1"
    echo "[INFO]   Dlon: $init_dlon"
    echo "[INFO]   Dlat: $init_dlat"

    $mptrac_dir/src/atm_init $work_dir/data/config.ctl \
                $work_dir/data/data.$i/init.nc \
                INIT_T0 "$t0" INIT_T1 "$t0" \
                INIT_Z0 10.0 INIT_Z1 10.0 \
                INIT_LON0 $init_lon0 \
                INIT_LON1 $init_lon1 \
                INIT_DLON $init_dlon \
                INIT_LAT0 $init_lat0 \
                INIT_LAT1 $init_lat1 \
                INIT_DLAT $init_dlat \
                INIT_IDX_OFFSET $((particles_num_subdomain * i)) > $work_dir/data/atm_generation.log 2>&1
    
    cat > data/data.$i/config.ctl <<EOF
MET_PRESS_LEVEL_DEF = ${met_press_level_def}
MET_CLAMS = ${met_clams}
MET_VERT_COORD = ${met_vert_coord}
MET_DT_OUT = ${met_dt_out}
ATM_TYPE = ${atm_type}
ATM_TYPE_OUT = 0
ADVECT = 2
ADVECT_VERT_COORD = 0
TURB_DX_TROP = 0
TURB_DX_STRAT = 0
TURB_DZ_TROP = 0
TURB_DZ_STRAT = 0.0
TURB_MESOX = 0.0
TURB_MESOZ = 0.0
DIRECTION = 1
METBASE = ${metbase}
MET_TROPO = 3
TDEC_TROP = 259200
TDEC_STRAT = 259200
DT_MOD = ${dt_mod}
DT_MET = ${dt_met}
T_START = ${t0}
T_STOP = ${t1}
ATM_DT_OUT = ${met_dt_out}
GRID_DT_OUT = 25920000
GRID_LON0 = ${grid_lon_start}
GRID_LON1 = $((grid_lon_start + grid_lon_size))
GRID_LAT0 = ${grid_lat_start}
GRID_LAT1 = $((grid_lat_start + grid_lat_size))
GRID_NX = ${grid_nx}
GRID_NY = ${grid_ny}
ENSEMBLE = 1
NQ = 5
QNT_NAME[0] = idx
QNT_NAME[1] = zeta
QNT_NAME[2] = m
QNT_NAME[3] = subdomain
QNT_NAME[4] = destination
DD_HALOS_SIZE = 2
DD_SUBDOMAINS_MERIDIONAL = ${domains_lat}
DD_SUBDOMAINS_ZONAL = ${domains_lon}
NAP_MAX = ${npmax}
EOF
done

#################################################################################################################
# Generate wind data files
#################################################################################################################

cd $work_dir/data/wind_data
echo "[INFO] Generate wind data files"

cat > wind.ctl <<EOF
WIND_T0 = ${t0}
WIND_NX = ${met_grid_nx}
WIND_NY = ${met_grid_ny}
WIND_NZ = ${met_grid_nz}
WIND_Z0 = 0
WIND_Z1 = 60
WIND_U0 = ${wind_speed}
WIND_U1 = ${wind_speed}
WIND_W0 = 0
WIND_ALPHA = ${wind_alpha}
EOF

for i in {0..6}; do
    new_t0=$(echo "scale=6; $t0 + $i * $dt_met" | bc)
    
    # Update the T0 value in wind.ctl
    sed -i "s/^WIND_T0 = .*/WIND_T0 = $new_t0/" wind.ctl
    
    # Run the wind command
    $mptrac_dir/src/wind wind.ctl wind > $work_dir/data/wind_generation.log 2>&1
    
    echo "[INFO] Generated file $((i+1)) with T0 = $new_t0"
done

#################################################################################################################
# Running MPTRAC with domain decomposition
#################################################################################################################

echo "[INFO] Running MPTRAC with domain decomposition with $processes processes"
cd $work_dir/data
if [[ "$hpc_mode" == "hpc" ]]; then
    srun --nodes 1 --ntasks $processes --ntasks-per-node $processes --cpus-per-task 1 ../../../src/trac dirlist.tab config.ctl init.nc ATM_BASENAME atm > $work_dir/data/srun_output.log 2>&1
else
    mpirun -np $processes ../../../src/trac dirlist.tab config.ctl init.nc ATM_BASENAME atm > $work_dir/data/mpirun_output.log 2>&1
fi
cd $work_dir

#################################################################################################################
# Compare results
#################################################################################################################

echo -e "\nComparing .tab files recursively..."
error=0
if [ -d "data.ref" ]; then
    # Find all .tab files in data.ref and compare with corresponding files in data
    while IFS= read -r -d '' ref_file; do
        # Get relative path from data.ref
        rel_path="${ref_file#data.ref/}"
        data_file="data/$rel_path"
        
        if [ -f "$data_file" ]; then
            diff -q -s "$data_file" "$ref_file" || error=1
        else
            echo "Missing file: $data_file"
            error=1
        fi
    done < <(find data.ref -type f -name "*.tab" -print0)

    # Also check for extra .tab files in data that don't exist in data.ref
    while IFS= read -r -d '' data_file; do
        # Get relative path from data
        rel_path="${data_file#data/}"
        ref_file="data.ref/$rel_path"
        
        if [ ! -f "$ref_file" ]; then
            echo "Extra .tab file in data: $data_file"
            error=1
        fi
    done < <(find data -type f -name "*.tab" -print0)
else
    echo "Reference directory data.ref not found - skipping comparison"
fi

if [ $error -eq 0 ]; then
    echo "[INFO] Test completed successfully."
else
    echo "[ERROR] Some files do not match."
fi

exit $error
