#! /bin/bash
set -e

# Usage: ./run.sh [compile] [hpc] [skip_wind] [skip_run] [skip_compare]
# Arguments:
#   compile: Enable compilation step (default is to skip compilation)
#   hpc: Run in HPC mode (uses ml and srun), otherwise uses mpirun
#   skip_wind: Skip wind data generation
#   skip_run: Skip MPTRAC simulation run
#   skip_compare: Skip result comparison
# Arguments can be provided in any order

#################################################################################################################
# Setup directory structure, parse arguments
#################################################################################################################

work_dir=$(pwd)
mptrac_dir=$(cd ../../ && pwd)
meteo_data_dir=$work_dir/data/wind_data
metbase=$meteo_data_dir/wind    # path to meteo files
echo "[INFO] MPTRAC directory: $mptrac_dir"
echo "[INFO] Working directory: $work_dir"
echo "[INFO] Meteo data directory: $meteo_data_dir"
echo "[INFO] Meteo base path: $metbase"
cd $work_dir
rm -rf data
mkdir -p data
mkdir -p $meteo_data_dir

# Parse arguments - check all arguments for keywords
compile=""
hpc_mode=""
skip_wind=""
skip_run=""
skip_compare=""
for arg in "$@"; do
    case "$arg" in
        "compile")
            compile="compile"
            ;;
        "hpc")
            hpc_mode="hpc"
            ;;
        "skip_wind")
            skip_wind="skip_wind"
            ;;
        "skip_run")
            skip_run="skip_run"
            ;;
        "skip_compare")
            skip_compare="skip_compare"
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

#################################################################################################################
# Define simulation parameters
#################################################################################################################

ntasks_per_node=9       # fixed number of tasks per node

use_gpu=0                           # enable GPU support in MPTRAC compilation (0=off, 1=on)
cpus_per_task=1                     # number of CPUs per task

domains_lon=3                       # number of subdomains in longitudinal direction
domains_lat=3                       # number of subdomains in latitudinal direction

# Set parameters for atmospheric initialization
use_evenly_distributed=0            # use INIT_EVENLY=1 for cosine-weighted particle distribution (0=off, 1=on)
particles_subdomain_lon_stride=4    # number of particles in longitudinal direction per subdomain
particles_subdomain_lat_stride=4    # number of particles in latitudinal direction per subdomain
grid_lon_start=0                    # starting longitude of the grid for the atmospheric initialization
grid_lon_d=1                        # longitude direction (1 or -1)
grid_lon_size=360                   # size of the grid for the atmospheric initialization in longitudinal direction
grid_lat_start=90                   # starting latitude of the grid for the atmospheric initialization
grid_lat_d=-1                       # latitude direction (1 or -1)
grid_lat_size=180                   # size of the grid for the atmospheric initialization in latitudinal direction

# Set meteo file parameters
met_press_level_def=-1  # pressure level definition
met_clams=0             # follow CLaMS filename guideline
met_vert_coord=0        # original vertical coordinate type
dt_met=3600             # read time interval for meteo files
dt_mod=600              # model time step
met_dt_out=3600         # output time interval for meteo files

# Set simulation time parameters
simulation_hours=6      # duration of simulation in hours
ntimesteps=$((simulation_hours * 3600 / dt_met + 1))  # number of time steps (including initial time)
atm_type=3              # original atmospheric data type
met_grid_nx=360         # number of grid points in longitudinal direction of meteo files
met_grid_ny=181         # number of grid points in latitudinal direction of meteo files
met_grid_nz=60          # number of vertical levels of meteo files
dex_glob=364            # allocated size of meteo grid in longitudinal direction
dey_glob=186            # allocated size of meteo grid in latitudinal direction
dez_glob=90             # allocated size of meteo grid in vertical direction
wind_alpha=90           # angle for wind rotation  in degrees
wind_speed=50           # wind speed in m/s

# Calculate derived parameters
ntasks=$((domains_lon * domains_lat))
nnodes=$(( (ntasks + ntasks_per_node - 1) / ntasks_per_node ))
particles_num_subdomain=$((particles_subdomain_lon_stride * particles_subdomain_lat_stride))
npmax=$((particles_num_subdomain * 2))
nparticles=$((particles_num_subdomain * ntasks))
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
echo "[INFO] ntasks_per_node: $ntasks_per_node"
echo "[INFO] cpus_per_task: $cpus_per_task"
echo "[INFO] ntasks: $ntasks"
echo "[INFO] nnodes: $nnodes"
echo "[INFO] use_evenly_distributed: $use_evenly_distributed"
echo "[INFO] use_gpu: $use_gpu"
echo "[INFO] simulation_hours: $simulation_hours"
echo "[INFO] ntimesteps: $ntimesteps"
echo "[INFO] particles_num_subdomain: $particles_num_subdomain"
echo "[INFO] npmax: $npmax"
echo "[INFO] nparticles: $nparticles"
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

# Set environment
export LD_LIBRARY_PATH=$mptrac_dir/libs/build/lib/:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=$cpus_per_task
export OMP_PROC_BIND=true
export OMP_PLACES=sockets
#export LANG=C
#export LC_ALL=C

defs="-DNP=${npmax} -DEX_GLOB=${dex_glob} -DEY_GLOB=${dey_glob} -DEP=${dez_glob} -DEX=${dex} -DEY=${dey} -DNQ=5"

echo "[INFO] defs: $defs"
echo "[INFO] LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

# Set compiler based on mode
if [[ "$hpc_mode" == "hpc" ]]; then
    compiler="nvc"
    echo "[INFO] Using nvc compiler (HPC mode)"
else
    compiler="gcc"
    echo "[INFO] Using gcc compiler (non-HPC mode)"
fi

# Check if compilation should be skipped
if [[ "$compile" != "compile" ]]; then
    echo "[INFO] Skipping compilation (default behavior - use 'compile' argument to enable)"
else
    echo "[INFO] Compile MPTRAC with domain decomposition using $compiler"
    cd $mptrac_dir/src
    make clean
    
    # First compile all executables with GPU support (including trac)
    echo "[INFO] Compiling all executables"
    make -j DD=1 MPI=1 STATIC=0 GPU=$use_gpu THRUST=1 COMPILER=$compiler DEFINES="$defs" || exit 1
    
    # Then recompile atm_init without GPU support (atm_init doesn't work with GPU=1)
    echo "[INFO] Recompiling atm_init without GPU support"
    rm -f $mptrac_dir/src/atm_init
    rm -f $mptrac_dir/src/wind
    make atm_init wind DD=1 MPI=1 STATIC=0 GPU=0 THRUST=1 COMPILER=$compiler DEFINES="$defs" || exit 1
fi

#################################################################################################################
# Create atmospheric initialization files & control files for each subdomain
#################################################################################################################

cd $work_dir

# Set timerange of simulation
t0=$($mptrac_dir/src/time2jsec 2022 6 2 0 0 0 0)
t1=$($mptrac_dir/src/time2jsec 2022 6 2 $simulation_hours 0 0 0)

echo "[INFO] write dirlist.tab and config.ctl to data directory"
> data/dirlist.tab  # Clear the file first
for i in $(seq 0 $((ntasks-1))); do
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
if [[ $use_evenly_distributed -eq 1 ]]; then
    echo "[WARNING] Using INIT_EVENLY=1: Particle count per subdomain may vary due to cosine weighting"
    echo "[WARNING] This could affect load balancing between subdomains"
fi
for i in $(seq 0 $((ntasks-1))); do
    mkdir -p data/data.$i

    ilat=$((i % domains_lat))
    ilon=$((i / domains_lon))

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
                INIT_EVENLY $use_evenly_distributed \
                INIT_IDX_OFFSET $((particles_num_subdomain * i)) >> $work_dir/data/atm_generation.log 2>&1
    
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

# Check if wind data generation should be skipped
if [[ "$skip_wind" == "skip_wind" ]]; then
    echo "[INFO] Skipping wind data generation (argument: $skip_wind)"
else
    cd $meteo_data_dir
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
WIND_LAT_REVERSE = 1
EOF

    for i in $(seq 0 $((ntimesteps-1))); do
        new_t0=$(echo "scale=6; $t0 + $i * $dt_met" | bc)
        
        # Update the T0 value in wind.ctl
        sed -i "s/^WIND_T0 = .*/WIND_T0 = $new_t0/" wind.ctl
        
        # Run the wind command
        $mptrac_dir/src/wind wind.ctl wind >> $work_dir/data/wind_generation.log 2>&1
        
        echo "[INFO] Generated file $((i+1))/$ntimesteps with T0 = $new_t0"
    done
fi

#################################################################################################################
# Running MPTRAC with domain decomposition
#################################################################################################################

# Check if MPTRAC run should be skipped
if [[ "$skip_run" == "skip_run" ]]; then
    echo "[INFO] Skipping MPTRAC simulation run (argument: $skip_run)"
else
    echo "[INFO] Running MPTRAC with domain decomposition with $ntasks tasks"
    cd $work_dir/data
    if [[ "$hpc_mode" == "hpc" ]]; then
        if [[ $use_gpu -eq 1 ]]; then
            echo "[INFO] Running with GPU support on HPC (partition=booster, gres=gpu:$ntasks_per_node, cpus-per-task=$cpus_per_task)"
            srun --partition=booster --gres=gpu:$ntasks_per_node --nodes $nnodes --ntasks $ntasks --account=gsp25 --time=00:59:00 \
                --ntasks-per-node $ntasks_per_node --cpus-per-task $cpus_per_task \
                --output=$work_dir/data/mptrac_gpu_${domains_lon}x${domains_lat}_n${nnodes}_t${ntasks}.out \
                --error=$work_dir/data/mptrac_gpu_${domains_lon}x${domains_lat}_n${nnodes}_t${ntasks}.err \
                $mptrac_dir/src/trac dirlist.tab config.ctl init.nc ATM_BASENAME atm
        else
            echo "[INFO] Running with CPU-only on HPC (cpus-per-task=$cpus_per_task)"
            srun --nodes $nnodes --ntasks $ntasks --account=gsp25 --time=00:59:00 \
                --ntasks-per-node $ntasks_per_node --cpus-per-task $cpus_per_task \
                --output=$work_dir/data/mptrac_cpu_${domains_lon}x${domains_lat}_n${nnodes}_t${ntasks}.out \
                --error=$work_dir/data/mptrac_cpu_${domains_lon}x${domains_lat}_n${nnodes}_t${ntasks}.err \
                $mptrac_dir/src/trac dirlist.tab config.ctl init.nc ATM_BASENAME atm
        fi
    else
        mpirun --oversubscribe -np $ntasks $mptrac_dir/src/trac dirlist.tab config.ctl init.nc ATM_BASENAME atm > $work_dir/data/mpirun_output.log 2>&1
    fi
fi

#################################################################################################################
# Compare results
#################################################################################################################

# Check if result comparison should be skipped
if [[ "$skip_compare" == "skip_compare" ]]; then
    echo "[INFO] Skipping result comparison (argument: $skip_compare)"
    echo "[INFO] Test completed (comparison skipped)."
    exit 0
fi

cd $work_dir
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
