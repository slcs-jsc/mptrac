#! /bin/bash
set -e

# Usage: ./run.sh [compile] [hpc] [skip_wind] [skip_run] [skip_compare]
# Arguments:
#   compile: Enable compilation step (default is to skip compilation)
#   hpc: Optional. Run in HPC mode (uses module and srun); default is local mode via mpirun
#   skip_wind: Skip wind data generation
#   skip_run: Skip MPTRAC simulation run
#   skip_compare: Skip result comparison
# Set DD_EXCHANGE_MODE=0 or 1 to select the exchange path (default: 0).
# Set DD_USE_GPU=1 for an OpenACC build on a GPU node (default: 0).
# DD_DOMAINS_LON/LAT, DD_PARTICLES_LON_STRIDE/LAT_STRIDE, and DD_NP_MAX
# configure non-default layouts and particle counts for benchmark runs.
# DD_SORT_DT, DD_SIMULATION_HOURS, DD_THRUST, DD_COMPILER, and DD_WIND_DIR
# select the run duration, sorting backend, toolchain, and input wind fields.
# DD_COMPARE_DIR selects reference .tab files (default: data.ref).
# DD_SAVE_DIR saves this run's .tab files for a later CPU/GPU comparison.
# DD_ISOSURF=3 adds the iso_var quantity for the isosurface migration test.
# DD_INIT_Z0/Z1/DZ and DD_INIT_IDX_STRIDE set its two heights and unique IDs;
# DD_WIND_ALPHA=0 makes particles cross the zonal subdomain boundary.
# From tests/dd_test, for example:
#   DD_EXCHANGE_MODE=1 DD_SAVE_DIR=data.cpu ./run.sh compile skip_wind
#   DD_EXCHANGE_MODE=1 DD_USE_GPU=1 DD_COMPARE_DIR=data.cpu ./run.sh compile skip_wind
# Add "hpc" on a GPU node using the configured Slurm launcher.
# Arguments can be provided in any order

#################################################################################################################
# Setup directory structure, parse arguments
#################################################################################################################

work_dir=$(pwd)
script_dir=$(cd "$(dirname "$0")" && pwd)
mptrac_dir=$(cd "$script_dir/../.." && pwd)
meteo_data_dir=${DD_WIND_DIR:-$work_dir/data/wind_data}
metbase=$meteo_data_dir/wind    # path to meteo files
echo "[INFO] MPTRAC directory: $mptrac_dir"
echo "[INFO] Working directory: $work_dir"
echo "[INFO] Meteo data directory: $meteo_data_dir"
echo "[INFO] Meteo base path: $metbase"
cd "$work_dir"

# Parse arguments - check all arguments for keywords
compile=""
hpc_mode="local"
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

# Keep generated wind fields when skip_wind is requested. This lets a CPU
# baseline and a GPU run use the same meteorological input without rebuilding it.
if [[ "$skip_wind" != "skip_wind" ]]; then
    rm -rf data
fi
mkdir -p data "$meteo_data_dir"

# Check if running on HPC system
if [[ "$hpc_mode" == "hpc" ]]; then
    echo "[INFO] Running in HPC mode"
    module purge
    module load Stages/2026
    module load GCC/14.3.0 ParaStationMPI/5.13.0-1
    module load GSL
    module load netCDF/4.9.3
else
    echo "[INFO] Running in non-HPC mode"
fi

#################################################################################################################
# Define simulation parameters
#################################################################################################################

ntasks_per_node=${DD_TASKS_PER_NODE:-9}

use_gpu=${DD_USE_GPU:-0}            # enable GPU support in MPTRAC compilation (0=off, 1=on)
cpus_per_task=1                     # number of CPUs per task
mpi_cc=${DD_MPICC:-mpicc}
mpi_run=${DD_MPIRUN:-mpirun}

domains_lon=${DD_DOMAINS_LON:-3}    # number of zonal subdomains
domains_lat=${DD_DOMAINS_LAT:-3}    # number of meridional subdomains
dd_exchange_mode=${DD_EXCHANGE_MODE:-0}
use_thrust=${DD_THRUST:-1}
isosurf=${DD_ISOSURF:-0}
if [[ $isosurf -eq 3 ]]; then
    nq=6
    iso_quantity='QNT_NAME[5] = iso_var'
else
    nq=5
    iso_quantity=''
fi

# Set parameters for atmospheric initialization
use_evenly_distributed=0            # use INIT_EVENLY=1 for cosine-weighted particle distribution (0=off, 1=on)
particles_subdomain_lon_stride=${DD_PARTICLES_LON_STRIDE:-4}
particles_subdomain_lat_stride=${DD_PARTICLES_LAT_STRIDE:-4}
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
simulation_hours=${DD_SIMULATION_HOURS:-6}
ntimesteps=$((simulation_hours * 3600 / dt_met + 1))  # number of time steps (including initial time)
atm_type=3              # original atmospheric data type
met_grid_nx=360         # number of grid points in longitudinal direction of meteo files
met_grid_ny=181         # number of grid points in latitudinal direction of meteo files
met_grid_nz=60          # number of vertical levels of meteo files
dex_glob=364            # allocated size of meteo grid in longitudinal direction
dey_glob=186            # allocated size of meteo grid in latitudinal direction
dez_glob=90             # allocated size of meteo grid in vertical direction
wind_alpha=${DD_WIND_ALPHA:-90}  # angle for wind rotation in degrees
wind_speed=50           # wind speed in m/s

# Calculate derived parameters
ntasks=$((domains_lon * domains_lat))
if (( ntasks_per_node > ntasks )); then
    ntasks_per_node=$ntasks
fi
nnodes=$(( (ntasks + ntasks_per_node - 1) / ntasks_per_node ))
particles_num_subdomain=$((particles_subdomain_lon_stride * particles_subdomain_lat_stride))
particle_id_stride=${DD_INIT_IDX_STRIDE:-$particles_num_subdomain}
npmax=${DD_NP_MAX:-$((particles_num_subdomain * 2))}
nparticles=$((particles_num_subdomain * ntasks))
subgrid_lon_size=$(echo "scale=6; $grid_lon_size / $domains_lon" | bc)
subgrid_lat_size=$(echo "scale=6; $grid_lat_size / $domains_lat" | bc)
subgrid_pos_dlon=$(echo "scale=6; $subgrid_lon_size / $particles_subdomain_lon_stride" | bc)
subgrid_pos_dlat=$(echo "scale=6; $subgrid_lat_size / $particles_subdomain_lat_stride" | bc)
subgrid_pos_lon_offset=$(echo "scale=6; $subgrid_pos_dlon / 2" | bc)
subgrid_pos_lat_offset=$(echo "scale=6; $subgrid_pos_dlat / 2" | bc)
dex=$((met_grid_nx/domains_lon + 10))
dey=$((met_grid_ny/domains_lat + 10))
grid_nx=$dex
grid_ny=$dey

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

defs="-DNP=${npmax} -DEX_GLOB=${dex_glob} -DEY_GLOB=${dey_glob} -DEP=${dez_glob} -DEX=${dex} -DEY=${dey} -DNQ=${nq}"

echo "[INFO] defs: $defs"
echo "[INFO] LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

# Set compiler based on mode
if [[ "$hpc_mode" == "hpc" ]]; then
    compiler=${DD_COMPILER:-nvc}
    echo "[INFO] Using $compiler compiler (HPC mode)"
else
    compiler=${DD_COMPILER:-gcc}
    echo "[INFO] Using $compiler compiler (non-HPC mode)"
fi

# Check if compilation should be skipped
if [[ "$compile" != "compile" ]]; then
    echo "[INFO] Skipping compilation (default behavior - use 'compile' argument to enable)"
else
    echo "[INFO] Compile MPTRAC with domain decomposition using $compiler"
    cd "$mptrac_dir"/src
    make clean
    
    # Build input generators on the CPU, then force a separate trac object
    # build with the requested accelerator setting.
    make -j atm_init wind time2jsec DD=1 MPI=1 STATIC=0 GPU=0 THRUST=0 COMPILER=$compiler MPICC="$mpi_cc" DEFINES="$defs"
    make -B trac DD=1 MPI=1 STATIC=0 GPU=$use_gpu THRUST=$use_thrust COMPILER=$compiler MPICC="$mpi_cc" DEFINES="$defs"
fi

#################################################################################################################
# Create atmospheric initialization files & control files for each subdomain
#################################################################################################################

cd "$work_dir"

# Set timerange of simulation
t0=$("$mptrac_dir"/src/time2jsec 2022 6 2 0 0 0 0)
t1=$("$mptrac_dir"/src/time2jsec 2022 6 2 "$simulation_hours" 0 0 0)

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
MET_CONVENTION = 0
ATM_TYPE = ${atm_type}
ATM_TYPE_OUT = 4
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
GRID_LAT0 = $((grid_lat_start - grid_lat_size))
GRID_LAT1 = ${grid_lat_start}
GRID_NX = ${grid_nx}
GRID_NY = ${grid_ny}
ENSEMBLE = 1
NQ = ${nq}
QNT_NAME[0] = idx
QNT_NAME[1] = zeta
QNT_NAME[2] = m
QNT_NAME[3] = current_subdomain
QNT_NAME[4] = target_subdomain
${iso_quantity}
ISOSURF = ${isosurf}
DD_HALOS_SIZE = 2
DD_SUBDOMAINS_MERIDIONAL = ${domains_lat}
DD_SUBDOMAINS_ZONAL = ${domains_lon}
DD_EXCHANGE_MODE = ${dd_exchange_mode}
DD_SORT_DT = ${DD_SORT_DT:-1800}
NAP_MAX = ${npmax}
EOF

echo "[INFO] Create initialization files for each subdomain"
if [[ $use_evenly_distributed -eq 1 ]]; then
    echo "[WARNING] Using INIT_EVENLY=1: Particle count per subdomain may vary due to cosine weighting"
    echo "[WARNING] This could affect load balancing between subdomains"
fi
for i in $(seq 0 $((ntasks-1))); do
    mkdir -p data/"data.$i"

    ilat=$((i % domains_lat))
    ilon=$((i / domains_lat))

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

    # Leave a fraction of one grid step beyond the last point so rounding
    # in the decimal step size cannot silently omit a whole particle row.
    init_lon1=$(echo "scale=6; $init_lon1 + $subgrid_pos_dlon / 4" | bc)
    init_lat1=$(echo "scale=6; $init_lat1 + $subgrid_pos_dlat / 4" | bc)

    # print initialization box
    echo "[INFO] Domain $i (ilon=$ilon, ilat=$ilat):"
    echo "[INFO]   Lon0: $init_lon0"
    echo "[INFO]   Lon1: $init_lon1"
    echo "[INFO]   Lat0: $init_lat0"
    echo "[INFO]   Lat1: $init_lat1"
    echo "[INFO]   Dlon: $init_dlon"
    echo "[INFO]   Dlat: $init_dlat"

    "$mptrac_dir"/src/atm_init "$work_dir"/data/config.ctl \
                 "$work_dir"/data/"data.$i"/init.nc \
                 INIT_T0 "$t0" INIT_T1 "$t0" \
                 INIT_Z0 "${DD_INIT_Z0:-10.0}" INIT_Z1 "${DD_INIT_Z1:-10.0}" \
                 INIT_DZ "${DD_INIT_DZ:-1.0}" \
                 INIT_LON0 "$init_lon0" \
                 INIT_LON1 "$init_lon1" \
                 INIT_DLON "$init_dlon" \
                 INIT_LAT0 "$init_lat0" \
                 INIT_LAT1 "$init_lat1" \
                 INIT_DLAT "$init_dlat" \
                 INIT_EVENLY "$use_evenly_distributed" \
                 INIT_IDX_OFFSET $((particle_id_stride * i)) >> "$work_dir"/data/atm_generation.log 2>&1
    
    cat > data/"data.$i"/config.ctl <<EOF
MET_PRESS_LEVEL_DEF = ${met_press_level_def}
MET_CLAMS = ${met_clams}
MET_VERT_COORD = ${met_vert_coord}
MET_DT_OUT = ${met_dt_out}
MET_CONVENTION = 0
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
GRID_LAT0 = $((grid_lat_start - grid_lat_size))
GRID_LAT1 = ${grid_lat_start}
GRID_NX = ${grid_nx}
GRID_NY = ${grid_ny}
ENSEMBLE = 1
NQ = ${nq}
QNT_NAME[0] = idx
QNT_NAME[1] = zeta
QNT_NAME[2] = m
QNT_NAME[3] = current_subdomain
QNT_NAME[4] = target_subdomain
${iso_quantity}
ISOSURF = ${isosurf}
DD_HALOS_SIZE = 2
DD_SUBDOMAINS_MERIDIONAL = ${domains_lat}
DD_SUBDOMAINS_ZONAL = ${domains_lon}
DD_EXCHANGE_MODE = ${dd_exchange_mode}
DD_SORT_DT = ${DD_SORT_DT:-1800}
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
    cd "$meteo_data_dir"
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
        "${WIND_EXEC:-$mptrac_dir/src/wind}" wind.ctl wind >> "$work_dir"/data/wind_generation.log 2>&1
        
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
    cd "$work_dir"/data
    if [[ "$hpc_mode" == "hpc" ]]; then
        if [[ $use_gpu -eq 1 ]]; then
            echo "[INFO] Running with GPU support on HPC (partition=booster, gres=gpu:$ntasks_per_node, cpus-per-task=$cpus_per_task)"
            srun --partition=booster --gres=gpu:$ntasks_per_node --nodes $nnodes --ntasks $ntasks --account=gsp25 --time=00:59:00 \
                --ntasks-per-node $ntasks_per_node --cpus-per-task $cpus_per_task \
                --output="$work_dir"/data/mptrac_gpu_${domains_lon}x${domains_lat}_n${nnodes}_t${ntasks}.out \
                --error="$work_dir"/data/mptrac_gpu_${domains_lon}x${domains_lat}_n${nnodes}_t${ntasks}.err \
                "$mptrac_dir"/src/trac dirlist.tab config.ctl init.nc ATM_BASENAME atm
        else
            echo "[INFO] Running with CPU-only on HPC (cpus-per-task=$cpus_per_task)"
            srun --nodes $nnodes --ntasks $ntasks --account=clams-esm --time=00:59:00 \
                --ntasks-per-node $ntasks_per_node --cpus-per-task $cpus_per_task \
                --output="$work_dir"/data/mptrac_cpu_${domains_lon}x${domains_lat}_n${nnodes}_t${ntasks}.out \
                --error="$work_dir"/data/mptrac_cpu_${domains_lon}x${domains_lat}_n${nnodes}_t${ntasks}.err \
                "$mptrac_dir"/src/trac dirlist.tab config.ctl init.nc ATM_BASENAME atm
        fi
    else
        "$mpi_run" --oversubscribe -np $ntasks "$mptrac_dir"/src/trac dirlist.tab config.ctl init.nc ATM_BASENAME atm > "$work_dir"/data/mpirun_output.log 2>&1
    fi
fi

#################################################################################################################
# Compare results
#################################################################################################################

cd "$work_dir"
if [ -n "${DD_SAVE_DIR:-}" ]; then
    mkdir -p "$DD_SAVE_DIR"
    while IFS= read -r -d '' data_file; do
        rel_path=${data_file#data/}
        mkdir -p "$DD_SAVE_DIR/$(dirname "$rel_path")"
        cp "$data_file" "$DD_SAVE_DIR/$rel_path"
    done < <(find data -type f -name '*.tab' -print0)
    echo "[INFO] Saved .tab files to $DD_SAVE_DIR"
fi

# Check if result comparison should be skipped
if [[ "$skip_compare" == "skip_compare" ]]; then
    echo "[INFO] Skipping result comparison (argument: $skip_compare)"
    echo "[INFO] Test completed (comparison skipped)."
    exit 0
fi

compare_dir=${DD_COMPARE_DIR:-data.ref}
if [ -d "$compare_dir" ]; then
    bash "$script_dir/check.sh" data "$compare_dir"
else
    echo "[ERROR] Reference directory $compare_dir not found"
    exit 1
fi

echo "[INFO] Test completed successfully."
