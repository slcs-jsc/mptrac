#! /bin/bash

# Check arguments...
if [ $# -ne 18 ] ; then
    cat <<EOF
MPTRAC benchmarking script

usage: $0 <system> <account> <runtime> <mptrac> <libs> <compiler> <nodes> <ntasks-per-node> <cpus-per-task> <gpu> <npmin> <npmax> <npfac> <meteo> <rng> <sort> <phys> <cache>

parameter choices:

  system: specify target system for benchmarking
    jwb = JUWELS Booster
    jwc_cpu = JUWELS Cluster CPU-only
    jwc_gpu = JUWELS Cluster V100
    jrc_cpu = JURECA-DC CPU-only
    jrc_gpu = JURECA-DC A100
    jrc_gh200 = JURECA Grace Hopper 200
    jrc_h100 = JURECA Hopper 100

  account: specify compute time budget
    jsc = JSC compute budget
    slmet = SDLCS compute budget
    exaww = WarmWorld compute budget

  runtime: specify maximum runtime (e.g. 00:30:00)

  mptrac: specify whether MPTRAC should be recompiled
    avail = use available copy of MPTRAC as is
    clone = remove existing copy of MPTRAC and make new clone from git repository
    pull = update existing copy of MPTRAC via pull request to git repository

  libs: specify whether libraries should be recompiled
    avail = assume libraries for MPTRAC are already available and compiled
    compute = recompile required libraries on compute node
    login = recompile required libraries on login node

  compiler: select C compiler
    gcc = GNU C compiler
    nvc = NVIDIA C compiler

  nodes | ntasks-per-node | cpus-per-task:
    provide slurm job configuration for the test job

  gpu: specify whether GPU offloading should be used
    0 = none
    1 = OpenACC
    2 = OpenACC + pinning

  npmin | npmax | npfac:
    provide minimum, maximum, and scaling factor for number of air parcels

  meteo: select meteorological data
    era5 = ERA5 reanalysis
    erai = ERA-Interim reanalysis

  rng: select random number generator
    0 = GSL
    1 = Squares
    2 = cuRAND

  sort: sorting time interval (0=disable sorting)

  phys: select simulation type
    full = full-physics simulation
    traj = trajectory-only test case

  cache: select disk cache
    0 = none
    1 = disk cache
    2 = async-I/O
EOF
    exit
fi

# Set variables...
trac=$(pwd)/mptrac
system=$1
account=$2
runtime=$3
mptrac=$4
libs=$5
compiler=$6
nodes=$7
ntasks_per_node=$8
cpus_per_task=$9
gpu=${10}
npmin=${11}
npmax=${12}
npfac=${13}
meteo=${14}
rng=${15}
sort=${16}
phys=${17}
cache=${18}

# Update MPTRAC repository...
echo -e "\nUpdating MPTRAC repository..."
if [ $mptrac = "clone" ] ; then
    rm -rf mptrac ; git clone https://github.com/slcs-jsc/mptrac.git || exit    
elif [ $mptrac = "pull" ] ; then
    cd mptrac && git pull || exit ; cd -
else
    echo "Use as is!"
fi

# Compile libraries...
echo -e "\nCompile libraries..."
if [ $libs = "login" ] ; then
    cd $trac/libs && ./build.sh -a ; cd -
else
    echo "Compile on compute node or assume already available!"
fi

# Get meteo data...
echo -e "\nDownloading meteo data..."
[ -s meteo ] \
    || wget --mirror --no-parent --no-host-directories --execute robots=off --reject="index.html*" --cut-dirs=5 https://datapub.fz-juelich.de/slcs/mptrac/data/projects/benchmarking/meteo/ \
	&& echo "Already available!"
du -h meteo || exit

# Execute batch job...
echo -e "\nExecute batch job..."
ntasks=$(echo $nodes $ntasks_per_node | awk '{print $1*$2}')
slurmset="--account=$account --time=$runtime --nodes=$nodes --ntasks=$ntasks --ntasks-per-node=$ntasks_per_node --cpus-per-task=$cpus_per_task --disable-dcgm --hint=nomultithread"
[ $gpu != "0" ] && slurmset+=" --gres=gpu:1"
if [ $system = "jwb" ] ; then
    slurmset+=" --partition=booster"
elif [ $system = "jwc_cpu" ] ; then
    slurmset+=" --partition=batch"
elif [ $system = "jwc_gpu" ] ; then
    slurmset+=" --partition=gpus"
elif [ $system = "jrc_cpu" ] ; then
    slurmset+=" --partition=dc-cpu"
elif [ $system = "jrc_gpu" ] ; then
    slurmset+=" --partition=dc-gpu"
elif [ $system = "jrc_gh200" ] ; then
    slurmset+=" --partition=dc-gh"
elif [ $system = "jrc_h100" ] ; then
    slurmset+=" --partition=dc-h100"
else
    echo "error: system \"$system\" unknown!"
    exit
fi
sbatch $slurmset -v --wait ./batch.sh $system $trac $libs $compiler $gpu $npmin $npmax $npfac $meteo $rng $sort $phys $cache
