#! /bin/bash

# Check arguments...
if [ $# -ne 12 ] ; then
    echo "usage: $0 <system> <libs> <compiler> <gpu> <npmin> <npmax> <npfac> <meteo> <rng> <sort> <phys> <cache>"
    exit
fi

# Set variables...
trac=../../
system=$1
libs=$2
compiler=$3
gpu=$4
npmin=$5
npmax=$6
npfac=$7
meteo=$8
rng=$9
sort=${10}
phys=${11}
cache=${12}

# Write info:
echo "job started: $(date)"
echo -e "\nuname: $(uname -a)"

# Compile libraries...
echo -e "\nCompile libraries..."
if [ "$libs" = "compute" ] ; then
    cd $trac/libs && ./build.sh -a ; cd -
fi
du -h $trac/libs/build/lib/ || exit

# Load modules...
echo -e "\nLoad modules..."
ml purge
if [ "$compiler" = "gcc" ] ; then
    ml GCC ParaStationMPI
    echo -e "\ngcc version:" ; gcc --version
elif [ "$compiler" = "nvc" ] ; then
    ml NVHPC ParaStationMPI
    echo -e "\nnvc version:" ; nvc --version
else
    echo "error: compiler \"$compiler\" unknown!"
    exit
fi
ml git
ml list

# Write info...
echo -e "\nMPI information..."
mpicc -v

# Set environment variables...
echo -e "\nSet environment variables..."
export LD_LIBRARY_PATH=$trac/libs/build/lib/:$LD_LIBRARY_PATH
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=true
export OMP_PLACES=sockets
#export OMP_DISPLAY_ENV=true
#export OMP_WAIT_POLICY=active
#export OMP_DYNAMIC=false
#export NV_ACC_NOTIFY=31
#export NV_ACC_TIME=1
	
# Create data directory...
rm -rf data && mkdir -p data || exit

# Loop over number of particles...
for np in $(echo "$npmin" "$npmax" "$npfac" | awk '{for(np=$1; np<=$2; np*=$3) print np}') ; do
    
    # Set logfile...
    mkdir -p logs || exit
    log=log.$system.$compiler.nodes_${SLURM_JOB_NUM_NODES}.tasks_${SLURM_NTASKS_PER_NODE}.threads_${SLURM_CPUS_PER_TASK}.gpu_$gpu.meteo_$meteo.rng_$rng.sort_$sort.np_$np.phys_${phys}.cache_${cache}
    echo -e "\nProcessing $log ..."
    {
	
	# Set compile flags...
	echo -e "\nCompile MPTRAC..."
 	np_comp=$(echo "$np" | awk '{if($1<100) print 100; else print $1}')
 	[ "$meteo" = "erai" ] && defs="-DNP=$np_comp -DNQ=8 -DEX=482 -DEY=242 -DEP=62"
 	[ "$meteo" = "era5" ] && defs="-DNP=$np_comp -DNQ=8 -DEX=1202 -DEY=602 -DEP=140"
 	flags="COMPILER=$compiler MPI=1 STATIC=0 "
	[ "$gpu" = "1" ] && flags+=" GPU=1"
	[ "$gpu" = "2" ] && flags+=" GPU=1 GPU_PIN=1"
 	[ "$rng" = "2" ] && flags+=" CURAND=1"
	[ "$sort" != "0" ] && flags+=" THRUST=1"
	[ "$cache" = "2" ] && flags+=" ASYNCIO=1"
	
	# Compile...
	cd $trac/src && make clean && make -j DEFINES="$defs" "$flags" || exit
	cd -
	
	# MPTRAC setup...
	qnt="NQ 8 QNT_NAME[0] theta QNT_NAME[1] pv QNT_NAME[2] h2o QNT_NAME[3] o3 QNT_NAME[4] Cn2o QNT_NAME[5] Cccl4 QNT_NAME[6] Cccl3f QNT_NAME[7] Cccl2f2"
	[ "$meteo" = "erai" ] && metbase="METBASE meteo/erai_pck/ei DT_MET 21600 DT_MOD 360 MET_TYPE 2"
	[ "$meteo" = "era5" ] && metbase="METBASE meteo/era5_pck/era5 DT_MET 3600 DT_MOD 180 MET_TYPE 2"
	[ "$cache" = "1" ] && metbase+=" MET_CACHE 1"
	if [ "$phys" = "full" ] ; then
	    param="CONV_CAPE 0.0 BOUND_P0 1e100 BOUND_P1 -1e100 BOUND_LAT0 -90 BOUND_LAT1 90 BOUND_DPS 150 MIXING_TROP 0.1 MIXING_STRAT 0.1"
	else
	    param="TURB_DX_TROP 0 TURB_DX_STRAT 0 TURB_DZ_TROP 0 TURB_DZ_STRAT 0 TURB_MESOX 0 TURB_MESOZ 0"
	fi
	t0=$($trac/src/time2jsec 2017 1 1 0 0 0 0)
	t1=$($trac/src/time2jsec 2017 1 2 0 0 0 0)
	
	# Create data directory...
	dir=data/$meteo/phys_${phys}/$np
	rm -rf "$dir" && mkdir -p "$dir" || exit
	
	# Create init file...
	$trac/src/atm_init - $dir/atm_init.tab "$qnt" ATM_TYPE 1 \
			   INIT_T0 "$t0" INIT_T1 "$t0" \
			   INIT_ULON 360 INIT_ULAT 180 \
			   INIT_Z0 30 INIT_Z1 30 INIT_UZ 60 \
			   INIT_EVENLY 1 INIT_REP "$np" INIT_MASS 1e9 || exit
	
	# Create data directories for tasks...
	for task in $(seq ${SLURM_NTASKS}) ; do
	    mkdir -p $dir/$task && cp $dir/atm_init.tab $dir/$task/ || exit
	    echo "$dir/$task" >> $dir/dirlist
	done
	
	# Calculate trajectories...
	srun $trac/src/trac $dir/dirlist - atm_init.tab "$qnt" T_STOP "$t1" \
	     $metbase MET_DT_OUT 3600 $param RNG_TYPE $rng SORT_DT $sort \
	     ATM_BASENAME atm ATM_DT_OUT 21600 ATM_TYPE 1 \
	     GRID_BASENAME grid GRID_DT_OUT 21600 \
	     GRID_NX 36 GRID_NY 36 GRID_Z0 -2 GRID_Z1 62 GRID_NZ 64
	
	# Calculate transport deviations...
	for var in mean stddev absdev median mad ; do
	    $trac/src/atm_dist - $dir/dist.tab $var \
			       $dir/1/atm_2017_01_01_00_00.bin reference/$meteo/phys_${phys}/$np/1/atm_2017_01_01_00_00.bin \
			       $dir/1/atm_2017_01_01_06_00.bin reference/$meteo/phys_${phys}/$np/1/atm_2017_01_01_06_00.bin \
			       $dir/1/atm_2017_01_01_12_00.bin reference/$meteo/phys_${phys}/$np/1/atm_2017_01_01_12_00.bin \
			       $dir/1/atm_2017_01_01_18_00.bin reference/$meteo/phys_${phys}/$np/1/atm_2017_01_01_18_00.bin \
			       $dir/1/atm_2017_01_02_00_00.bin reference/$meteo/phys_${phys}/$np/1/atm_2017_01_02_00_00.bin \
			       $qnt ATM_TYPE 1
	    echo -e "\nTransport deviations: $var ..."
	    cat $dir/dist.tab
	done
	
	# Compare grid files...
	echo -e "\nGrid differences..."
	paste $dir/1/grid_2017_01_02_00_00.tab reference/$meteo/phys_${phys}/$np/1/grid_2017_01_02_00_00.tab | awk '{
	  if(NF>0 && $1!="#") {
            cmax=NF/2
            for(c=1; c<=cmax; c++) {
              diff=1.0*$(c)-1.0*$(c+cmax)
              mean=0.5*$(c)+0.5*$(c+cmax)
              bias_abs[c]+=diff
              rmse_abs[c]+=diff^2
              if(mean!=0) {
                bias_rel[c]+=diff/mean
                rmse_rel[c]+=(diff/mean)^2
              }
              n[c]++
            }
          }
        }END{
          for(c=1; c<=cmax; c++)
            print "  column=", c, "| bias_abs=", bias_abs[c]/n[c], "| bias_rel=", 100.*bias_rel[c]/n[c], "% | rmse_abs=", sqrt(rmse_abs[c]/n[c]), "| rmse_rel=", 100.*sqrt(rmse_rel[c]/n[c]), "% | n=", n[c]
        }'
	
    } > logs/$log 2>&1
    
done

# Write info:
echo -e "\njob finished: $(date)"
