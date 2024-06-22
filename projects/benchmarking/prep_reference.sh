#! /bin/bash

# Prepare reference data...
rm -rf reference && mkdir -p reference || exit
for meteo in era5 erai ; do
    for phys in traj full ; do
	./benchmark.sh jwb exaww 03:00:00 ../.. avail gcc 1 1 48 0 100 100000000 10 $meteo meteo 0 0 $phys 0 reference
	rsync -av data/ reference/
    done
done

# Remove temporary files...
rm $(find reference -name "atm_init.tab") \
   $(find reference -name "dist.tab") \
   $(find reference -name "dirlist")

# Copy data...
refdir=/p/data1/slmet/datapub/hoffmann1/mptrac/projects/benchmarking/reference
rm -rf $refdir && mkdir -p $refdir && rsync -av reference/ $refdir/
