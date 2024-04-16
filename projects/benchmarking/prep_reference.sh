#! /bin/bash

# Prepare reference data...
rm -rf reference && mkdir -p reference || exit
./benchmark.sh jwb exaww 03:00:00 avail avail gcc 1 1 48 0 100 100000000 10 era5 0 0 traj 0
rsync -av data/ reference/
./benchmark.sh jwb exaww 03:00:00 avail avail gcc 1 1 48 0 100 100000000 10 era5 0 0 full 0
rsync -av data/ reference/
./benchmark.sh jwb exaww 03:00:00 avail avail gcc 1 1 48 0 100 100000000 10 erai 0 0 traj 0
rsync -av data/ reference/
./benchmark.sh jwb exaww 03:00:00 avail avail gcc 1 1 48 0 100 100000000 10 erai 0 0 full 0
rsync -av data/ reference/

# Remove temporary files...
rm $(find reference -name "atm_init.tab") \
   $(find reference -name "dist.tab") \
   $(find reference -name "dirlist")

# Copy data...
refdir=/p/largedata/slmet/slmet111/datapub/hoffmann1/mptrac/projects/benchmarking/reference
rm -rf $refdir && mkdir -p $refdir && rsync -av reference/ $refdir/
