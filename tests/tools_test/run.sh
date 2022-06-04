#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4

# Setup...
trac=../../src

# Create directories...
rm -rf data && mkdir -p data

echo "time..."
(for year in 1900 1980 2000 2020 2100 ; do
    for mon in 1 7 12 ; do
	for day in 1 15 31 ; do
	    for hour in 0 12 24 ; do
		t0=$($trac/time2jsec $year $mon $day $hour 0 0 0)
		echo "$($trac/jsec2time $t0) = $t0"
		d0=$($trac/day2doy $year $mon $day)
		echo "$($trac/doy2day $d0) = $d0"
	    done
	done
    done
done) | tee data/time.tab

echo "sedi..."
(for p in 1000 100 10 1 ; do
    for t in 200 250 300 ; do
	for rp in 0.1 1 10 100 ; do
	    for rho in 500 1000 2000 ; do
		echo
		$trac/sedi $p $t $rp $rho
	    done
	done
    done
done) | tee data/sedi.tab

echo "tnat..."
(for p in 1000 100 10 1 ; do
    for h2o in 2e-6 4e-6 6e-6 ; do
	for hno3 in 1e-9 4e-9 9e-9 ; do
	    echo
	    $trac/tnat $p $h2o $hno3
	done
    done
done) | tee data/tnat.tab

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/$(basename "$f") "$f"
    [ $? -ne 0 ] && error=1
done
exit $error
