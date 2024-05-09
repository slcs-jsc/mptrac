#! /bin/bash

# Check arguments...
if [ $# -ne 6 ] ; then
    cat <<EOF
ERA5 data download

usage: $0 <year> <month> <day> <dir> <res> <dt>

This script retrieves ECMWF's ERA5 reanalysis data
(https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5/) and
converts them for use with MPTRAC. The data are retrieved from the
Climate Data Store (https://cds.climate.copernicus.eu/).

This script retrieves the data for a given day (<year>, <month>,
<day>) and places the converted data files in the archive directory
(<dir>). The data are retrieved on 137 model levels and a regular
longitude-latitude grid with grid box size <res> (in degrees). Also
specify the timestep <dt> (in hours) of the retrieved data.

You will need Python and the CDS API to access and retrieve data from
the CDS (https://cds.climate.copernicus.eu/api-how-to). The Climate
Data Operators (CDO, https://code.mpimet.mpg.de/projects/cdo) are used
to convert the netCDF files for MPTRAC.
EOF
    exit
fi

# Set variables...
trac=../../src
year=$1
mon=$(echo $2 | awk '{printf("%02d", $1)}')
day=$(echo $3 | awk '{printf("%02d", $1)}')
dir=$4/${year}/${mon}
res=$5
dt=$6

# Set time string...
tstr=$(seq 0 $dt 23 | awk '{printf("%02d:00:00/", $1)}' | awk '{print substr($1,1,length($1)-1)}')

# Create directories...
mkdir -p $dir || exit
tmp=$(mktemp -d tmp.XXXXXXXX)

# Download data...
cd $tmp || exit
python <<EOF
import cdsapi
c = cdsapi.Client()

c.retrieve("reanalysis-era5-complete", {
    "class": "ea",
    "date": "${year}-${mon}-${day}",
    "expver": "1",
    "levtype": "sfc",
    "param": "34.128/59.128/129.128/134.128/159.128/165.128/166.128/167.128/172.128",
    "stream": "oper",
    "time": "$tstr",
    "type": "an",
    "grid": "$res/$res"
}, "sf.grib")

c.retrieve("reanalysis-era5-complete", {
    "class": "ea",
    "date": "${year}-${mon}-${day}",
    "expver": "1",
    "levelist": "1/to/137",
    "levtype": "ml",
    "param": "75/76/130/131/132/133/135/203/246/247/248",
    "stream": "oper",
    "time": "$tstr",
    "type": "an",
    "grid": "$res/$res"
}, "ml.grib")
EOF
cd - || exit

# Convert to netCDF...
cdo -f nc copy $tmp/ml.grib $tmp/ml.nc && cdo -f nc copy $tmp/sf.grib $tmp/sf.nc || exit

# Merge netCDF files...
cdo merge $tmp/sf.nc $tmp/ml.nc $tmp/all.nc || exit

# Calculate pressure...
cdo -pressure_fl $tmp/all.nc $tmp/p.nc && cdo merge $tmp/all.nc $tmp/p.nc $tmp/allp.nc || exit

# Split timesteps...
let n=0
for t in $(seq 0 $dt 23) ; do
    let n=$n+1
    hour=$(echo $t | awk '{printf("%02d", $1)}')
    cdo seltimestep,$n $tmp/allp.nc $dir/era5_${year}_${mon}_${day}_${hour}.nc || exit
done

# Remove temporary directory...
rm -rf $tmp
