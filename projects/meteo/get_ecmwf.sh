#! /bin/bash
set -euo pipefail

# Check arguments...
if [ $# -ne 5 ] ; then
    cat <<EOF
ECMWF forecast data download

usage: $0 <year> <month> <day> <dir> <ifs|aifs>

This script retrieves ECMWF's IFS or AIFS forecast data
(https://www.ecmwf.int/en/forecasts/dataset/) and converts them for
use with MPTRAC. The data are retrieved from Amazon AWS.

This script retrieves the forecast for 00z for a given day (<year>,
<month>, <day>) and places the converted data files in the archive
directory (<dir>). You can set <year>, <mon>, <day> to "today"
to retrieve the most recent dara.

The data are retrieved on 13 pressure levels and a
regular longitude-latitude grid with grid box size of 0.25 degrees.
The time step is 3 hourly until 144 h and 6 hourly until 360 h.

The ECMWF open data are provided under the CC BY-4 license. The
Climate Data Operators (CDO, https://code.mpimet.mpg.de/projects/cdo)
are used to convert the netCDF files for MPTRAC.
EOF
    exit
fi

# Set variables...
trac=../../src
year=$1
mon=$(echo "$2" | awk '{printf("%02d", $1)}')
day=$(echo "$3" | awk '{printf("%02d", $1)}')
dir=$4
model=$5

# Determine date...
date=${year}${mon}${day}
hour="00"
if [ "$year" = "today" ] ; then
    date=$(date -u -d "${hour} hour ago" +%Y%m%d)
    year=$(echo $date | awk '{print substr($1,1,4)}')
    mon=$(echo $date | awk '{print substr($1,5,2)}')
    day=$(echo $date | awk '{print substr($1,7,2)}')
fi

# Create directories...
outdir="${dir}/${model}/${date}/${hour}z"
rm -rf tmp && mkdir -p tmp "$outdir" || exit

# Loop over forecast steps...
for step in $(seq 0 3 144) $(seq 150 6 360) ; do
    
    # Write info...
    echo "Downloading HRES forecast: ${date}${hour}Z + ${step}h"
    
    # Download...
    root="https://ecmwf-forecasts.s3.eu-central-1.amazonaws.com"
    url="${root}/${date}/${hour}z/${model}/0p25/oper/${date}${hour}0000-${step}h-oper-fc.grib2"
    wget -c --show-progress -O "tmp/${step}.grib2" "$url" || exit
    
    # Get time...
    t2=$($trac/time2jsec $year $mon $day $hour 0 0 0 | awk -v step=$step '{printf("%.2f", $1 + step * 3600.)}')
    year2=$($trac/jsec2time $t2 | awk '{print $1}')
    mon2=$($trac/jsec2time $t2 | awk '{printf("%02d", $2)}')
    day2=$($trac/jsec2time $t2 | awk '{printf("%02d", $3)}')
    hour2=$($trac/jsec2time $t2 | awk '{printf("%02d", $4)}')
    
    # Convert to netCDF...
    cdo -f nc4 copy tmp/${step}.grib2 $outdir/${model}_${year2}_${mon2}_${day2}_${hour2}.nc || exit
    
done

# Remove temporary data...
rm -rf tmp || exit

# Write info...
echo "Download complete."
