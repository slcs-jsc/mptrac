#! /bin/bash
set -euo pipefail

# Check arguments...
if [ $# -ne 5 ] ; then
    cat <<EOF
ECMWF Forecast Data Download

Usage:
  $0 <year> <month> <day> <dir> <ifs|aifs-single>

Description:
  This script retrieves ECMWF IFS or AIFS forecast data
  (https://www.ecmwf.int/en/forecasts/dataset/) and converts it for
  use with MPTRAC. The data are downloaded from Amazon AWS.

  The script retrieves the forecast for 00Z on the specified date
  (<year>, <month>, <day>) and stores the converted data files in the
  output directory <dir>. You can set <year>, <month>, and <day> to
  "today" to retrieve the most recent available data.

  The forecast data include 13 pressure levels on a regular
  longitude-latitude grid with a spatial resolution of 0.25 degrees.
  Time steps are 3-hourly up to 144 hours and 6-hourly up to 360 hours
  for IFS. AIFS has 6-hourly steps up to 360 hours.

Requirements:
  - wget (https://www.gnu.org/software/wget) for data retrieval
  - Climate Data Operators (CDO, https://code.mpimet.mpg.de/projects/cdo)
    for data processing

License:
  - The ECMWF open data are provided under the Creative Commons BY 4.0 license.
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

# Set step sequence...
if [ "$model" = "ifs" ]; then
    step_seq="$(seq 0 3 144) $(seq 150 6 360)"
elif [ "$model" = "aifs-single" ]; then
    step_seq=$(seq 0 6 360)
else
    echo "Cannot select model $model!"
    exit
fi

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
outdir=${dir}/"${model}_${year}_${mon}_${day}"
tmp=$(mktemp -d tmp.XXXXXXXX)
rm -rf $tmp && mkdir -p $tmp "$outdir" || exit

# Loop over forecast steps...
for step in $step_seq ; do
    
    # Write info...
    echo "Downloading HRES forecast: ${date}${hour}Z + ${step}h"
    
    # Download...
    root="https://ecmwf-forecasts.s3.eu-central-1.amazonaws.com"
    url="${root}/${date}/${hour}z/${model}/0p25/oper/${date}${hour}0000-${step}h-oper-fc.grib2"
    wget -c --show-progress -O "$tmp/${step}.grib2" "$url" || exit
    
    # Get time...
    t2=$($trac/time2jsec $year $mon $day $hour 0 0 0 | awk -v step=$step '{printf("%.2f", $1 + step * 3600.)}')
    year2=$($trac/jsec2time $t2 | awk '{print $1}')
    mon2=$($trac/jsec2time $t2 | awk '{printf("%02d", $2)}')
    day2=$($trac/jsec2time $t2 | awk '{printf("%02d", $3)}')
    hour2=$($trac/jsec2time $t2 | awk '{printf("%02d", $4)}')
    
    # Convert to netCDF...
    outfile="$outdir/${model}_${year2}_${mon2}_${day2}_${hour2}.nc"
    cdo -f nc4 copy "$tmp/${step}.grib2" "$outfile" || exit
    
    # Copy geopotential...
    if [ "$step" -eq 0 ]; then
	
	# Fix for AIFS data...
	if [ "$model" = "aifs-single" ] ; then
	    
	    # Rename z_2 to z and delete original z...
	    tmpfile="$tmp/tmp_${model}_${year2}_${mon2}_${day2}_${hour2}.nc"
            cdo -f nc4 chname,z_2,z -delvar,z "$outfile" "$tmpfile" || exit
            mv "$tmpfile" "$outfile" || exit
	fi
	
        # Save clean z to zfile for future merging...
	zfile="$tmp/z_${model}_${year2}_${mon2}_${day2}_00.nc"
        cdo -f nc4 -selvar,z "$outfile" "$zfile" || exit
	
    elif [ "$step" -gt 0 ] ; then
	
        # Remove existing z from current file...
	cleanedfile="$tmp/cleaned_${model}_${year2}_${mon2}_${day2}_${hour2}.nc"
        cdo -f nc4 -delvar,z "$outfile" "$cleanedfile" || exit
	
        # Merge cleaned file with step-0 z...
	tmpfile="$tmp/tmp_${model}_${year2}_${mon2}_${day2}_${hour2}.nc"
        cdo -f nc4 merge "$zfile" "$cleanedfile" "$tmpfile" || exit
        mv "$tmpfile" "$outfile" || exit
    fi
    
done

# Remove temporary data...
rm -rf $tmp || exit

# Write info...
echo "Download complete."
