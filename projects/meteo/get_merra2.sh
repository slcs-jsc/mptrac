#! /bin/bash

# Check arguments...
if [ $# -ne 4 ] ; then
    cat <<EOF
MERRA-2 Data Download

Usage:
  $0 <year> <month> <day> <dir>

Description:
  This script retrieves instantaneous 3-dimensional, 3-hourly MERRA-2
  data from NASA's GES DISC archive
  (https://doi.org/10.5067/WWQSXQ8IVFW8) and converts it for use
  with MPTRAC.

  It downloads the data for the specified date (<year>, <month>, <day>)
  and stores the converted NetCDF files in the specified output
  directory (<dir>).

  Note: The data are provided on model levels. You will need to
  interpolate to pressure levels using MPTRAC's internal interpolation
  (via the 'MET_NP' and 'MET_P' control parameters).

Requirements:
  - wget (https://www.gnu.org/software/wget) to retrieve data.
  - NetCDF Operators (NCO, https://nco.sourceforge.net) for file conversion.
  - A NASA EarthData account (https://urs.earthdata.nasa.gov) linked to GES DISC.
  - Configuration of '~/.netrc' and '~/.urs_cookies' for authenticated wget access:
    See https://disc.gsfc.nasa.gov/data-access#windows_wget for setup instructions.

License:
  MERRA-2 data are provided by NASA under the Earth Science Data and Information System (ESDIS) policy.
EOF
    exit
fi

# Set variables...
year="$1"
mon=$(echo "$2" | awk '{printf("%02d", $1)}')
day=$(echo "$3" | awk '{printf("%02d", $1)}')
dir="$4/$year"

# Create directories...
mkdir -p "$dir" || exit
tmp=$(mktemp -d tmp.XXXXXXXX)

# Download data...
cd "$tmp" || exit
vers=400
[ "$year" -le 2010 ] && vers=300
[ "$year" -le 2000 ] && vers=200
[ "$year" -le 1991 ] && vers=100
wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies "https://data.gesdisc.earthdata.nasa.gov/data/MERRA2/M2I3NVASM.5.12.4/$year/$mon/MERRA2_${vers}.inst3_3d_asm_Nv.${year}${mon}${day}.nc4" || exit
cd - || exit

# Convert data...
for f in $(ls "$tmp"/MERRA2*"${year}""${mon}"*nc4) ; do
    
    # Write info...
    echo "Convert $f ..."
    
    # Get time string...
    tstr=$(basename "$f" .nc4 | awk '{print substr($1,length($1)-7,4)"_"substr($1,length($1)-3,2)"_"substr($1,length($1)-1,2)}')
    
    # Loop over hours...
    for h in $(seq 0 1 7) ; do
	
	# Set hour string...
	hstr=$(echo "$h" | awk '{printf("%02d", 3*$1)}')
	
	# Extract data...
	ncks -3 -O -v PS,PHIS,PL,T,U,V,OMEGA,O3,QV,QL,QI -d time,"$h" "$f" "merra2_${tstr}_${hstr}.nc" || exit
	
	# Rename variables...
	ncrename -v PHIS,Z -v OMEGA,W -v QV,Q -v QL,CLWC -v QI,CIWC "merra2_${tstr}_${hstr}.nc" || exit
	
	# Move file...
	mv "merra2_${tstr}_${hstr}.nc" "$dir" || exit
	
    done
done

# Remove temporary directory...
rm -rf "$tmp"
