#! /bin/bash

# Check arguments...
if [ $# -ne 6 ] ; then
    cat <<EOF
ERA5 Reanalysis Data Download

Usage:
  $0 <year> <month> <day> <dir> <res> <dt>

Description:
  This script retrieves ECMWF's ERA5 reanalysis data
  (https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5/)
  and converts it for use with MPTRAC. The data are accessed through
  the Climate Data Store (CDS): https://cds.climate.copernicus.eu/

  The script retrieves data for the specified date (<year>, <month>, <day>)
  and stores the converted NetCDF files in the output directory (<dir>).

  Data are provided on 137 model levels and a regular
  latitude-longitude grid with horizontal resolution <res> (in degrees).
  The temporal resolution is specified via <dt> (in hours).

Requirements:
  - Python and the CDS API (https://cds.climate.copernicus.eu/api-how-to)
    must be installed and configured to retrieve data from the CDS.
  - Climate Data Operators (CDO, https://code.mpimet.mpg.de/projects/cdo)
    are required to convert GRIB files to NetCDF format for MPTRAC.

License:
  ERA5 data are distributed under the Copernicus open data license.
EOF
    exit
fi

# Set variables...
year=$1
mon=$(echo "$2" | awk '{printf("%02d", $1)}')
day=$(echo "$3" | awk '{printf("%02d", $1)}')
dir=$4/${year}/${mon}
res=$5
dt=$6

# Create Python environment...
python3 -m venv venv
source ./venv/bin/activate
pip install "cdsapi>=0.7.2"

# Set time string...
tstr=$(seq 0 "$dt" 23 | awk '{printf("%02d:00:00/", $1)}' | awk '{print substr($1,1,length($1)-1)}')

# Create directories...
mkdir -p "$dir" || exit
tmp=$(mktemp -d tmp.XXXXXXXX)

# Download data...
source ./venv/bin/activate
cd "$tmp" || exit
python3 <<EOF
import cdsapi
c = cdsapi.Client()

"""
ECMWF Parameter Database:
https://codes.ecmwf.int/grib/param-db/

Surface variables (analysis):
   34.128 = Sea surface temperature [K]
   59.128 = Convective available potential energy [J kg**-1]
  129.128 = Geopotential [m**2 s**-2]
  134.128 = Surface pressure [Pa]
  159.128 = Boundary layer height [m]
  160.128 = Standard deviation of sub-gridscale orography [m]
  164.128 = Total cloud cover (0 - 1)
  165.128 = 10 metre U wind component [m s**-1]
  166.128 = 10 metre V wind component [m s**-1]
  167.128 = 2 metre temperature [K]
  168.128 = 2 metre dewpoint temperature [K]
  172.128 = Land-sea mask (0 - 1)
  186.128 = Low cloud cover (0 - 1)
  187.128 = Medium cloud cover (0 - 1)
  188.128 = High cloud cover (0 - 1)
  229.128 = Instantaneous eastward turbulent surface stress [N m**-2]
  230.128 = Instantaneous northward turbulent surface stress [N m**-2]
  231.128 = Instantaneous surface sensible heat net flux [W m**-2]

# Surface variables (only in forecasts):
#   50.128 = Large-scale precipitation fraction [s]
#  142.128 = Large-scale precipitation [m]
#  143.128 = Convective precipitation [m]
#  146.128 = Time-integrated surface sensible heat net flux [J m**-2]
#  147.128 = Time-integrated surface latent heat net flux [J m**-2]
#  228.128 = Total precipitation [m]

Model levels:
   75 = Specific rain water content [kg kg**-1]
   76 = Specific snow water content [kg kg**-1]
  130 = Temperature [K]
  131 = U component of wind [m s**-1]
  132 = V component of wind [m s**-1]
  133 = Specific humidity [kg kg**-1]
  135 = Vertical velocity [Pa s**-1]
  203 = Ozone mass mixing ratio [kg kg**-1]
  246 = Specific cloud liquid water content [kg kg**-1]
  247 = Specific cloud ice water content [kg kg**-1]
  248 = Fraction of cloud cover (0 - 1)

# Model levels (optional):
#   77 = Eta-coordinate vertical velocity [s**-1]
"""

c.retrieve("reanalysis-era5-complete", {
    "class": "ea",
    "date": "${year}-${mon}-${day}",
    "expver": "1",
    "levtype": "sfc",
    "param": "34.128/59.128/129.128/134.128/159.128/160.128/164.128/165.128/166.128/167.128/168.128/172.128/186.128/187.128/188.128/229.128/230.128/231.128",
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
cdo -f nc copy "$tmp/ml.grib" "$tmp/ml.nc" && cdo -f nc copy "$tmp/sf.grib" "$tmp/sf.nc" || exit

# Merge netCDF files...
cdo merge "$tmp/sf.nc" "$tmp/ml.nc" "$tmp/all.nc" || exit

# Calculate pressure...
cdo -pressure_full "$tmp/all.nc" "$tmp/p.nc" && cdo merge "$tmp/all.nc" "$tmp/p.nc" "$tmp/allp.nc" || exit

# Split timesteps...
n=0
for t in $(seq 0 "$dt" 23) ; do
    let n=$n+1
    hour=$(echo "$t" | awk '{printf("%02d", $1)}')
    cdo seltimestep,$n "$tmp/allp.nc" "$dir/era5_${year}_${mon}_${day}_${hour}.nc" || exit
done

# Remove temporary directory...
rm -rf "$tmp"
