#! /bin/bash

# Check arguments...
if [ $# -ne 2 ] ; then
    cat <<EOF
NCEP–NCAR Reanalysis 1 Data Download

Usage:
  $0 <year> <dir>

Description:
  This script downloads 6-hourly reanalysis data from the NCEP–NCAR
  Reanalysis 1 project and converts them for use with MPTRAC.

  The data are retrieved for the specified year (<year>) and saved in
  the specified output directory (<dir>).

  Data source: NOAA PSL
  https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html

Requirements:
  - wget (https://www.gnu.org/software/wget) for data retrieval
  - Climate Data Operators (CDO, https://code.mpimet.mpg.de/projects/cdo)
    for data processing
  - NetCDF Operators (NCO, https://nco.sourceforge.net) for format conversion

License:
  The NCEP–NCAR Reanalysis 1 data are provided by NOAA and are in the public domain.
EOF
    exit
fi

# Set variables...
trac=../../src
year="$1"
dir="$2/$year"

# Create directories...
rm -rf "$dir" && mkdir -p "$dir" || exit
tmp=$(mktemp -d tmp.XXXXXXXX)

# Download surface data...
cd "$tmp" || exit
wget https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/surface/hgt.sfc.nc
wget https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/surface/pres.sfc."$year".nc || exit

# Download pressure level data...
for var in air omega shum uwnd vwnd ; do
    wget https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/pressure/"$var"."$year".nc || exit    
done
cd - || exit

# Loop over time steps...
for ts in $(seq 0 1 1463) ; do
    
    # Remove temporary files...
    rm -f "$tmp"/tmp*nc
    
    # Get 3D data...
    for var in air omega shum uwnd vwnd ; do
	ncks -d time,"$ts","$ts" "$tmp/$var.$year.nc" "$tmp/tmp.$var.nc" || continue
	ncwa -O -a time -d time,0,0 "$tmp/tmp.$var.nc" "$tmp/tmp.mean.$var.nc" || exit
	cdo intlevelx,1000,925,850,700,600,500,400,300,250,200,150,100,70,50,30,20,10 "$tmp/tmp.mean.$var.nc" "$tmp/tmp.ip.$var.nc" || exit
    done
    
    # Get surface pressure...
    ncks -d time,"$ts" "$tmp/pres.sfc.$year.nc" "$tmp/tmp.pres.sfc.nc" || exit
    ncwa -O -a time -d time,0,0 "$tmp/tmp.pres.sfc.nc" "$tmp/tmp.mean.pres.sfc.nc" || exit
    
    # Get geopotential height...
    ncks -d time,0 "$tmp/hgt.sfc.nc" "$tmp/tmp.hgt.sfc.nc" || exit
    ncwa -O -a time -d "$tmp"/time,0,0 "$tmp"/tmp.hgt.sfc.nc "$tmp"/tmp.mean.hgt.sfc.nc || exit
    
    # Merge and rename variables...
    cdo -f nc merge "$tmp"/tmp.ip.*.nc "$tmp/tmp.mean.pres.sfc.nc" "$tmp/tmp.mean.hgt.sfc.nc" "$tmp/tmp.merge.nc" || exit
    ncrename -d level,lev -v level,lev "$tmp/tmp.merge.nc" || exit
    ncrename -v air,T -v uwnd,U -v vwnd,V -v omega,W -v shum,Q -v pres,PS -v hgt,ZM "$tmp/tmp.merge.nc" || exit
    
    # Get time...
    jsec=$(ncks --trd -v time "$tmp/tmp.air.nc" | tr "=" " " | awk '{if($1=="time[0]") printf("%.0f", ($2-1753152)*3600.)}')
    tstr=$($trac/jsec2time "$jsec" | awk '{printf("%d_%02d_%02d_%02d", $1, $2, $3, $4)}')
    
    # Convert mbar to Pa...
    cdo chlevel,1000,100000,925,92500,850,85000,700,70000,600,60000,500,50000,400,40000,300,30000,250,25000,200,20000,150,15000,100,10000,70,7000,50,5000,30,3000,20,2000,10,1000 "$tmp/tmp.merge.nc" "$dir/ncep_$tstr.nc" || exit
    ncatted -a units,lev,m,c,"Pa" "$dir/ncep_$tstr.nc" || exit
    
    # Remove temporary files...
    rm -f "$tmp"/tmp*nc
    
done

# Remove temporary directory...
rm -rf "$tmp"
