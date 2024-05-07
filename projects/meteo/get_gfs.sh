#! /bin/bash

## Load modules...
#ml purge
#ml Stages/2022
#ml GCC ParaStationMPI CDO NCO

# Check arguments...
if [ $# -ne 4 ] ; then
    cat <<EOF
GFS Data Download

usage: $0 <year> <month> <day> <dir>

This script retrieves 3-hourly 10-day forecasts from NCEP's Global
Forecast System (GFS) at 0.5 degree horizontal resolution and converts
them for use with MPTRAC. The data are retrieved from NOAA's
Operational Model Archive and Distribution System (NOMADS,
https://nomads.ncep.noaa.gov/).

This script retrieves the GFS data for a given day (<year>, <month>,
<day>) and places the converted data files in the <dir> directory.

You will need wget (https://www.gnu.org/software/wget) to retrieve the
GFS data, and the Climate Data Operators (CDO,
https://code.mpimet.mpg.de/projects/cdo) and netCDF Operators (NCO,
https://nco.sourceforge.net) to convert the netCDF files for MPTRAC.
EOF
    exit
fi

# Set variables...
trac=../../src
year=$1
mon=$(echo $2 | awk '{printf("%02d", $1)}')
day=$(echo $3 | awk '{printf("%02d", $1)}')
dir=$4/data_${year}_${mon}_${day}

# Set time string...
tstr=$(echo $year $mon $day | awk '{printf("%04d%02d%02d%%2F00", $1, $2, $3)}')
echo "Get GFS analysis for $tstr..."

# Create data directory...
rm -rf $dir && mkdir -p $dir || exit

# Download and convert GFS forecast data...
for fh in $(seq 0 3 240 | awk '{printf(" %03d", $1)}') ; do
    
    # Remove temporary files...
    rm -f $dir/tmp.*
    
    # Download grib file...
    echo "  Download forecast for hour $fh..."
    wget -nv -O $dir/tmp.grib2 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p50.pl?file=gfs.t00z.pgrb2full.0p50.f'$fh'&&lev_0.01_mb=on&lev_0.02_mb=on&lev_0.04_mb=on&lev_0.07_mb=on&lev_0.1_mb=on&lev_0.2_mb=on&lev_0.4_mb=on&lev_0.7_mb=on&lev_1_mb=on&lev_2_mb=on&lev_3_mb=on&lev_5_mb=on&lev_7_mb=on&lev_10_mb=on&lev_15_mb=on&lev_20_mb=on&lev_30_mb=on&lev_40_mb=on&lev_50_mb=on&lev_70_mb=on&lev_100_mb=on&lev_150_mb=on&lev_200_mb=on&lev_250_mb=on&lev_300_mb=on&lev_350_mb=on&lev_400_mb=on&lev_450_mb=on&lev_500_mb=on&lev_550_mb=on&lev_600_mb=on&lev_650_mb=on&lev_700_mb=on&lev_750_mb=on&lev_800_mb=on&lev_850_mb=on&lev_900_mb=on&lev_925_mb=on&lev_950_mb=on&lev_975_mb=on&lev_1000_mb=on&lev_surface=on&var_TMP=on&var_UGRD=on&var_VGRD=on&var_VVEL=on&var_SPFH=on&var_O3MR=on&var_PRES=on&var_HGT=on&leftlon=0&rightlon=360&toplat=90&bottomlat=-90&dir=%2Fgfs.'$tstr'%2Fatmos' || continue
    
    # Set netCDF output filename...
    out=$($trac/jsec2time $($trac/time2jsec $year $mon $day 0 0 0 0 | awk -v fh=$fh '{print $1+3600*fh}') | awk '{printf("gfs_%d_%02d_%02d_%02d.nc", $1, $2, $3, $4)}')
    
    # Convert from grib to netCDF...
    cdo -f nc copy $dir/tmp.grib2 $dir/tmp.nc || exit
    ncrename -d plev,lev -v plev,lev -v orog,zm -v sp,ps -v o3mr,o3 $dir/tmp.nc || exit
    cdo invertlev $dir/tmp.nc $dir/$out || exit
    
    # Remove temporary files...
    rm -f $dir/tmp.*

done
