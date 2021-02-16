#! /bin/bash

# Setup...
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../libs/build/lib
trac=../src
if [ $# -eq 0 ] ; then
    metfile=meteo/ei_2011_06_05_00.nc
    p0=240
else
    if [ $# -ne 2 ] ; then
	echo "usage: $0 <met.nc> <plev>"
	exit
    fi
    metfile=$1
    p0=$2
fi
tstr=$(basename $metfile | awk 'BEGIN{FS="_"}{printf("%d-%02d-%02d, %02d:00 UTC", $2, $3, $4, $5)}')
tstr2=$(basename $metfile | awk 'BEGIN{FS="_"}{printf("%d_%02d_%02d_%02d", $2, $3, $4, $5)}')
z0=$(echo $p0 | awk '{print 7.*log(1013.25/$1)}')

# Create directories...
mkdir -p data plots

# Sample meteo data...
$trac/met_map - data/map_${tstr2}_${p0}hPa.tab $metfile MAP_Z0 $z0

# Plot meteo map data...
gnuplot <<EOF
set term png truecolor crop linewidth 2 font "Helvetica" 24 size 1440,900
set size ratio 0.75
#set pm3d map interp 2,2
set pm3d map
set pal def
set xla "longitude [deg]"
set yla "latitude [deg]"
set xtics 60
set ytics 30
set mxtics 6
set mytics 5
set xra [-180:180]
set yra [-90:90]
set title "ERA-Interim | $tstr | $p0 Pa"

set out "plots/map_${tstr2}_${p0}hPa_t.png"
set cbla "temperature [K]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$6) t "", "wcl.tab" u 1:2:(250) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_u.png"
set cbla "zonal wind [m/s]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$7) t "", "wcl.tab" u 1:2:(0) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_v.png"
set cbla "meridional wind [m/s]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$8) t "", "wcl.tab" u 1:2:(0) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_w.png"
set cbla "vertical velocity [hPa/s]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$9) t "", "wcl.tab" u 1:2:(0) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_h2o.png"
set cbla "water vapor vmr [ppm]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1e6*\$10) t "", "wcl.tab" u 1:2:(0) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_o3.png"
set cbla "ozone vmr [ppb]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1e9*\$11) t "", "wcl.tab" u 1:2:(0) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_z.png"
set cbla "geopotential height [km]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$12) t "", "wcl.tab" u 1:2:(10) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_pv.png"
set cbla "potential vortivity [PVU]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$13) t "", "wcl.tab" u 1:2:(0) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_lwc.png"
set cbla "cloud liquid water content [g/kg]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1e3*\$23) t "", "wcl.tab" u 1:2:(0) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_iwc.png"
set cbla "cloud ice water content [g/kg]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1e3*\$24) t "", "wcl.tab" u 1:2:(0) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_rh.png"
set cbla "relative humidity [%]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$31) t "", "wcl.tab" u 1:2:(0) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_rhice.png"
set cbla "relative humidity over ice [%]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$32) t "", "wcl.tab" u 1:2:(0) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_tdew.png"
set cbla "dew point temperature [K]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$33) t "", "wcl.tab" u 1:2:(230) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_tice.png"
set cbla "frost point temperature [K]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$34) t "", "wcl.tab" u 1:2:(230) w l lt -1 lw .5 lc rgbcolor "gray" t ""

set title "ERA-Interim | $tstr | surface"
set out "plots/map_${tstr2}_${p0}hPa_ps.png"
set cbla "pressure [hPa]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$14) t "", "wcl.tab" u 1:2:(1000) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_ts.png"
set cbla "temperature [K]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$15) t "", "wcl.tab" u 1:2:(270) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_zs.png"
set cbla "geopotential height [km]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$16) t "", "wcl.tab" u 1:2:(0) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_us.png"
set cbla "zonal wind [m/s]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$17) t "", "wcl.tab" u 1:2:(0) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_vs.png"
set cbla "meridional wind [m/s]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$18) t "", "wcl.tab" u 1:2:(0) w l lt -1 lw .5 lc rgbcolor "gray" t ""

set title "ERA-Interim | $tstr | tropopause"
set out "plots/map_${tstr2}_${p0}hPa_pt.png"
set cbla "pressure [hPa]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$19) t "", "wcl.tab" u 1:2:(300) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_zt.png"
set cbla "geopotential height [km]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$20) t "", "wcl.tab" u 1:2:(10) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_tt.png"
set cbla "temperature [K]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$21) t "", "wcl.tab" u 1:2:(220) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_h2ot.png"
set cbla "water vapor vmr [ppm]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1e6*\$22) t "", "wcl.tab" u 1:2:(0) w l lt -1 lw .5 lc rgbcolor "gray" t ""

set title "ERA-Interim | $tstr"
set out "plots/map_${tstr2}_${p0}hPa_tccw.png"
set cbla "total column cloud water [kg/m^2]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$25) t "", "wcl.tab" u 1:2:(0) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_pc.png"
set cbla "cloud top pressure [hPa]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$26) t "", "wcl.tab" u 1:2:(300) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_plcl.png"
set cbla "pressure at lifted condensation level [hPa]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$27) t "", "wcl.tab" u 1:2:(300) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_plfc.png"
set cbla "pressure at level of free convection [hPa]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$28) t "", "wcl.tab" u 1:2:(300) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_pel.png"
set cbla "pressure at equilibrium level [hPa]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$29) t "", "wcl.tab" u 1:2:(300) w l lt -1 lw .5 lc rgbcolor "gray" t ""
set out "plots/map_${tstr2}_${p0}hPa_cape.png"
set cbla "convective available potential energy [J/kg]"
splot "data/map_${tstr2}_${p0}hPa.tab" u 3:4:(1.*\$30) t "", "wcl.tab" u 1:2:(300) w l lt -1 lw .5 lc rgbcolor "gray" t ""
EOF

# Check deviations...
echo -e "\nDifferences wrt to reference data:"
diff -s data/map_${tstr2}_${p0}hPa.tab data.org/map_${tstr2}_${p0}hPa.tab
