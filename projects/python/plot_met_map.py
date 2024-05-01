######################################################################
## Author: Farahnaz Khosrawi
## Data created: 14.09.2023
## Last modified: 11.01.2024
## Purpose: Plots MPTRAC met_map output - global fields (parameter vs
## longitude/latitude)
######################################################################

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import xarray as xr
import os

# MPTRAC met_map ERA5 ouput - 47 species are written out. Note newer MPTRAC version write 49 species out!
# $1 = time [s]
# $2 = altitude [km]
# $3 = longitude [deg]
# $4 = latitude [deg]
# $5 = pressure [hPa]
# $6 = temperature [K]
# $7 = zonal wind [m/s]
# $8 = meridional wind [m/s]
# $9 = vertical velocity [hPa/s]
# $10 = H2O volume mixing ratio [ppv]
# $11 = O3 volume mixing ratio [ppv]
# $12 = geopotential height [km]
# $13 = potential vorticity [PVU]
# $14 = surface pressure [hPa]
# $15 = surface temperature [K]
# $16 = surface geopotential height [km]
# $17 = surface zonal wind [m/s]
# $18 = surface meridional wind [m/s]
# $19 = land-sea mask [1]
# $20 = sea surface temperature [K]
# $21 = tropopause pressure [hPa]
# $22 = tropopause geopotential height [km]
# $23 = tropopause temperature [K]
# $24 = tropopause water vapor [ppv]
# $25 = cloud liquid water content [kg/kg]
# $26 = cloud ice water content [kg/kg]
# $27 = cloud cover [1]
# $28 = total column cloud water [kg/m^2]
# $29 = cloud top pressure [hPa]
# $30 = cloud bottom pressure [hPa]
# $31 = pressure at lifted condensation level (LCL) [hPa]
# $32 = pressure at level of free convection (LFC) [hPa]
# $33 = pressure at equilibrium level (EL) [hPa]
# $34 = convective available potential energy (CAPE) [J/kg]
# $35 = convective inhibition (CIN) [J/kg]
# $36 = relative humidity over water [%]
# $37 = relative humidity over ice [%]
# $38 = dew point temperature [K]
# $39 = frost point temperature [K]
# $40 = NAT temperature [K]
# $41 = HNO3 volume mixing ratio [ppv]
# $42 = OH concentration [molec/cm^3]
# $43 = H2O2 concentration [molec/cm^3]
# $44 = boundary layer pressure [hPa]
# $45 = number of data points
# $46 = number of tropopause data points
# $47 = number of CAPE data points

data = 'ERA5'
level = '10_km'    #10_km, 5_km, 2_km

file = 'data/map_era5_2017010817_2_2.tab'

header=['time', 'altitude', 'longitude', 'latitude', 'pressure', 'temperature', 'u', 'v', 'w', 'H2O', 'O3', 'gph', 'pv', 'psurf', 'tsurf', 'gph_surf' , 'u_surf', 'v_surf', 'land_sea_mask', 'SST', 'p_trop', 'gph_trop', 'temp_trop', 'H2O_trop', 'cloud_lw', 'cloud_iwc', 'cloud_cover', 'tot_col_water', 'cloud_top press', 'cloud_bottom_press', 'pressure_lifted_CL', 'press_lev_free_conv', 'press_EL', 'CAPE', 'CIN', 'RH_w', 'RH_ice', 'T_dew', 'T_frost', 'T_NAT', 'HNO3', 'OH', 'H2O2', 'press_bl', 'np', 'np_trop', 'np_CAPE']

print(header)
print()
print('reading file:   ', file)

data1 = pd.read_csv(file, delim_whitespace=True, names=header, comment='#', dtype=np.float64, na_values="NAN")

time = xr.DataArray(data1.time)
altitude = xr.DataArray(data1.altitude)
longitude = xr.DataArray(data1.longitude)
latitude =  xr.DataArray(data1.latitude)
H2O =  xr.DataArray(data1.H2O)
O3 = xr.DataArray(data1.O3)
temp = xr.DataArray(data1.temperature)

# ERA5 2° x 2°
nlon = 181
nlat = 91

lat_new = np.array(latitude).reshape(nlat,nlon)
lon_new = np.array(longitude).reshape(nlat,nlon)
o3_new = np.array(O3).reshape(nlat,nlon)*1.e9    #ppbv
h2o_new = np.array(H2O).reshape(nlat,nlon)*1.e6  # ppmv
temp_new =  np.array(temp).reshape(nlat,nlon)

print()
print('level = ',level, 'altitude = ',np.array(altitude[0]),'km')
print()

if level == '10_km':
  height = '10 km'
elif level == '5_km':
  height = '5 km'
elif level == '2_km':
  height = '2 km'

if (os.path.isdir('plots')==0):
  os.mkdir('plots')

plt.figure(1, figsize=(10,6))
plt.contourf(lon_new[0,:], lat_new[:,0], o3_new, cmap='viridis')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
cbar = plt.colorbar(label = 'O$_3$ (ppbv)')
plt.title('ECMWF ERA5 2° x 2° ('+height+')')
plt.savefig('plots/'+data+'_O3.png', bbox_inches='tight')

plt.figure(2, figsize=(10,6))
plt.contourf(lon_new[0,:], lat_new[:,0], h2o_new, cmap='viridis')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
cbar = plt.colorbar(label='H$_2$O (ppmv)')
plt.title('ECMWF ERA5 2° x 2° ('+height+')')
plt.savefig('plots/'+data+'_H2O.png', bbox_inches='tight')

plt.figure(3, figsize=(10,6))
plt.contourf(lon_new[0,:], lat_new[:,0], temp_new, cmap='viridis')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
cbar = plt.colorbar(label='Temperature (K)')
plt.title('ECMWF ERA5 2° x 2° ('+height+')')
plt.savefig('plots/'+data+'_temp.png', bbox_inches='tight')

plt.show()
