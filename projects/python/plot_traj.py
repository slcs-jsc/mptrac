######################################################################
## Author: Farahnaz Khosrawi
## Data created: 11.10.2023
## Last modified: 11.12.2023
## Purpose: Plots single air parcel trajectories on a map - MPTRAC
## atm_select output
######################################################################

from mpl_toolkits.basemap import Basemap
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import sys,os

traj_no = '5450'

file = 'data/traj_'+traj_no+'.tab'

header = ['time', 'altitude', 'longitude', 'latitude']
print(header)
print('reading file:   ', file)

# read trajectory data...
data1 = pd.read_csv(file, delim_whitespace=True, names=header, comment='#', dtype=np.float64, na_values="NAN")

time = xr.DataArray(data1.time)
altitude = xr.DataArray(data1.altitude)
longitude = xr.DataArray(data1.longitude)
latitude = xr.DataArray(data1.latitude)

# set map projection...
m = Basemap(llcrnrlon=-180.,llcrnrlat=-90.,urcrnrlon=180.,urcrnrlat=90.,\
            resolution='l',projection='cyl',\
            lat_0=0.,lon_0=180.)

m.drawcoastlines()
m.fillcontinents()
# draw parallels
m.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
# draw meridians
m.drawmeridians(np.arange(-180,180,60),labels=[0,0,0,1])

# plot trajectory...
plt.plot(longitude[:-1], latitude[:-1], scaley=True, scalex=True, color = 'b')

# plot trajectory start point...
plt.plot(longitude[0], latitude[0], marker = 'x', color = 'k')

start_alt = np.round(np.array(altitude[0]), 2)
plot_title = ('Traj. No. '+str(traj_no)+', '+str(start_alt)+' km')
plt.title(plot_title)

if (os.path.isdir('plots')==0):
  os.mkdir('plots')
plt.savefig('plots/traj.png', format = 'png')

plt.show()
