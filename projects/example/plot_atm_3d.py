"""
@author: Mingzhao Liu
"""

import matplotlib.pyplot as plt
import matplotlib.colors
import cartopy.crs as ccrs
import numpy as np
import pandas as pd
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import glob,sys,os, cv2

datadir = 'data'
plotdir = 'pythonplots'
if (os.path.isdir(plotdir)==0):
  os.mkdir(plotdir)

mycmap = matplotlib.colors.LinearSegmentedColormap.from_list( \
        'mycmap', [(0 , 'whitesmoke'), (1 / 6  , 'blue'), \
        (2 / 6  , 'cyan'),(3 / 6  , 'green'),(4 / 6  , 'yellow'), \
        (5 / 6  , 'orange'),(1, 'red')])
lon0 = -90
lon1 = 60
lat0 = -65
lat1 = -15
alt0 = 5
alt1 = 15
vlon = -72.117 
vlat = -40.59
filelist = sorted(glob.glob(f'{datadir}/atm_20*.tab'))
for FILE in filelist :
  minute = FILE.split('_')[-1][:2]
  hour = FILE.split('_')[-2]
  day = FILE.split('_')[-3]
  month = FILE.split('_')[-4]
  year = FILE.split('_')[-5]
  data = pd.read_csv(FILE, comment='#', sep=' ', header=None).values
  lon = data[:,2]
  lat = data[:,3]
  alt = data[:,1]
  mass = data[:,4]

  fig = plt.figure(dpi=200)
  ax = fig.add_subplot(projection='3d')
  ax.set_zlabel('Altitude [km]')
  ax.set_xlim(lon0,lon1)
  ax.set_ylim(lat0,lat1)
  ax.set_zlim(alt0,alt1)
  ax.set_title(f'MPTRAC | {year}-{month}-{day} {hour}:{minute} UTC',fontsize=10 )
  ax.xaxis.set_major_formatter(LongitudeFormatter())
  ax.yaxis.set_major_formatter(LatitudeFormatter())
  ax.set_xticks(np.arange(lon0,lon1,30))
  ax.set_yticks(np.arange(lat0,lat1,10))
  c=ax.scatter(lon,lat,alt,c=mass,s=1,marker=',', 
       linewidth=0, vmin=200, vmax=250,
       cmap = mycmap)
  ax.plot([vlon,vlon],[vlat,vlat],[alt0,alt1], c='r',alpha=1)
  ax.tick_params(labelsize=10)
  fig.colorbar(c,ax=ax,shrink=0.75,pad=0.12).set_label('Mass [kg]',fontsize=10  )
  plt.savefig(f'{plotdir}/atm3d_{year}_{month}_{day}_{hour}_{minute}.png', bbox_inches='tight')
  plt.close(fig)
