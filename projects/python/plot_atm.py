"""
@author: Mingzhao Liu
"""

# Import modules...
import matplotlib.pyplot as plt
import matplotlib.colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import pandas as pd
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import glob,sys,os

# Check arguments...
if (len(sys.argv)-1 != 2):
  print('Please input <datadir> <plotdir>')
  exit()

# Setup...
datadir = sys.argv[1]
plotdir = sys.argv[2]
if (os.path.isdir(plotdir)==0):
  os.mkdir(plotdir)

# Set plot range...
lon0 = -90
lon1 = 60
lat0 = -65
lat1 = -15

# Set position of volcano...
vlon = -72.117
vlat = -40.59

# Set colorbar...
mycmap = matplotlib.colors.LinearSegmentedColormap.from_list( \
        'mycmap', [(0 , 'whitesmoke'), (1 / 6  , 'blue'), \
        (2 / 6  , 'cyan'),(3 / 6  , 'green'),(4 / 6  , 'yellow'), \
        (5 / 6  , 'orange'),(1, 'red')])

# Loop over files...
filelist = sorted(glob.glob(f'{datadir}/atm_20*.tab'))
for FILE in filelist :

  # Write info...
  print('Plot '+FILE+' ...')
  
  # Set time...
  minute = FILE.split('_')[-1][:2]
  hour = FILE.split('_')[-2]
  day = FILE.split('_')[-3]
  month = FILE.split('_')[-4]
  year = FILE.split('_')[-5]

  # Read data...
  data = pd.read_csv(FILE, comment='#', sep=' ', header=None).values
  lon = data[:,2]
  lat = data[:,3]
  alt = data[:,1]

  # Plot...
  fig = plt.figure(dpi=200)
  ax = fig.add_subplot(projection=ccrs.PlateCarree(central_longitude=0))
  ax.set_aspect(2)
  ax.set_title(f'MPTRAC | {year}-{month}-{day} {hour}:{minute} UTC', fontsize=12)
  ax.add_feature(cfeature.COASTLINE, lw=0.5)
  ax.set_xticks(np.arange(-180,180,30), crs=ccrs.PlateCarree(central_longitude=0))
  ax.set_yticks(np.arange(-90,90,10), crs=ccrs.PlateCarree())
  ax.xaxis.set_major_formatter(LongitudeFormatter())
  ax.yaxis.set_major_formatter(LatitudeFormatter())
  ax.set_extent([lon0, lon1, lat0, lat1], crs=ccrs.PlateCarree(central_longitude=0))
  c=ax.scatter(lon,lat,c=alt,s=1,marker=',', vmin=5, vmax=15, linewidth=0,
               cmap = mycmap,transform=ccrs.PlateCarree(central_longitude=0))
  ax.scatter(vlon,vlat, c='r',marker='^', \
             transform=ccrs.PlateCarree(central_longitude=0))
  ax.tick_params(labelsize=12)
  fig.colorbar(c,ax=ax,shrink=0.75,pad=0.01).set_label('Altitude [km]',fontsize=12)
  plt.savefig(f'{plotdir}/atm_{year}_{month}_{day}_{hour}_{minute}.png', bbox_inches='tight')
  plt.close(fig)
