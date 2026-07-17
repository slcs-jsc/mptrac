"""Plot MPTRAC column density output on a regional map."""

# Import modules...
import glob,sys,os
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

# Expect input and output directories on the command line.
if (len(sys.argv)-1 != 2):
  print('Please input <datadir> <plotdir>')
  exit()

# Read the source data directory and create the output directory if needed.
datadir = sys.argv[1]
plotdir = sys.argv[2]
if (os.path.isdir(plotdir)==0):
  os.mkdir(plotdir)

# Define the regional plotting domain and the regular grid dimensions.
lon0 = -90
lon1 = 60
lat0 = -65
lat1 = -15
nlon = 300
nlat = 100

# Mark the volcano location for geographic reference.
vlon = -72.117
vlat = -40.59

# Use a nonlinear normalization to keep weak and strong column signals visible.
mycmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'column_density', ['#f7fbff', '#d0e1f2', '#73b3d8', '#2879b9', '#1fa187', '#8bd646', '#fde725', '#f8961e'])
norm = matplotlib.colors.PowerNorm(gamma=0.65, vmin=0, vmax=1)

# Process all gridded output files in the selected input directory.
filelist = sorted(glob.glob(f'{datadir}/grid_20*.tab'))
for FILE in filelist:

  print(f"  [grid] {Path(FILE).name}")

  # Extract the timestamp from the standard MPTRAC filename.
  second = FILE.split('_')[-1].split('.')[0]
  minute = FILE.split('_')[-2]
  hour = FILE.split('_')[-3]
  day = FILE.split('_')[-4]
  month = FILE.split('_')[-5]
  year = FILE.split('_')[-6]

  # Read the gridded column density field and reshape it to a 2-D map.
  data = pd.read_csv(FILE, comment='#', sep=r'\s+', header=None).values
  lat = data[:,3][0:nlat]
  lon = np.zeros((nlon))
  for i in range(nlon):
    lon[i] = data[:,2][nlat*i]
  col = (data[:,6]*1000).reshape(nlat,nlon,order='F')

  # Render one regional map per input file.
  fig = plt.figure(figsize=(11, 5.8), dpi=220)
  ax = fig.add_subplot(projection=ccrs.PlateCarree())
  ax.set_title(f'MPTRAC column density | {year}-{month}-{day}, {hour}:{minute} UTC', fontsize=13, pad=12)
  ax.set_extent([lon0, lon1, lat0, lat1], crs=ccrs.PlateCarree())
  ax.add_feature(cfeature.OCEAN.with_scale('50m'), facecolor='#dceffd', zorder=0)
  ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor='#f4f1e8', zorder=0)
  ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.35, edgecolor='#707070')
  ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.6, edgecolor='#3a3a3a')
  ax.set_xticks(np.arange(-90, 61, 30), crs=ccrs.PlateCarree())
  ax.set_yticks(np.arange(-60, -14, 10), crs=ccrs.PlateCarree())
  ax.xaxis.set_major_formatter(LongitudeFormatter())
  ax.yaxis.set_major_formatter(LatitudeFormatter())
  ax.gridlines(linestyle='--', linewidth=0.35, color='#7a7a7a', alpha=0.45)
  c = ax.pcolormesh(
      lon,
      lat,
      col,
      cmap=mycmap,
      norm=norm,
      shading='auto',
      transform=ccrs.PlateCarree(),
      zorder=1,
  )
  ax.scatter(
      vlon,
      vlat,
      c='#b11226',
      marker='^',
      s=90,
      edgecolors='white',
      linewidth=0.7,
      transform=ccrs.PlateCarree(),
      zorder=3,
  )
  ax.tick_params(labelsize=11)
  cbar = fig.colorbar(c, ax=ax, orientation='horizontal', shrink=0.9, pad=0.08, fraction=0.06)
  cbar.set_label('Column density [g/m$^2$]', fontsize=11)
  plt.savefig(f'{plotdir}/grid_{year}_{month}_{day}_{hour}_{minute}_{second}.png', bbox_inches='tight')
  plt.close(fig)
