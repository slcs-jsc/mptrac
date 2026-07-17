"""Plot MPTRAC air parcel snapshots on a regional map."""

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

# Focus the regional map on the 2011 Puyehue-Cordon Caulle plume domain.
lon0 = -90
lon1 = 60
lat0 = -65
lat1 = -15

# Mark the volcano location for geographic reference.
vlon = -72.117
vlat = -40.59

# Use a warm-to-cool palette to emphasize parcel altitude.
mycmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'parcel_altitude', ['#d8ecf8', '#65b6e8', '#2b8cbe', '#2ca25f', '#ffd166', '#f8961e', '#d62828'])

# Process all parcel snapshot files in the selected input directory.
filelist = sorted(glob.glob(f'{datadir}/atm_20*.tab'))
for FILE in filelist:

  print(f"  [atm] {Path(FILE).name}")

  # Extract the timestamp from the standard MPTRAC filename.
  second = FILE.split('_')[-1].split('.')[0]
  minute = FILE.split('_')[-2]
  hour = FILE.split('_')[-3]
  day = FILE.split('_')[-4]
  month = FILE.split('_')[-5]
  year = FILE.split('_')[-6]

  # Read longitude, latitude, and altitude from the parcel table.
  data = pd.read_csv(FILE, comment='#', sep=r'\s+', header=None).values
  lon = data[:,2]
  lat = data[:,3]
  alt = data[:,1]

  # Render one regional map per input file.
  fig = plt.figure(figsize=(11, 5.8), dpi=220)
  ax = fig.add_subplot(projection=ccrs.PlateCarree())
  ax.set_title(f'MPTRAC air parcels | {year}-{month}-{day}, {hour}:{minute} UTC', fontsize=13, pad=12)
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
  c = ax.scatter(
      lon,
      lat,
      c=alt,
      s=4,
      marker='o',
      vmin=5,
      vmax=15,
      linewidth=0,
      alpha=0.78,
      cmap=mycmap,
      transform=ccrs.PlateCarree(),
      zorder=2,
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
  cbar.set_label('Altitude [km]', fontsize=11)
  plt.savefig(f'{plotdir}/atm_{year}_{month}_{day}_{hour}_{minute}_{second}.png', bbox_inches='tight')
  plt.close(fig)
