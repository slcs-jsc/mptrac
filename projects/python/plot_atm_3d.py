"""Plot MPTRAC air parcel snapshots in a stylized 3-D view."""

# Import modules...
import glob,sys,os
from pathlib import Path

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

# Use a warm-to-cool palette to emphasize parcel mass in the plume.
mycmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'parcel_mass', ['#d8ecf8', '#67b7dc', '#1d6996', '#2ca25f', '#ffd166', '#f8961e', '#c1121f'])

# Define the regional domain and altitude range shown in the 3-D view.
lon0 = -90
lon1 = 60
lat0 = -65
lat1 = -15
alt0 = 5
alt1 = 15

# Mark the volcano location for geographic reference.
vlon = -72.117
vlat = -40.59

# Process all parcel snapshot files in the selected input directory.
filelist = sorted(glob.glob(f'{datadir}/atm_20*.tab'))
for FILE in filelist:

  print(f"  [atm3d] {Path(FILE).name}")

  # Extract the timestamp from the standard MPTRAC filename.
  second = FILE.split('_')[-1].split('.')[0]
  minute = FILE.split('_')[-2]
  hour = FILE.split('_')[-3]
  day = FILE.split('_')[-4]
  month = FILE.split('_')[-5]
  year = FILE.split('_')[-6]

  # Read longitude, latitude, altitude, and mass from the parcel table.
  data = pd.read_csv(FILE, comment='#', sep=r'\s+', header=None).values
  lon = data[:,2]
  lat = data[:,3]
  alt = data[:,1]
  mass = data[:,4]

  # Render one 3-D scene per input file.
  fig = plt.figure(figsize=(9.6, 5.8), dpi=220)
  ax = fig.add_subplot(projection='3d')
  ax.set_position([0.06, 0.16, 0.88, 0.76])
  ax.set_title(f'MPTRAC air parcels | {year}-{month}-{day}, {hour}:{minute} UTC', fontsize=13, pad=4)
  ax.set_xlim(lon0,lon1)
  ax.set_ylim(lat0,lat1)
  ax.set_zlim(alt0,alt1)
  ax.set_xlabel('Longitude', labelpad=10)
  ax.set_ylabel('Latitude', labelpad=10)
  ax.set_zlabel('Altitude [km]', labelpad=8)
  ax.xaxis.set_major_formatter(LongitudeFormatter())
  ax.yaxis.set_major_formatter(LatitudeFormatter())
  ax.set_xticks(np.arange(lon0, lon1 + 1, 30))
  ax.set_yticks(np.arange(lat0, lat1 + 1, 10))
  ax.set_zticks(np.arange(alt0, alt1 + 1, 2))
  ax.view_init(elev=24, azim=-128)
  ax.set_box_aspect((2.35, 1.3, 0.82))

  # Use matching pane colors so the bounding box stays visually unobtrusive.
  pane_color = (0.92, 0.94, 0.96, 0.5)
  ax.xaxis.pane.set_facecolor(pane_color)
  ax.yaxis.pane.set_facecolor(pane_color)
  ax.zaxis.pane.set_facecolor(pane_color)
  ax.xaxis.pane.set_edgecolor((0.7, 0.7, 0.7, 0.45))
  ax.yaxis.pane.set_edgecolor((0.7, 0.7, 0.7, 0.45))
  ax.zaxis.pane.set_edgecolor((0.7, 0.7, 0.7, 0.15))
  ax.grid(True)
  for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
    axis._axinfo['grid']['linewidth'] = 0.45
    axis._axinfo['grid']['linestyle'] = '--'
    axis._axinfo['grid']['color'] = (0.45, 0.45, 0.45, 0.35)

  # Add a light ground projection to improve depth perception.
  ax.scatter(
      lon,
      lat,
      np.full_like(alt, alt0),
      c='#7f8c8d',
      s=2,
      marker='o',
      linewidth=0,
      alpha=0.06,
      depthshade=False,
      zorder=1,
  )

  c = ax.scatter(
      lon,
      lat,
      alt,
      c=mass,
      s=6,
      marker='o',
      linewidth=0,
      alpha=0.82,
      depthshade=False,
      vmin=200,
      vmax=250,
      cmap=mycmap,
      zorder=3,
  )

  ax.scatter([vlon], [vlat], [alt0], c='#b11226', marker='^', s=70, edgecolors='white', linewidth=0.7, zorder=5)
  ax.tick_params(labelsize=10, pad=2)

  cbar = fig.colorbar(c, ax=ax, orientation='horizontal', shrink=0.76, pad=0.018, fraction=0.04)
  cbar.set_label('Mass [kg]', fontsize=11)
  plt.savefig(f'{plotdir}/atm3d_{year}_{month}_{day}_{hour}_{minute}_{second}.png', bbox_inches='tight')
  plt.close(fig)
