"""Plot a single MPTRAC trajectory on a global Cartopy map."""

import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# Select the example trajectory file that will be rendered.
traj_no = "5450"
file = f"data/traj_{traj_no}.tab"
header = ["time", "altitude", "longitude", "latitude"]

print(f"  [traj] {file}")

# Read the trajectory table written by MPTRAC atm_select.
data = pd.read_csv(
    file,
    sep=r"\s+",
    names=header,
    comment="#",
    dtype=np.float64,
    na_values="NAN",
)

altitude = data["altitude"].to_numpy()
longitude = data["longitude"].to_numpy()
latitude = data["latitude"].to_numpy()

if not os.path.isdir("plots"):
    os.mkdir("plots")

fig = plt.figure(figsize=(11, 5.5), dpi=200)
ax = fig.add_subplot(projection=ccrs.PlateCarree())
ax.set_title(f"Trajectory {traj_no} | Start altitude {altitude[0]:.2f} km", fontsize=13)
ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.OCEAN.with_scale("110m"), facecolor="#dceffd", zorder=0)
ax.add_feature(cfeature.LAND.with_scale("110m"), facecolor="#f4f1e8", zorder=0)
ax.add_feature(cfeature.COASTLINE.with_scale("110m"), linewidth=0.6, edgecolor="#3a3a3a")
ax.add_feature(cfeature.BORDERS.with_scale("110m"), linewidth=0.3, edgecolor="#666666")
ax.set_xticks(np.arange(-180, 181, 60), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(-90, 91, 30), crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
ax.gridlines(linestyle="--", linewidth=0.35, color="#7a7a7a", alpha=0.5)

# Split the line at dateline crossings to avoid artificial wrap-around segments.
jump_idx = np.where(np.abs(np.diff(longitude)) > 180)[0] + 1
segments = np.split(np.column_stack((longitude, latitude)), jump_idx)
for segment in segments:
    ax.plot(
        segment[:, 0],
        segment[:, 1],
        color="#1f78b4",
        linewidth=1.6,
        transform=ccrs.PlateCarree(),
        zorder=2,
    )

# Mark the start and end points so the trajectory direction remains visible.
ax.scatter(
    longitude[0],
    latitude[0],
    marker="^",
    s=48,
    color="#d7301f",
    edgecolor="white",
    linewidth=0.5,
    transform=ccrs.PlateCarree(),
    zorder=3,
)
ax.scatter(
    longitude[-1],
    latitude[-1],
    marker="o",
    s=24,
    color="#253494",
    edgecolor="white",
    linewidth=0.4,
    transform=ccrs.PlateCarree(),
    zorder=3,
)

plt.savefig("plots/traj.png", bbox_inches="tight")
plt.close(fig)
