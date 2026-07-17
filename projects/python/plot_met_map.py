######################################################################
## Author: Farahnaz Khosrawi
## Data created: 14.09.2023
## Last modified: 17.07.2026
## Purpose: Plot MPTRAC met_map output on geographic maps
######################################################################

import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# Define the MPTRAC met_map column layout used by the example input file.
header = [
    "time", "altitude", "longitude", "latitude", "pressure", "temperature",
    "u", "v", "w", "H2O", "O3", "gph", "pv", "psurf", "tsurf", "gph_surf",
    "u_surf", "v_surf", "land_sea_mask", "SST", "p_trop", "gph_trop", "temp_trop",
    "H2O_trop", "cloud_lw", "cloud_iwc", "cloud_cover", "tot_col_water",
    "cloud_top press", "cloud_bottom_press", "pressure_lifted_CL",
    "press_lev_free_conv", "press_EL", "CAPE", "CIN", "RH_w", "RH_ice",
    "T_dew", "T_frost", "T_NAT", "HNO3", "OH", "H2O2", "press_bl",
    "np", "np_trop", "np_CAPE",
]

# Use the shipped ERA5 example file and write all output to plots/.
data_name = "ERA5"
level = "10_km"
file = "data/map_era5_2017010817_2_2.tab"

print(f"  [met_map] {file}")

# Read the gridded meteorological fields written by MPTRAC met_map.
data = pd.read_csv(
    file,
    sep=r"\s+",
    names=header,
    comment="#",
    dtype=np.float64,
    na_values="NAN",
)

nlon = 181
nlat = 91
lon = data["longitude"].to_numpy().reshape(nlat, nlon)
lat = data["latitude"].to_numpy().reshape(nlat, nlon)
o3 = data["O3"].to_numpy().reshape(nlat, nlon) * 1.0e9
h2o = data["H2O"].to_numpy().reshape(nlat, nlon) * 1.0e6
temp = data["temperature"].to_numpy().reshape(nlat, nlon)
altitude = data["altitude"].to_numpy()[0]

print(f"  [met_map] altitude={altitude:.1f} km, level={level}")

height = {"10_km": "10 km", "5_km": "5 km", "2_km": "2 km"}.get(level, level)

if not os.path.isdir("plots"):
    os.mkdir("plots")


def plot_field(values, label, filename, cmap):
    # Render one global map for the selected variable.
    fig = plt.figure(figsize=(11, 5.5), dpi=200)
    ax = fig.add_subplot(projection=ccrs.PlateCarree())
    ax.set_title(f"ECMWF ERA5 2° x 2° | {height}", fontsize=13)
    ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.OCEAN.with_scale("110m"), facecolor="#dceffd", zorder=0)
    ax.add_feature(cfeature.LAND.with_scale("110m"), facecolor="#f4f1e8", zorder=0)
    ax.add_feature(cfeature.COASTLINE.with_scale("110m"), linewidth=0.6, edgecolor="#3a3a3a")
    ax.set_xticks(np.arange(-180, 181, 60), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(-90, 91, 30), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.gridlines(linestyle="--", linewidth=0.35, color="#7a7a7a", alpha=0.5)
    mesh = ax.contourf(lon[0, :], lat[:, 0], values, levels=18, cmap=cmap, transform=ccrs.PlateCarree())
    cbar = fig.colorbar(mesh, ax=ax, shrink=0.82, pad=0.03)
    cbar.set_label(label)
    plt.savefig(filename, bbox_inches="tight")
    plt.close(fig)


# Create one figure each for ozone, water vapour, and temperature.
plot_field(o3, "O$_3$ (ppbv)", f"plots/{data_name}_O3.png", "viridis")
plot_field(h2o, "H$_2$O (ppmv)", f"plots/{data_name}_H2O.png", "YlGnBu")
plot_field(temp, "Temperature (K)", f"plots/{data_name}_temp.png", "magma")
