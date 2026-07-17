#! /bin/bash

set -e

# Create a local virtual environment on first use and reuse it afterwards.
if [ ! -d venv ]; then
  python3 -m venv venv
fi
source ./venv/bin/activate

# Install the plotting dependencies required by the example scripts.
pip install -q cartopy matplotlib numpy pandas scipy xarray

# Render all plots in batch mode without opening GUI windows.
export MPLBACKEND=Agg
datadir="../example/data.ref"
plotdir="plots"

echo "Generate parcel maps..."
python3 ./plot_atm.py "$datadir" "$plotdir"

echo "Generate 3D parcel maps..."
python3 ./plot_atm_3d.py "$datadir" "$plotdir"

echo "Generate column density maps..."
python3 ./plot_grid.py "$datadir" "$plotdir"

echo "Generate meteorological maps..."
python3 ./plot_met_map.py

echo "Generate trajectory map..."
python3 ./plot_traj.py
