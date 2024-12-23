#! /bin/bash

# Create Python environment...
python3 -m venv venv
source ./venv/bin/activate
pip install basemap basemap-data cartopy matplotlib numpy pandas scipy xarray

# Setup...
datadir="../example/data.ref"
plotdir="plots"

# Plot...
python3 ./plot_atm.py $datadir $plotdir
python3 ./plot_atm_3d.py $datadir $plotdir
python3 ./plot_grid.py $datadir $plotdir
python3 ./plot_met_map.py
python3 ./plot_traj.py
