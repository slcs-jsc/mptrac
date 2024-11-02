#! /bin/bash

# Activate venv...
source ~/venv/bin/activate

# Setup...
datadir="../example/data.ref"
plotdir="plots"

# Plot...
python3 ./plot_atm.py $datadir $plotdir
python3 ./plot_atm_3d.py $datadir $plotdir
python3 ./plot_grid.py $datadir $plotdir
python3 ./plot_met_map.py
python3 ./plot_traj.py
