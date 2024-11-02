#! /bin/bash

datadir="../example/data.ref"
plotdir="plots"

python3 ./plot_atm.py $datadir $plotdir
python3 ./plot_atm_3d.py $datadir $plotdir
python3 ./plot_grid.py $datadir $plotdir
