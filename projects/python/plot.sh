#! /bin/bash

datadir="../example/data.ref"
plotdir="plots"

python ./plot_atm.py $datadir $plotdir
python ./plot_atm_3d.py $datadir $plotdir
python ./plot_grid.py $datadir $plotdir
