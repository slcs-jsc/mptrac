# Python Script  Examples 

Here we provide some example Python scripts that can be used to plot MPTRAC output.

## Plotting of atm and grid output of MPTRAC

Call the python plot scripts with parameter <datadir> and <plotdir>

      python ./plot_atm.py ../example/data pythonplots

Library required

* scipy-stack
* Cartopy

      sudo apt-get install python3-scipy python3-cartopy

## Plotting of atm_select and met_map output of MPTRAC

To plot the data created with atm_select call the python script atm_traj.py using Python3

      python3 plot_traj.py

Here, a single trajectory on a map is plotted from the datafile traj_5450.tab provided in the directory data. 

To plot data created with met_map call the python script plot_met_map.py using Python3  

      python3 plot_met_map.py

Here, three parameters (O3, H2O and temperature) are plotted vs longitude and latitude from the datafile map_era5_2017010817_2_2.tab provided in the directory data.       

The python libraries required to run these scripts are:

* numpy
* pandas
* xarray
* matplotlib
* Basemap
