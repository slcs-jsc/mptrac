# Python Script  Examples 

Here we provide some example Python scripts that can be used to plot MPTRAC output.

## Plotting of atm and grid output of MPTRAC

Three programs are provided plot_atm.py, plot_atm_3d.py and plot_grid.py.

Call the python plot scripts with parameter _datadir_ (the directory of the data files) and _plotdir_ (the directory where the plots will be written to), e.g for plot_atm.py it would be:

      python ./plot_atm.py ../example/data pythonplots

The libraryies required for running these python scripts are:

* scipy-stack
* Cartopy

Under Linux OS missing libraries can be downloaded and installed by using sudo apt-get install, e.g. for cartopy and scipy-stack it would be: 

      sudo apt-get install python3-scipy python3-cartopy

Additionally a shell script plot.sh is provided where all three programs can be run at once. The input and output files are then provided in the shell script and the script itself can be run with

       ./plot.sh

## Plotting of atm_select and met_map output of MPTRAC

To plot the data created with atm_select call the python script atm_traj.py using Python3

      python3 plot_traj.py

Here, a single trajectory on a map is plotted from the datafile traj_5450.tab provided in the directory data. 

To plot data created with met_map call the python script plot_met_map.py using Python3  

      python3 plot_met_map.py

Here, three parameters (O3, H2O and temperature) are plotted vs longitude and latitude from the datafile map_era5_2017010817_2_2.tab provided in the directory data. The met_map file contains data at 10 km altitude. In the program other altitudes can be chosen by changing the parameter "level". However, in order to plot other altitudes new files would be needed to be created.       

The python libraries required to run these scripts are:

* numpy
* pandas
* xarray
* matplotlib
* Basemap
