# Applications and workflows

## Applications

The MPTRAC model comes with a number of individual programs or applications. The most important app of MPTRAC is the tool [trac](https://github.com/slcs-jsc/mptrac/blob/documentation/docs/manual/apps/trac.md), which is used to conduct the trajectory calculations.

The apps met_map, met_prof, and met_zm can be used to extract global maps, vertical profiles, and zonal means from meteorological data. The app met_sample can be used to sample the meteo data at individual locations in space and time.

The app [atm_conv](https://github.com/slcs-jsc/mptrac/blob/documentation/docs/manual/apps/atm_conv.md) can be used to convert between different file formats of the particle data (ASCII, binary, netCDF). The app atm_dist can be used to calculate transport deviations between trajectory sets. The app atm_init can be used to create particle data files with initial trajectory seeds. The app atm_select can extract subsets of the particle data, like individual trajectories. The app atm_split can split sets of particles into larger sets, retaining their total mass. The app atm_stat calculates trajectory statistics, for example, the mean position.

The tools day2doy, doy2day, jsec2time, and time2jsec are used for time conversion. In particular, they can be used to determine the day of the year (doy) for a given date and convert between a UTC time (YYYY-MM-DD, HH:MM:SS) and the absolute time in seconds since 2000-01-01, 00:00 UTC (the internal time coordinate of MPTRAC).

The tools lapse, tropo, and tropo_sample can be used to determine lapse rate statistics and to prepare and sample tropopause data files.

Please see the [list of files in the doxygen manual](https://slcs-jsc.github.io/mptrac/doxygen/files.html) for more information.

## Workflows

The individual apps can be connected to more comprehensive simulation workflows by means of bash scripts.

This is an example showing how atm_init and atm_split are used to initialize a simulation, trac is used to calculate the trajectories, and atm_dist is used to calculate transport deviations:

```
#! /bin/bash

# Setup...
trac=../src

# Create directories...
mkdir -p data plots

# Set time range of simulations...
t0=$($trac/time2jsec 2011 6 5 0 0 0 0)
t1=$($trac/time2jsec 2011 6 8 0 0 0 0)

# Set initial air parcel positions...
$trac/atm_init trac.ctl data/atm_init.tab \
               INIT_T0 $t0 INIT_T1 $t0 \
               INIT_Z0 10.0 INIT_Z1 10.0 \
               INIT_LON0 -72.117 INIT_LON1 -72.117 \
               INIT_LAT0 -40.59 INIT_LAT1 -40.59

# Split air parcels...
$trac/atm_split trac.ctl data/atm_init.tab data/atm_split.tab \
                SPLIT_N 10000 SPLIT_M 1e9 SPLIT_DX 30.0 SPLIT_DZ 1.0

# Calculate trajectories...
echo "data" > data/dirlist
$trac/trac data/dirlist trac.ctl atm_split.tab \
           ATM_BASENAME atm GRID_BASENAME grid
```

Please see the wiki page on [Control parameters](https://github.com/slcs-jsc/mptrac/blob/documentation/docs/manual/control-parameters.md) for more information on the control file trac.ctl and the control parameters.
