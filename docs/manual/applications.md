# Applications and workflows

## Applications

The MPTRAC model comes with a number of individual programs or
applications. The most important app of MPTRAC is the tool
[trac](apps/trac.md), which is used to conduct the trajectory
calculations. The tables below give an overview of the available apps
and tools by application field.

All apps and tools accept `-h` and `--help` to print a short command-line
help message.

### Trajectory calculations

| App | Purpose |
| --- | --- |
| [trac](apps/trac.md) | Run forward or backward trajectory calculations. |

### Atmospheric particle data

| App | Purpose |
| --- | --- |
| [atm_init](apps/atm_init.md) | Create atmospheric particle files with initial trajectory seeds. |
| [atm_split](apps/atm_split.md) | Split particle sets into larger sets while retaining total mass. |
| [atm_select](apps/atm_select.md) | Extract subsets of atmospheric particle data, such as individual trajectories. |
| [atm_stat](apps/atm_stat.md) | Calculate statistics for atmospheric particle data. |
| [atm_dist](apps/atm_dist.md) | Calculate transport deviations between trajectory sets. |
| [atm_conv](apps/atm_conv.md) | Convert atmospheric particle data between supported file formats. |
| [atm2grid](apps/atm2grid.md) | Convert atmospheric particle data to gridded output. |

### Meteorological data

| App | Purpose |
| --- | --- |
| [met_conv](apps/met_conv.md) | Convert meteorological data between supported file formats. |
| [met_check_dt](apps/met_check_dt.md) | Check model time-step constraints for meteorological data. |
| [met_lapse](apps/met_lapse.md) | Calculate lapse-rate statistics from meteorological data. |
| [met_map](apps/met_map.md) | Extract maps from meteorological data. |
| [met_prof](apps/met_prof.md) | Extract vertical profiles from meteorological data. |
| [met_sample](apps/met_sample.md) | Sample meteorological data at atmospheric particle locations. |
| [met_spec](apps/met_spec.md) | Perform spectral analysis of meteorological temperature fields. |
| [met_subgrid](apps/met_subgrid.md) | Calculate subgrid-scale wind and vertical velocity statistics. |
| [met_zm](apps/met_zm.md) | Extract zonal means from meteorological data. |
| [cape](apps/cape.md) | Add CAPE data to meteorological netCDF files. |
| [sedi](apps/sedi.md) | Calculate sedimentation velocities for aerosol or cloud particles. |
| [tnat](apps/tnat.md) | Calculate polar stratospheric cloud formation temperatures. |
| [wind](apps/wind.md) | Create synthetic meteorological wind fields for test cases. |

### Tropopause data

| App | Purpose |
| --- | --- |
| [tropo](apps/tropo.md) | Create tropopause data sets from meteorological data. |
| [tropo_clim](apps/tropo_clim.md) | Calculate climatological statistics from tropopause data. |
| [tropo_sample](apps/tropo_sample.md) | Sample tropopause data at atmospheric particle locations. |
| [tropo_zm](apps/tropo_zm.md) | Extract zonal means from tropopause data. |

### Time conversion

| App | Purpose |
| --- | --- |
| [day2doy](apps/day2doy.md) | Convert a calendar date to day of year. |
| [doy2day](apps/doy2day.md) | Convert day of year to a calendar date. |
| [jsec2time](apps/jsec2time.md) | Convert MPTRAC seconds since 2000-01-01, 00:00 UTC to UTC time. |
| [time2jsec](apps/time2jsec.md) | Convert UTC time to MPTRAC seconds since 2000-01-01, 00:00 UTC. |

Please see the
[list of files in the Doxygen manual](https://slcs-jsc.github.io/mptrac/doxygen/files.html)
for the source-code reference of all apps and tools.

## Workflows

The individual apps can be connected to more comprehensive simulation
workflows by means of bash scripts.

This is an example showing how `atm_init` and `atm_split` are used to
initialize a simulation, and `trac` is used to calculate the trajectories:

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

Please see the page on [control parameters](control-parameters.md) for more information.
