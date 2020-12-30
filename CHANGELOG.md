# Changelog
All notable changes to this project will be documented in this file.

## [v2.0] - 2020-12-30
* Revised GPU code.
* Minor fix to allow for compilation on JURECA-DC and JUWELS.
* Update to gsl-2.6 and netcdf-c-4.7.4.
* Added HDF5 and zlib libraries to enable netCDF-4.
* Modified read_met for ERA5 netcdf files created by new cdo version.
* Added NVTX markers.
* Minor changes in Makefile.
* Publication of releases on zenodo.

## [v1.9] - 2020-03-11
* Added module to calculate wet deposition.
* Added option to read surface geopotential for NCEP reanalysis.
* Revised output of cloud parameters.
* Bugfix geopot function.
* Added new function to read surface data.
* Added cloud liquid water content and ice water content to meteo data.
* Added simple scheme for SO2 chemistry.
* Added OH climatology.
* Revised interpolation schemes for meteo data.
* Added code to calculate relative humidity.
* Revised code for multi-GPU runs.
* Replaced extract tool by atm_select.
* Modified netCDF output of tropo tool.
* Fixed spline interpolation in tropopause code.
* Modified interpolation functions to deal with missing values.

## [v1.8] - 2019-07-14
* Added new tool to sample tropopause data.
* Added tropopause water vapor to meteo tools.
* Added stride for particle output.
* Added water vapor to tropo tool.
* Revised and added test for decay module.
* Added GPU code.
* Revised example.
* Added fix for polar latitudes to PV code.
* Added OpenMP parallelization to tropopause code.
* Added tool to create tropopause climatology.
* Modified calculation of lapse rate.
* Faster calculation of geopotential heights and tropopause.
* Use virtual temperature to calculate geopotential heights.
* Revised smoothing filter in met_geopot.
* Added initial check of particle position.
* Added vertical interpolation to met_prof and met_zm.
* Added split based on kernel function to atm_split.
* Removed tools smago and wind.
* Added more flexible scheme for meteo filenames.
* Changed parameters of some tools.
* Added spatial interpolation to met_map.
* Added number of particles to grid output.
* Removed calculation of tracer conservation errors from atm_dist tool.
* Added filters to atm_stat tool.
* Removed ensemble filter from atm_dist tool.
* Update air parcel statistics tool.
* Renamed some of the tools for consistency.
* Added check that meteo grids are consistent.
* Added periodic boundary conditions in met_sample.
* Update PV code.
* Update GSL and netCDF libraries.
* New control parameter MET_DT_OUT.
* New tool "invert" to estimate volcanic emissions.
* Added day2doy and doy2day utilities.
* Added code to calculate dynamical tropopause.
* Revised time step control.
* Modifications of grid and profile output.
* Revised downsampling code.
* Enabled separate calculation of horizontal and vertical diffusion.

## [v1.7] - 2018-12-21
* Modified met_sample to work with geopotential heights.
* Added quantity flags to extraxt horizontal wid and vertical velocities.
* Allow for different grid spacing of geopotential data.
* Added code to calculate tropopause.
* Added low-pass filter to downsampling code.
* Changed standard pressure for calculation of potential temperature.
* Fixed initialization of randum number generators.
* Revised code to calculate PV on pressure levels.
* Modified approach to calculate RQTDs and RTCEs in dist tool.
* Added coordinate check when reading netCDF files.
* Added code to read CLaMS netCDF files.
* Added code for downsampling of meteo data.
* Added tool to convert between atmospheric data file formats.
* Added code to calculate tracer conservation errors in dist tool.
* Added relative transport deviations for quantities in dist tool.
* Tuning of dist and trac tool.
* Minor performance improvement of sedi module.
* Bugfix dist, met_zm and met_prof tools.
* Removed unused quantities.
* Extended Makefile to allow for compilation on KNL.
* Added option to strip binaries to Makefile.

## [v1.6] - 2018-07-29
* Added binary I/O for atmospheric data.
* Removed position checks from init tool.
* Fixed check of CFL criterion.
* Added code to stage meteorological data.
* Extended read_met function for ERA5 data.
* Added check for CFL criterion.
* Added low level wind components to meteo diagnostics.
* Added new tool for k-means clustering of trajectories.
* Added boundary conditions and options for evenly distribution to init tool.
* Added HNO3 climatology for PSC diagnostics.
* Added gravity wave diagnostics to meteo module.
* Modification of timers in the trac tool.
* Updates center, dist, meteo, trac and wind tool.

## [v1.5] - 2017-12-31
* Added module to gather PSC data.
* Added code to estimate T_ice and T_NAT.
* Added time filter for air parcel data.
* Added code to calculate potential vorticity.
* Output of surface pressure in profile data.

## [v1.4] - 2017-06-30
* Added code to read meteo data on model levels.
* Allow for meteo data at -180..+180 or 0..360 degrees.
* Implemented different options to get surface pressure.
* Removed nearest neighbour and cubic interpolation schemes.
* Added profile output.
* Revised grid and sample output.
* Use fixed units for known quantities.
* Removed sampling code.
* Added example script.

## [v1.3] - 2016-12-01
* Read time from meteo data filenames in read_met.
* Added different interpolation schemes for meteorological data.
* Added new tool to calculate zonal means of meteo data.
* Bugfix isosurface module.
* Update Smagorinsky code.
* Revised grid output.

## [v1.2] - 2016-06-27
* Added new tool to sample meteorological data files.
* Added new tool to calculate deviations of trajectories.
* Added code to force air parcels to balloon pressure levels.
* Added code to restrict trajectories to isosurfaces.
* Added code to create periodic boundary conditions for meteorological data.
* Added code to get zonal mean tropopause pressure from NCAR/NCEP climatology.
* Moidifed parameter names for the split tool.

## [v1.1] - 2015-12-23
* Added tool to create synthetic wind fields.
* Revisions based on Thomas's thesis.

## [v1.0] - 2015-09-21 
* This is an initial version of the MPTRAC model.
