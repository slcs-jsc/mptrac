# Meteorological input data

## Introduction

The MPTRAC model needs meteorological input data (winds, temperature, ...) to perform the trajectory calculations and for other tasks.

 The meteo data should be provided as netCDF files on a regular grid (pressure, longitude, latitude). Each time step of the reanalysis data should be provided as a separate netCDF file. The time step of the corresponding meteo data file is inferred from the filename following this naming convention:

```
     <basename>_<YYYY>_<MM>_<DD>_<HH>.nc
```

At a minimum, MPTRAC needs temperature (T), zonal wind (U), meridional wind (V), and vertical velocity (W). Additionally, many applications need specific humidity (Q), ozone (O3), cloud liquid water content (CLWC), and cloud ice water content (CIWC).

It is also necessary to download the surface pressure (PS or LNSP) and the surface geopotential height (Z). Optionally, surface temperature (T2M) and winds (U10M, V10M) are read in.

It is possible to provide meteo data on model levels rather than pressure levels. In this case, pressure (PL) on model levels is required as an additional input variable in order to perform a vertical interpolation from model levels to pressure levels.

See the example/meteo/ directory in the MPTRAC repository for an example of the meteo data netCDF files.

## The ECMWF reanalyses

### Overview

As an example, we discuss the use of [ERA-Interim](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era-interim) and [ERA5](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5) reanalysis data provided by the European Centre for Medium-Range Weather Forecasts (ECMWF) to run the MPTRAC model.

ERA-Interim is a global atmospheric reanalysis that is available from 1 January 1979 to 31 August 2019. The data assimilation system used to produce ERA-Interim is based on a 2006 release of the IFS (Cy31r2). The system includes a 4-dimensional variational analysis (4D-Var) with a 12-hour analysis window. The spatial resolution of the data set is approximately 80 km (T255 spectral) on 60 levels in the vertical from the surface up to 0.1 hPa.

A first segment of the ERA5 dataset is now available for public use (1979 to within 5 days of real time). ERA5 provides hourly estimates of a large number of atmospheric, land and oceanic climate variables. The data cover the Earth on a 30 km grid and resolve the atmosphere using 137 levels from the surface up to a height of 80 km. ERA5 includes information about uncertainties for all variables at reduced spatial and temporal resolutions.

### Download of ERA-Interim data

The first step is to download the ERA-Interim data from [ECMWF's data archive](https://apps.ecmwf.int/datasets/data/interim-full-daily). It is recommended to download the data in grib file format at 0.75°x0.75° horizontal resolution, on all 60 model levels, and at 6-hourly time intervals.

*Alternatively, download ERA5 from the [Copernicus Climate Change Service (C3S) Climate Date Store](https://cds.climate.copernicus.eu/#!/search?text=ERA5&type=dataset). It is recommended to download the data in grib file format at 0.3°x0.3° horizontal resolution, on all 137 model levels, and at hourly time intervals.

### Interpolation to pressure levels

The second step is to interpolate the ERA-Interim data from model levels to pressure levels. Here, a set of target pressure levels needs to be selected. It is recommended to use a set of 60 pressure levels that matches the vertical sampling of the model level data (as defined by a fixed surface pressure of 1013.25 hPa and [ECMWF's a and b level coefficients](https://www.ecmwf.int/en/forecasts/documentation-and-support/60-model-levels)).

For ERA5 a different set of [a and be level coeffients](https://www.ecmwf.int/en/forecasts/documentation-and-support/137-model-levels) should be considered. Potentially, the vertical sampling of the target pressure levels in the lower troposphere and the number of levels in the mesosphere can be reduced to reduce the amount of data.

The CDO tools can be used to perform the vertical interpolation from model levels to pressure levels and to convert from grib to netCDF file format. For example:

```
    # merge surface and model level grib files...
    cdo -P 8 merge 2020010100_sf.grb 2020010100_ml.grb merge.grib

    # interpolate to pressure levels and create netCDF file (for ERA-Interim)...
    cdo -P 8 -a -f nc -t ecmwf ml2pl,101325,101085,100726,100202,99474.7,98514,97298.7,95814.8,94055.1,92018.9,89711.2,87142,84326.3,81283,78034.6,74606.3,71026.3,67324,63530.6,59677.7,55797.3,51920.9,48079.1,44300.9,40613.3,37040.7,33604.4,30321.7,27205.9,24265.2,21502.5,18914.7,16508.9,14290.2,12261.4,10422.9,8772.74,7306.63,6018.02,4906.71,3960.29,3196.42,2579.89,2082.27,1680.64,1356.47,1094.83,883.66,713.22,575.65,464.62,373.97,298.5,234.78,180.58,134.48,95.64,63.65,38.43,20,10 merge.grib ei_2020_01_01_00.nc

    # interpolate to pressure levels and create netCDF file (for ERA5)...
    cdo -P 8 -a -f nc -t ecmwf ml2pl,101325,98660.4,93476.7,84269.6,81208.5,79564,77846.6,76060,74208.6,72297.9,70334.7,68326.2,66280.8,64207.6,62116.2,60016.7,57919.3,55834.3,53772,51742,49758.4,47831,45963.2,44153.9,42401.9,40705.8,39064.5,37476.7,35941.1,34456.6,33022,31636.1,30297.6,29005.5,27758.5,26555.6,25395.5,24277.2,23199.5,22161.5,21161.9,20199.7,19273.9,18383.4,17527.3,16704.5,15914,15154.9,14426.2,13727,13056.4,12413.4,11797.1,11206.8,10641.5,10100.5,9582.8,9087.74,8614.5,8161.82,7728.1,7311.87,6911.87,6526.95,6156.07,5798.34,5452.99,5119.9,4799.15,4490.82,4194.93,3911.49,3640.47,3381.74,3135.12,2900.39,2677.35,2465.77,2265.43,2076.1,1897.52,1729.45,1571.62,1423.77,1285.61,1156.85,1037.2,926.34,823.97,729.74,643.34,564.41,492.62,427.59,368.98,316.42,269.54,227.97,191.34,159.28,131.43,107.42,86.9,69.52,54.96,42.88,32.99,24.99,18.61,13.61,9.75 -selname,lnsp,z,t,u,v,w,q,o3 merge.grib ea5_2020_01_01_00.nc
```

## Other meteorological data sets

The MPTRAC model has also been tested with a number of other meteorological data sets.

The [Modern-Era Retrospective analysis for Research and Applications, Version 2](https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/) (MERRA-2) provides data beginning in 1980. It was introduced to replace the original MERRA dataset because of the advances made in the assimilation system that enable the assimilation of modern hyperspectral radiance and microwave observations, along with GPS-Radio Occultation datasets. It also uses NASA's ozone profile observations that began in late 2004. Additional advances in both the GEOS model and the GSI assimilation system are included in MERRA-2. The spatial resolution remains about the same (about 50 km in the latitudinal direction) as in MERRA.

The [NCEP/NCAR Reanalysis 1](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html) project is using a state-of-the-art analysis/forecast system to perform data assimilation using past data from 1948 to the present. A large subset of this data is available from PSL in its original 4 times daily format and as daily averages. However, the data from 1948-1957 is a little different, in the regular (non-Gaussian) gridded data. That data was done at 8 times daily in the model because the inputs available in that era were available at 3Z, 9Z, 15Z, and 21Z, whereas the 4x daily data has been available at 0Z, 6Z, 12Z, and 18Z. These latter times were forecasted and the combined result for this early era is 8x daily. The local ingestion process took only the 0Z, 6Z, 12Z, and 18Z forecasted values, and thus only those were used to make the daily time series and monthly means here.

The [Global Forecast System (GFS)](https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs) is a weather forecast model produced by the National Centers for Environmental Prediction (NCEP). Dozens of atmospheric and land-soil variables are available through this dataset, from temperatures, winds, and precipitation to soil moisture and atmospheric ozone concentration. The entire globe is covered by the GFS at a base horizontal resolution of 18 miles (28 kilometers) between grid points, which is used by the operational forecasters who predict weather out to 16 days in the future. Horizontal resolution drops to 44 miles (70 kilometers) between grid points for forecasts between one week and two weeks.

## MeteoCloud data archive in Jülich

For registered users and collaboration partners at the Jülich Supercomputing Centre, Germany, a number of meteorological data sets are readily available for MPTRAC in netCDF file format. The data sets are stored in the [MeteoCloud data archive](https://datapub.fz-juelich.de/slcs/meteocloud) in Jülich. Please contact us in order to inquire about getting access to the data.

## Interoperability with CLaMS data 

MPTRAC can read meteorological data that follows the data format guidelines for the CLaMS. This option can be set true with the parameter: CLAMS_MET_DATA = 1.
