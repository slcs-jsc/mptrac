# Tropopause data

This page provides information on the methods used to calculate and extract tropopause data with MPTRAC.

## Determination of the tropopause

In general, the tropopause data in MPTRAC are calculated from the given meteorological input data, i.e, the data are not read from the meteo data files.

For the determination of the height level of the tropopause, the meteorological data are first interpolated to a fine vertical grid by means of cubic spline interpolation. The tropopause level is estimated based on the interpolated data. This approach improves the vertical resolution of the tropopause data. However, the vertical resolution of the meteorological input data will also play an important role.

The method used to extract the height level of the tropopause is selected by means of the parameter MET_TROPO. The following options to determine the tropopause level have been implemented:

| MET_TROPO   | Description       |
| ----------- | ----------------- |
| 0           | no tropopause is calculated |
| 1           | use NCEP/NCAR Reanalysis-1 monthly mean zonal mean climatology |
| 2           | extract the cold point (temperature minimum along the vertical profile) |
| 3 (default) | use WMO definition based on thermal lapse rate |
| 4           | extract location of the WMO second tropopause |
| 5           | use dynamical definition (3.5 PVU in the extratropics and 380 K in the tropics) |

## Tropopause data along trajectories

The trac tool of MPTRAC can extract the tropopause pressure along the trajectories.

The tropopause pressure needs to be selected as quantity "pt", for example:

```
    NQ = 1
    QNT_NAME[0] = pt
```

The default unit is "hPa".

## Tropopause data files

MPTRAC provides the tool "tropo" that generates tropopause data files from given meteorological input data.

The tropopause data are provided as netCDF files. Each data file may contain one or more time steps of the reanalysis.

The data files provide geopotential height, pressure, temperature, and water vapor volume mixing ratios for the cold point, WMNO 1st and 2nd tropopause, and the dynamical tropopause.

The tool "tropo_sample" can be used to extract tropopause data from the netCDF files created by "tropo" for a set of locations specified as a particle data file (atm.tab file).

Tropopause data files for various meteorological reanalyses can be found in the meteocloud data archive in Jülich:

```
    /p/fastdata/slmet/slmet111/met_data/tropo
```

## References

WMO tropopause definition: https://library.wmo.int/doc_num.php?explnum_id=6960

NCEP/NCAR Reanalysis-1: https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.pressure.html
