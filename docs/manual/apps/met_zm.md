# met_zm

Extract zonal mean from meteorological data.

The calling sequence is:

```
# calling met_zm
$ ./met_zm  <ctl> <zm.tab> <met0> [<met1> ...]
```

The required arguments are:
* ctl: The file containing the control parameters.
* zm.tab: The output file where the zonal means are stored.
* met0 \[\<met1\> ...\]: The meteorological input files from which the zonal mean should be calculatd from.

Calling met_zm and creating a global zonal mean data with latitude band of 10 degree from ERA5 data would look like as follows:

```
# Calling met_zm
$./met_prof - zm_era5_2017_01_08_17_1_1.tab era5_2017_01_08_17.nc \
	ZM_DLAT 10  
```

The optional control parameters required are listed as follows:

| parameter | purpose | default | 
|:-----------|:---------|---------:|
| ZM_Z0 | Altitude range start [km] | -999 |
| ZM_Z1 | Altitude range end [km] | -999 |
| ZM_DZ | Altitude interval [km]| -999 |
| ZM_LON0 | Longitude range start | -180 |
| ZM_LON1 | Longitude range end | 180 |
| ZM_LAT0 | Latitude range start | -90 |
| ZM_LAT1 | Latitude range end | 90 |
| ZM_DLAT | Latitude interval | -999 |

**Note:** if ZM_DZ, ZM_Z0, ZM_Z1 are not defined, the output dimension of altitude will follow the meteorological data.