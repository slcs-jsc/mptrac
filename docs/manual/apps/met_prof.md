# met_prof

Extract vertical profile from meteorological data.

```
# calling met_prof
$ ./met_prof  <ctl> <prof.tab> <met0> [<met1> ...]
```

The required arguments are:

* ctl: The file containing the control parameters.
* prof.tab: The output file in which the profile data should be stored.
* \<met0\> \[\<met1\> ...]: The meteorological input data.

Calling met_prof and creating a 1 km resolution verical profile from a global ERA5 data sample of 1° x 1° would look like as follows:

```
# Calling met_map
$./met_prof - prof_era5_2017_01_08_17_1_1.tab era5_2017_01_08_17.nc PROF_LON0 -180 PROF_LON1 180\
	PROF_LAT0 -90 PROF_LAT1 90 PROF_DLON 1 PROF_DLAT 1 PROF_Z0 0 PROF_Z1 25 PROF_DZ 1
```

The optional control parameters required are listed as follows:

| parameter | purpose | default | 
|:-----------|:---------|---------:|
| PROF_Z0 | Altitude range start [km] | -999 |
| PROF_Z1 | Altitude range end [km] | -999 |
| PROF_DZ | Altitude interval [km]| -999 |
| PROF_LON0 | Longitude range start | 0 |
| PROF_LON1 | Longitude range end | 0 |
| PROF_DLON | Longitude interval | -999 |
| PROF_LAT0 | Latitude range start | 0 |
| PROF_LAT1 | Latitude range end | 0 |
| PROF_DLAT | Latitude interval | -999 |
