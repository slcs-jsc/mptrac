Extract vertical profile from meteorological data.

```
# calling met_prof
$ ./met_prof  <ctl> <prof.tab> <met0> [<met1> ...]
```

The required arguments are:

* ctl: The file containing the control parameters.
* prof.tab: The output file in which the profile data should be stored.
* \<met0\> \[\<met1\> ...]: The meteorological input data.

The optional control parameters required are listed as follows:

| parameter | purpose | default | 
|:-----------|:---------|---------:|
| PROF_Z0 | Altitude range start [km] | -999 |
| PROF_Z1 | Altitude range end [km] | -999 |
| PROF_Z0 | Altitude interval [km]| -999 |
| PROF_LON0 | Longitude range start | 0 |
| PROF_LON1 | Longitude range end | 0 |
| PROF_DLON | Longitude interval | -999 |
| PROF_LAT0 | Latitude range start | 0 |
| PROF_LAT1 | Latitude range end | 0 |
| PROF_DLAT | Latitude interval | -999 |
