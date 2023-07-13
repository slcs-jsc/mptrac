# tropo

This application creates a tropopause climatology from meteorological data.

```
# calling tropo
$ ./tropo  <ctl> <tropo.nc> <met0> [<met1> ...]
```
The required arguments are:
* ctl: The file containing the control parameters.
* tropo.nc: The output file in which the climatology should be stored.
* met0 \[\<met1\> ...\]: The meteorological data from which the climatology should be created.
