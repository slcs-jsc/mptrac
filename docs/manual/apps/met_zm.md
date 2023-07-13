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
