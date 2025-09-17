# wind 

This app creates a meteorological data file with synthetic wind fields suitable for testing advection using the method of Williamson et al. (1992, Sect. 3). The synthetic wind field describes a global rotation around the Earth with a given rotation axis. The output is written in form of MPTRAC's netCDF meteo file format.

* Williamson, D. L., Drake, J. B., Hack, J. J., Jakob, R., & Swarztrauber, P. N. (1992). A standard test set for numerical approximations to the shallow water equations in spherical geometry. Journal of computational physics, 102(1), 211-224.

```
# Calling wind
$ ./wind <ctl> <metbase>
```

The required parameters are:
* ctl: the control parameter file
* metbase: the basename of the output file

Optional control parameters:
* WIND_T0: time step of meteo data (default: 0 s)
* WIND_NX: number of longitudes (default: 360)
* WIND_NY: number of latitudes (default: 181)
* WIND_NZ: number of pressure levels (default: 61)
* WIND_Z0: log-pressure height of lowermost pressure level (default: 0 km)
* WIND_Z1: log-pressure height of uppermost pressure level (default: 60 km)
* WIND_U0: horizontal wind speed at lowermost level (default: 38.587660177302 m/s)
* WIND_U1: horizontal wind speed at uppermost level (default: 38.587660177302 m/s)
* WIND_W0: vertical velocity (default: 0 m/s)
* WIND_ALPHA: rotation angle (default: 0 deg)
* WIND_LAT_REVERSE: latitude direction flag (default: 0, set to 1 for reversed latitude direction)
