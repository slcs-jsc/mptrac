# wind

`wind` creates synthetic meteorological fields in MPTRAC's netCDF format.
The wind field combines global rotation with optional spatial perturbations.

## Global rotation

The base wind describes a global rotation around an axis set by
`WIND_ALPHA`, following Williamson et al. (1992, Sect. 3). `WIND_U0` and
`WIND_U1` set the speed at the lowest and highest levels, with linear
interpolation between them. At `WIND_ALPHA=0`, the flow is zonal:
`u = U cos(latitude)` and `v = 0`. `WIND_W0` sets the base vertical velocity.

## Optional spatial perturbations

`WIND_PERT_U`, `WIND_PERT_V`, and `WIND_PERT_W` add smooth, deterministic
variability to the three wind components. `WIND_PERT_K` sets their spatial
frequency. The perturbations vary with longitude, latitude, and altitude;
they are periodic at 0°/360° longitude and taper to zero at the poles.
Each amplitude can be set independently. All three default to zero, so the
global rotation remains unchanged unless perturbations are enabled.

The generator also writes configurable auxiliary meteorological fields
(`PS`, `BLH`, `ISHF`, `IEWS`, `INSS`, `Q`, etc.) through MPTRAC's standard
`write_met_nc()` path. `WIND_BLH` is converted to boundary-layer top
pressure. `WIND_Q` and `WIND_O3` are specified as mass mixing ratios and
converted to the representation used by the meteo writer.

Williamson, D. L., Drake, J. B., Hack, J. J., Jakob, R., & Swarztrauber,
P. N. (1992). A standard test set for numerical approximations to the
shallow water equations in spherical geometry. *Journal of Computational
Physics*, 102(1), 211–224.

## Usage

```bash
./wind <ctl> <metbase>
```

The required parameters are:

* ctl: the control parameter file
* metbase: the basename of the output file

## Optional control parameters

* WIND_T0: time step of meteo data (default: 0 s)
* WIND_NX: number of longitudes (default: 360, minimum: 2)
* WIND_NY: number of latitudes (default: 181, minimum: 2)
* WIND_NZ: number of pressure levels (default: 61, minimum: 2)
* WIND_Z0: log-pressure height of lowermost pressure level (default: 0 km)
* WIND_Z1: log-pressure height of uppermost pressure level (default: 60 km)
* WIND_U0: horizontal wind speed at lowermost level (default: 38.587660177302 m/s)
* WIND_U1: horizontal wind speed at uppermost level (default: 38.587660177302 m/s)
* WIND_W0: vertical velocity (default: 0 m/s)
* WIND_PERT_U: amplitude of the zonal wind perturbation (default: 0 m/s)
* WIND_PERT_V: amplitude of the meridional wind perturbation (default: 0 m/s)
* WIND_PERT_W: amplitude of the vertical wind perturbation (default: 0 m/s)
* WIND_PERT_K: positive integer number of perturbation cycles around the globe (default: 24, corresponding to a 15° wavelength)
* WIND_ALPHA: rotation angle (default: 0 deg)
* WIND_LAT_REVERSE: latitude direction flag (default: 0, set to 1 for reversed latitude direction)
* WIND_TEMP0: temperature at the lowermost level (default: 280 K)
* WIND_TEMP1: temperature at the uppermost level (default: 280 K)
* WIND_PS: surface pressure in hPa (default: 1013.25 hPa)
* WIND_ZS: surface height (default: 0 km)
* WIND_T2M: 2 m temperature (default: 280 K)
* WIND_IEWS: eastward turbulent surface stress (default: 0 N/m^2)
* WIND_INSS: northward turbulent surface stress (default: 0 N/m^2)
* WIND_ISHF: surface sensible heat flux in ECMWF convention (positive downward; default: 0 W/m^2)
* WIND_LSM: land-sea mask (default: 1)
* WIND_SST: sea surface temperature (default: 280 K)
* WIND_BLH: boundary layer height (default: 1 km). Converted to boundary layer pressure via the barometric formula (P(z) in hPa, where z = WIND_ZS + WIND_BLH)
* WIND_Q: specific humidity (default: 0 kg/kg)
* WIND_O3: ozone mass mixing ratio (default: 0 kg/kg)

## Perturbation formula

The perturbations are deterministic and are applied only when at least one
amplitude is nonzero. With longitude `lambda`, latitude `phi`, and level index
`i` (from 0 to `WIND_NZ - 1`), define `T = cos(phi)^2` and
`H = 0.75 + 0.25 cos(2 pi i / (WIND_NZ - 1))`. Angles are in radians. The
three-dimensional wind receives these additions before output:

```
u += WIND_PERT_U * H * T * sin(K * lambda) * cos(K * phi)
v += WIND_PERT_V * H * T * cos(K * lambda) * cos(K * phi + pi/4)
w += WIND_PERT_W * H * T * sin(K * lambda + pi/4) * sin(K * phi + pi/4)
```

The same horizontal perturbations, with `H = 1`, are added to surface wind.
The vertical perturbation is in m/s and is converted to pressure velocity
along with `WIND_W0`. The longitude functions are periodic at 0°/360°; the
latitude taper makes the perturbations vanish smoothly at both poles. The
default zero amplitudes preserve the original wind field.
