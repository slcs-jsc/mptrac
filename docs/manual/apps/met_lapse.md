# met_lapse

`met_lapse` calculates lapse-rate statistics from one or more
meteorological data files.

```bash
./met_lapse <ctl> <lapse.tab> <met0> [<met1> ...]
```

Required arguments:

- `ctl`: control parameter file.
- `lapse.tab`: output file for the lapse-rate statistics.
- `met0 [<met1> ...]`: meteorological input files.

Set optional `LAPSE_*` and `MET_*` parameters in the control file for
reliable operation. Appending them after the meteorological input files
is not recommended for this tool.

Relevant control parameters:

| parameter | purpose | default |
|:----------|:--------|--------:|
| `LAPSE_DZ` | vertical layer thickness for lapse-rate estimates [0.1 km levels] | 20 |
| `LAPSE_LAT0` | latitude range start [deg] | -90 |
| `LAPSE_LAT1` | latitude range end [deg] | 90 |
| `LAPSE_Z0` | altitude range start [km] | 0 |
| `LAPSE_Z1` | altitude range end [km] | 100 |
| `LAPSE_INTPOL` | use spline interpolation (`1`) or linear interpolation (`0`) | 1 |

The output columns are:

- `$1`: mean altitude [km]
- `$2`: mean latitude [deg]
- `$3`: lapse-rate bin [K/km]
- `$4`: counts of maxima per bin
- `$5`: total number of maxima
- `$6`: normalized frequency of maxima
- `$7`: counts of minima per bin
- `$8`: total number of minima
- `$9`: normalized frequency of minima
- `$10`: counts of means per bin
- `$11`: total number of means
- `$12`: normalized frequency of means
- `$13`: counts of standard deviations per bin
- `$14`: total number of standard deviations
- `$15`: normalized frequency of standard deviations

Example:

```bash
./met_lapse lapse.ctl lapse_2011_06_05.tab ei_2011_06_05_00.nc
```
