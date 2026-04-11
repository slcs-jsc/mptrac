# tropo_clim

`tropo_clim` calculates tropopause climatology statistics from one or
more tropopause data files.

```bash
./tropo_clim <ctl> <clim.tab> <var> <tropo.nc> [<tropo2.nc> ...]
```

Required arguments:

- `ctl`: control parameter file.
- `clim.tab`: output file for the tropopause climatology.
- `var`: tropopause variable prefix to process.
- `tropo.nc [<tropo2.nc> ...]`: tropopause input files.

Supported values for `var` are the variable prefixes written by
`tropo`, for example:

| var | description |
|:----|:------------|
| `clp` | cold point tropopause |
| `dyn` | dynamical tropopause |
| `wmo_1st` | WMO first lapse-rate tropopause |
| `wmo_2nd` | WMO second lapse-rate tropopause |

For a selected variable prefix such as `wmo_1st`, the input file must
contain the variables `wmo_1st_z`, `wmo_1st_p`, and `wmo_1st_t`.
Water vapor and ozone variables with the suffixes `_q` and `_o3` are
used if they are present.

The output columns are:

- `$1`: longitude [deg]
- `$2`: latitude [deg]
- `$3`: tropopause height mean [km]
- `$4`: tropopause pressure mean [hPa]
- `$5`: tropopause temperature mean [K]
- `$6`: tropopause water vapor mean [ppv]
- `$7`: tropopause ozone mean [ppv]
- `$8`: tropopause height standard deviation [km]
- `$9`: tropopause pressure standard deviation [hPa]
- `$10`: tropopause temperature standard deviation [K]
- `$11`: tropopause water vapor standard deviation [ppv]
- `$12`: tropopause ozone standard deviation [ppv]
- `$13`: number of valid data points
- `$14`: occurrence frequency [%]

Example:

```bash
for var in clp dyn wmo_1st wmo_2nd ; do
  ./tropo_clim tropo.ctl clim_${var}.tab ${var} tropo_2011_06_05.nc
done
```

