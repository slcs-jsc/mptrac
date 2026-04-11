# met_check_dt

`met_check_dt` checks model time-step constraints for a meteorological
data file. It reports height-dependent estimates for horizontal and
vertical advection and diffusion time steps.

```bash
./met_check_dt <ctl> <dt_file> <met> [KEY VALUE ...]
```

Required arguments:

- `ctl`: control parameter file.
- `dt_file`: output file for the time-step diagnostics.
- `met`: meteorological input file.

Optional trailing arguments can be used to override control parameters
on the command line as `KEY VALUE` pairs.

Relevant control parameters:

| parameter | purpose | default |
|:----------|:--------|--------:|
| `KX` | horizontal diffusivity [m2/s] | 50.0 |
| `KZ` | vertical diffusivity [m2/s] | 0.1 |
| `DX` | horizontal grid spacing [km] | required |
| `CMAX` | maximum Courant number for advection | 0.5 |
| `NMAX` | maximum stability factor for diffusion | 0.3 |

The output columns are:

- `$1`: height [km]
- `$2`: time step for horizontal advection [s]
- `$3`: time step for vertical advection [s]
- `$4`: time step for horizontal diffusion [s]
- `$5`: time step for vertical diffusion [s]

Example:

```bash
./met_check_dt trac.ctl dt_2011_06_05.tab ei_2011_06_05_00.nc DX 50
```

