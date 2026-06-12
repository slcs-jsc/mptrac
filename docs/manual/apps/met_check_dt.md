# met_check_dt

`met_check_dt` checks model time-step constraints for a meteorological
data file. It reports height-dependent estimates for horizontal and
vertical advection and diffusion time steps using the current MPTRAC
control parameters.

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
| `TURB_DX_PBL` | PBL horizontal diffusivity [m2/s] used for column `$4` | 50 |
| `TURB_DX_TROP` | free-tropospheric horizontal diffusivity [m2/s] used for column `$4` | 50 |
| `TURB_DX_STRAT` | stratospheric horizontal diffusivity [m2/s] used for column `$4` | 0 |
| `TURB_DZ_PBL` | PBL vertical diffusivity [m2/s] used for column `$5` | 0 |
| `TURB_DZ_TROP` | free-tropospheric vertical diffusivity [m2/s] used for column `$5` | 0 |
| `TURB_DZ_STRAT` | stratospheric vertical diffusivity [m2/s] used for column `$5` | 0.1 |
| `TURB_PBL_TRANS` | fractional depth of the PBL transition layer used for column `$6` | 0 |
| `DX` | horizontal grid spacing [km] | required |
| `CMAX` | maximum Courant number for advection | 0.5 |
| `NMAX` | maximum stability factor for diffusion | 0.3 |

The output columns are:

- `$1`: height [km]
- `$2`: time step for horizontal advection [s]
- `$3`: time step for vertical advection [s]
- `$4`: time step for horizontal diffusion [s] using layer-aware `TURB_DX_*`
- `$5`: time step for vertical diffusion [s] using layer-aware `TURB_DZ_*`
- `$6`: time step for resolving diffusion across the PBL transition layer [s]
- `$7`: time step for PBL diffusion on the full PBL-depth scale [s]

Columns `$4` and `$5` follow the same PBL/troposphere/stratosphere
weighting used by the fixed-`K` diffusion module. Column `$6` is the
stricter transition-layer criterion and is only finite within the
configured PBL transition layer. Column `$7` is a less conservative
PBL-scale criterion and is only finite within the PBL itself.

Example:

```bash
./met_check_dt trac.ctl dt_2011_06_05.tab ei_2011_06_05_00.nc DX 50 TURB_DZ_PBL 1.0 TURB_DZ_TROP 0.1 TURB_DZ_STRAT 0.1 TURB_PBL_TRANS 0.1
```
