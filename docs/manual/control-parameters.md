# Control parameters

The MPTRAC applications are controlled either through a control
parameter file or by specifying control parameters directly as
command-line arguments. A comprehensive
[list of control parameters](https://slcs-jsc.github.io/mptrac/doxygen/structctl__t.html)
is available in the Doxygen manual. The Doxygen page is the reference
for the `ctl_t` fields themselves, while this manual page explains the
file syntax, command-line overrides, and practical usage patterns.

By default, if no values are explicitly provided, the applications
will use predefined default values for most control
parameters. However, it's important to carefully review the log output
from the MPTRAC tools to ensure that the control parameters are indeed
set as intended.

## Control parameter file

Example of a control parameter file:

```
# Quantities...
NQ = 3
QNT_NAME[0] = t
QNT_NAME[1] = u
QNT_NAME[2] = v

# Meteo data...
METBASE = meteo/ei
DT_MET = 86400

# Grid output...
GRID_LON0 = -90
GRID_LON1 = 60
GRID_LAT0 = -60
GRID_LAT1 = -15
GRID_NX = 300
GRID_NY = 90
```

Note that the blanks before and after the equal sign are
mandatory. Array indices start counting from zero, i.e, `a[0]`, `a[1]`,
..., as in the C programming language. You can use a minus sign (`-`) to
indicate that no control parameter file is being used.

## Command line arguments

Control parameters can also be specified via command line arguments:

```
./atm_init trac.ctl data/atm_init.tab \
               INIT_LON0 -72.117 INIT_LON1 -72.117 \
               INIT_LAT0 -40.59 INIT_LAT1 -40.59
```

Command line arguments take precedence over values provided in the
control parameter file.

For enabled periodic particle and diagnostic output, `ATM_DT_OUT`,
`BUDGET_DT_OUT`, `DEPO_DT_OUT`, `CSI_DT_OUT`, `ENS_DT_OUT`, `GRID_DT_OUT`,
and `VTK_DT_OUT` must be finite and positive, in seconds. A basename of `-`
disables the corresponding output without requiring its interval. Positive
intervals should align with `DT_MOD`; enabled snapshot outputs also write at
`T_STOP`. `CSI_DT_OUT` selects CSI verification snapshot times, with a final
evaluation of the `T_STOP` observation window when points are available.
`MET_DT_OUT <= 0` instead disables sampling of meteorological particle
quantities; a nonfinite value is rejected
when such quantities are configured.

`SAMPLE_DZ > 0` specifies a vertical **half-width** in kilometres around each
sample observation, using log-pressure altitude. The sample VMR uses the air
column within that layer, clipped to surface pressure and the meteorological
model top. `SAMPLE_DZ <= 0` disables vertical filtering and leaves sample VMR
undefined.

This mechanism is also useful for short sensitivity studies. For example, a release layer can be defined in relative PBL coordinates without editing the control file permanently:

```
./atm_init trac.ctl data/atm_init.tab \
               INIT_T0 360547200 INIT_T1 360547200 \
               INIT_RELZ0 0.0 INIT_RELZ1 1.2 INIT_DRELZ 0.2
```

For `ISOSURF=1`, `2`, or `3`, include `iso_var` among the configured
quantities. It stores each particle's initial isosurface reference and is
included in particle output and MPI domain-decomposition exchanges. Its
units depend on `ISOSURF`: hPa for pressure (`1`), hPa/K for the
pressure-to-temperature ratio (`2`), and K for potential temperature (`3`).
Existing particle files without this quantity need a configuration that
provides a free quantity slot; `NQ` and other quantity indices are not
changed automatically. `ISOSURF=0` and `4` do not require it.

## PBL turbulence

With `DIFFUSION=1` and `TURB_PBL_SCHEME=1`, the Hanna-PBL closure uses
along-wind and cross-wind velocity fluctuations and cannot run with
mesoscale Langevin turbulence. In this configuration, `TURB_MESOX` and
`TURB_MESOZ` default to zero; an explicit nonzero value is rejected.
When Hanna-PBL is disabled, both retain their default of 0.16. Setting
`DIFFUSION=0` leaves these turbulence settings inactive. Hanna-PBL
remains compatible with fixed-`K` diffusion outside the PBL.
Hanna turbulence and vertical fixed-`K` diffusion reflect at the PBL top
from their respective sides. Changing PBL height, resolved advection,
and convection can still change a particle's regime.
`TURB_PBL_TRANS` retains its Fixed-`K`-only behavior; with Hanna enabled it
only blends Fixed-`K` diffusivity above the reflecting boundary.
`DT_PBL` sets the maximum internal Hanna-PBL timestep in seconds
(default `2`). It must be finite and positive when Hanna-PBL diffusion
is active. All particles use a common number of substeps, while each
retains its signed individual timestep. Smaller values increase the
temporal resolution and computational cost.

## Long-lived tracer chemistry and decay

`TRACER_CHEM` accepts `0` (default, disabled) or `1`. When enabled with
`Cccl4`, `Cccl3f`, `Cccl2f2`, or `Cn2o` configured, it requires the
photolysis file `CLIM_PHOTO` and O($^1$D) file `CLIM_O1D_FILENAME`.
Their default paths point to the bundled climatologies. It updates
those species-specific mixing ratios using signed particle timesteps.
Prescribed lower-boundary time series can reset these concentrations
where the boundary criteria apply, independently of chemistry.

`TDEC_TROP` and `TDEC_STRAT` are generic e-folding lifetimes in seconds.
Both default to zero, disabling generic decay. Set both to finite,
positive values to enable it; a single positive lifetime is invalid.
Generic decay updates `m` and/or `vmr` and, when configured, mass-loss
diagnostics. It is independent of `TRACER_CHEM` and `RADIO_DECAY`.

## Age-of-air reference clock

Configure `aoa` as a particle quantity to store a reference timestamp in
seconds since 2000-01-01 00:00 UTC. `BOUND_LAT0`, `BOUND_LAT1`, `BOUND_P0`,
and `BOUND_P1` define the region where it is reset to particle time; optional
`BOUND_DPS`, `BOUND_DZS`, `BOUND_ZETAS`, and `BOUND_PBL` further restrict that
region. Compute elapsed age externally as observation time minus `aoa`.

## OH and H2O2 chemistry

`OH_CHEM_REACTION=0` disables OH chemistry. Modes `1`, `2`, and `3`
select a constant bimolecular coefficient, an Arrhenius coefficient,
and termolecular falloff, respectively. `OH_CHEM[0]` is the coefficient
in mode 1; modes 2 and 3 also use `[1]`, and mode 3 uses `[2]` and `[3]`
for the high-pressure limit and its temperature exponent. The two
limiting coefficients in mode 3 must be positive. The active OH modes
require a valid `CLIM_OH_FILENAME`; `OH_CHEM_BETA` controls optional
diurnal variation. These kinetics may be configured for species other
than SO$_2$.

`H2O2_CHEM_REACTION=1` selects the SO$_2$-specific aqueous oxidation
approximation and requires a valid `CLIM_H2O2_FILENAME`; `0` disables it.
If quantity `Cx` is configured, the empirical high-SO$_2$ correction
applies for `SPECIES=SO2` or unspecified `SPECIES`. It does not apply
to an explicitly selected different species. `Cx` stores a dimensionless
mole fraction, not a numerical ppbv value. The correction coefficients
are calibrated for that mole-fraction input.

## Wet deposition

In-cloud and below-cloud scavenging can be configured independently.
For each branch, a positive `WET_DEPO_IC_A` or `WET_DEPO_BC_A` selects
empirical aerosol scavenging and takes precedence over the corresponding
positive Henry-law coefficient `WET_DEPO_IC_H[0]` or `WET_DEPO_BC_H[0]`.
If neither coefficient is positive, that branch is disabled. The Henry-law
temperature coefficients are `[1]`; `SPECIES` supplies defaults that can
be overridden. `WET_DEPO_PRE[0]=0.5` and `WET_DEPO_PRE[1]=0.36` define
the empirical precipitation proxy from column condensate. Radioactive
wet deposition uses these precipitation coefficients independently of
general wet-deposition activation. A positive `WET_DEPO_SO2_PH` applies
the SO$_2$ pH correction only to in-cloud Henry-law scavenging; when a
different species is named, this combination is rejected.

## Dry deposition

`DRY_DEPO_MODE` selects dry deposition: `0` disables it (default), `1`
uses the positive prescribed velocity `DRY_DEPO_VDEP` in m/s, and `2`
calculates aerosol velocity from the particle quantities `rp` and `rhop`.
`DRY_DEPO_VDEP` is ignored in mode 2 and does not activate deposition.
Both active modes and `RADIO_DEPO` use `DRY_DEPO_DZ` as the geometric
surface-layer thickness in metres (default 30 m). This replaces the old
pressure-based `DRY_DEPO_DP` parameter; old values cannot be reused as metres.

## Radionuclide deposition

Radionuclide activities are stored in Bq per air parcel. `RADIO_DECAY`
recognizes `Arn222`, `Apb210`, `Abe7`, `Acs137`, `Ai131`, and `Axe133`;
`RADIO_DEPO` applies only to `Apb210`, `Abe7`, `Acs137`, and `Ai131`.
Both switches accept only `0` or `1`. Deposition uses positive individual
timesteps only; decay also uses negative timesteps mathematically.

| Parameter | Meaning | Default |
|-----------|---------|---------|
| `RADIO_DECAY` | Apply radioactive decay to airborne and deposited activity | `0` |
| `RADIO_DEPO` | Apply dry and wet radionuclide deposition | `0` |
| `DEPO_BASENAME` | Basename of cumulative ground-inventory files; `-` disables output | `-` |
| `DEPO_DT_OUT` | Deposition output interval in seconds | `86400` |
| `DEPO_TYPE` | Output format: `0` for ASCII, `1` for netCDF | `0` |

The horizontal deposition grid is configured with `GRID_LON0`, `GRID_LON1`,
`GRID_NX`, `GRID_LAT0`, `GRID_LAT1`, and `GRID_NY`. `GRID_BASENAME` does not
need to be enabled. `GRID_NC_LEVEL` controls netCDF compression for deposition
files as it does for regular grid output.

`RADIO_DEPO` currently requires `MET_COORD_TYPE=0` and is not supported in
domain-decomposition (DD) builds. Invalid combinations are rejected while
reading the control parameters.

Example:

```text
NQ = 4
QNT_NAME[0] = Apb210
QNT_NAME[1] = Abe7
QNT_NAME[2] = Acs137
QNT_NAME[3] = Ai131

RADIO_DECAY = 1
RADIO_DEPO = 1
DEPO_BASENAME = depo
DEPO_DT_OUT = 43200
DEPO_TYPE = 1

GRID_LON0 = -90
GRID_LON1 = -50
GRID_NX = 80
GRID_LAT0 = -55
GRID_LAT1 = -25
GRID_NY = 60
```
