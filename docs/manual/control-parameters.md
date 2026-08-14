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

This mechanism is also useful for short sensitivity studies. For example, a release layer can be defined in relative PBL coordinates without editing the control file permanently:

```
./atm_init trac.ctl data/atm_init.tab \
               INIT_T0 360547200 INIT_T1 360547200 \
               INIT_RELZ0 0.0 INIT_RELZ1 1.2 INIT_DRELZ 0.2
```

## Radionuclide deposition

Radionuclide decay and deposition require the activity quantities to be part
of the air-parcel data. The recognized quantity names are `Apb210`, `Abe7`,
`Acs137`, and `Ai131`; their values are activities in Bq per air parcel.

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
