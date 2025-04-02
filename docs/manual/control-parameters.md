# Control parameters

The MPTRAC applications are controlled either through a control parameter file or by specifying control parameters directly as command-line arguments. A comprehensive [list of control parameters](https://slcs-jsc.github.io/mptrac/doxygen/structctl__t.html) is available in the Doxygen manual.

By default, if no values are explicitly provided, the applications will use predefined default values for most control parameters. However, it's important to carefully review the log output from the MPTRAC tools to ensure that the control parameters are set as intended.

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

Note that the blanks before and after the equal sign are mandatory. Array indices start counting from zero, i.e, a[0], a[1], ..., as in the C programming language. You can use a minus sign to indicate that no control parameter file is being used.

## Command line arguments

Control parameters can also be specified via command line arguments:

```
./atm_init trac.ctl data/atm_init.tab \
               INIT_LON0 -72.117 INIT_LON1 -72.117 \
               INIT_LAT0 -40.59 INIT_LAT1 -40.59
```

Command line arguments take precedence over values provided in the control parameter file.
