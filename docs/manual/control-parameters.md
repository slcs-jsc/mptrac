# Control parameters

The MPTRAC apps are controlled by means of a control parameter file or by specifying control parameters as command line arguments. A complete [list of control parameters](https://slcs-jsc.github.io/mptrac/structctl__t.html) can be found in the Doxygen manual.

In most cases, default values for control parameters will be used, if no value is explicitly specified. Please carefully check the log output of the MPTRAC tools to make sure the control parameter values are selected as desired.

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

Note that blanks before and after the equal sign are mandatory. Array indices start counting from zero, i.e, a[0], a[1], ..., like in C. You can use the minus sign to indicate that no control parameter file is being used.

## Command line arguments

Control parameters can also be specified via the command line:

```
./atm_init trac.ctl data/atm_init.tab \
               INIT_LON0 -72.117 INIT_LON1 -72.117 \
               INIT_LAT0 -40.59 INIT_LAT1 -40.59
```

Command line arguments have priority over the values given in the control parameter file.
