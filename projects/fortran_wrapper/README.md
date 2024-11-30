# Fortran Wrapper

This project shows the integration of C functions in Fortran. The project started as natESM sprint LAGOOn (url will be added, when published).

To integrate multi language programming with Fortran and C in a single program the ISO_C_BINDINGS <!-- https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fC_005fBINDING.html --> module is since Fortran 2003 the standard tool that faciliates interoperability between the two languages. In Fortran the Compiler is informed with the BIND(C) attribute that a symbol shall be interoperable with C.

To compile the Fortran wrapper, navigate to the `src` directory and run the following command:
```
$ make wrapper
```

To test the compiled wrapper, use the following command in the same directory:
```
$ make wrapper_test
```

## Application

### trac_fortran.f90

This program calculates trajectories based on MPTRAC's module_advect function. It is the Fortran counterpart to trac.c. The program is called with at least four (4) arguments:

```
# calling trac_fortran
$ ./trac_fortran  <dirlist> <ctl> <atm_in> <metbase>
```

The required arguments are:

* dirlist: A file containing directories to be processed. Each directory has to have an own control parameter file and a starting point file.

* ctl: In the control parameter file the configuration parameters can be set.

* atm_in: The starting point file contains a list of starting points for the trajectory calculation.

* metbase: Here, the path to the meteorological data files and their basename shall be given.

Example from the repository:

```
# minimal required input to run trac_fortran
$ ./trac_fortran data/dirlist trac.ctl atm_split.tab meteo/ei
```

In addition it is possible to append control parameters, giving the flag name first, followed by the parameter, all separated by a space. Example:

```
#
$ ./trac_fortran data/dirlist trac.ctl atm_split.tab meteo/ei ATM_BASENAME atm_diff GRID_BASENAME grid_diff
```

## Interface

### mptrac_fortran.f90

To ensure the interoperabilty between Fortran and C an interface is needed. This interface includes some general dimension variables, the structures atm_t, clim_photo_t, clim_t, clim_ts_t, clim_zm_t, clt_t, met_t and the functions mptrac_get_met, mptrac_module_advect, mptrac_module_timesteps, mptrac_read_atm, mptrac_read_clim, mptrac_read_ctl, mptrac_read_met, mptrac_write_output. The functions have always the prefix mptrac to indicate that they are only an interface calling the original MPTRAC function.

### Checking order and array sizes

It is crucial that the order and array sizes in the Fortran interface match those in the original C structure. The structures in MPTRAC consist of extensive variable lists (up to O100 for the control parameter list) and data types that may also depend on self-defined structures (struct of struct). A shell script find_vars.sh (in tests/wrapper_test/) can be used to check differences in the structures and variable dimensions. This script is automatically executed when the wrapper test is performed.

## Example

An example of using the Fortran wrapper for trajectory calculations can be found in tests/wrapper_test/run.sh

This file executes the trajectory calculation first with the original C code trac.c. The output is stored in the directory data.ref. Afterwards the Fortran code trac_fortran.f90 is called. This output is stored in the directory data. Both directories are compared to make sure that the Fortran version produces the same results as the C version.
