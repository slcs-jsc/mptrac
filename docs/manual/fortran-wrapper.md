# Fortran wrapper

MPTRAC has been equiped with a wrapper to access the C functions from
Fortran. The project started as a coding sprint in the
[natESM project LAGOOn](https://www.nat-esm.de/services/support-through-sprints/documentation).

To integrate multi language programming with Fortran and C in a single
program the [ISO_C_BINDINGS](https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fC_005fBINDING.html)
module is since Fortran 2003 the standard approach that facilitates
interoperability between the two languages. In Fortran, the compiler
is informed with the BIND(C) attribute that a symbol shall be
interoperable with C.

## Compilation and testing

To compile the Fortran wrapper, navigate to the `src` directory and
run the following command:

```
$ make wrapper
```

To test the compiled wrapper, use the following command in the same
directory:

```
$ make wrapper_test
```

## Application

### trac_fortran.f90

This program calculates trajectories by calling MPTRAC's time step
routine. It is the Fortran counterpart to trac.c. The program is called
with a directory list, a control file name, and an atmospheric input
file name:

```
# calling trac_fortran
$ ./trac_fortran <dirlist> <ctl> <atm_in> [KEY VALUE ...]
```

The required arguments are:

* `dirlist`: A file containing directories to be processed. Each
  directory has to have an own control parameter file and a starting
  point file.

* `ctl`: In the control parameter file the configuration parameters can
  be set.

* `atm_in`: The starting point file contains a list of starting points
  for the trajectory calculation.

Example from the repository:

```
# minimal required input to run trac_fortran
$ ./trac_fortran data/dirlist trac.ctl atm_split.tab
```

In addition it is possible to append control parameters, giving the
flag name first, followed by the parameter, all separated by a
space. Example:

```
#
$ ./trac_fortran data/dirlist trac.ctl atm_split.tab ATM_BASENAME atm_diff GRID_BASENAME grid_diff
```

## Interface

### mptrac_fortran.f90

To ensure the interoperability between Fortran and C an interface is
needed. This interface includes some general dimension variables, the
structures `atm_t`, `cache_t`, `clim_photo_t`, `clim_t`, `clim_ts_t`,
`clim_zm_t`, `ctl_t`, `dd_t`, `met_t` and the functions
`mptrac_alloc`, `mptrac_free`, `mptrac_get_met`, `mptrac_init`,
`mptrac_read_atm`, `mptrac_read_clim`, `mptrac_read_ctl`,
`mptrac_read_met`, `mptrac_run_timestep`, `mptrac_write_atm`,
`mptrac_write_met`, `mptrac_write_output`, `mptrac_update_device`, and
`mptrac_update_host`. The functions have the prefix `mptrac` to
indicate that they are interfaces to the original MPTRAC functions.

### Checking order and array sizes

It is crucial that the order and array sizes in the Fortran interface
match those in the original C structure. The structures in MPTRAC
consist of extensive variable lists (up to O100 for the control
parameter list) and data types that may also depend on self-defined
structures (struct of struct). A shell script `find-vars.sh` (in
`tests/wrapper_test/`) can be used to check differences in the
structures and variable dimensions. This script is automatically
executed when the wrapper test is performed.

## Example

An example of using the Fortran wrapper for trajectory calculations
can be found in `tests/wrapper_test/run.sh`.

This file executes the trajectory calculation first with the original C
code `trac.c`. The output is stored in the directory `data/c/`.
Afterwards, the Fortran code `trac_fortran.f90` is called. This output
is stored in the directory `data/fortran/`. Both directories are
compared to make sure that the Fortran version produces the same
results as the C version.
