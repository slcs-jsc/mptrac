# atm2grid

`atm2grid` converts an atmospheric particle data file to gridded output
data. It reads the output settings from the control file and writes a
grid file using the configured `GRID_BASENAME`.

```bash
./atm2grid <ctl> <atm_in> [KEY VALUE ...]
```

Required arguments:

- `ctl`: control parameter file.
- `atm_in`: atmospheric input file.

Optional trailing arguments can be used to override control parameters
on the command line as `KEY VALUE` pairs.

The control parameter `GRID_BASENAME` must be set. The output filename
is derived from `GRID_BASENAME` and the time stamp read from `atm_in`.
For tabular grid output, the file suffix is `.tab`; for netCDF grid
output, the file suffix is `.nc`.

Example:

```bash
./atm2grid trac.ctl atm_2011_06_05_00_00.tab GRID_BASENAME grid
```

