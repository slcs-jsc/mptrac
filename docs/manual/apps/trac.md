# trac

`trac` is the main application for forward and backward trajectory
calculations. It is called with a directory list, a control file name,
and an atmospheric input file name:

```bash
./trac <dirlist> <ctl> <atm_in> [KEY VALUE ...]
```

Required arguments:

- `dirlist`: text file containing the work directories to be processed.
- `ctl`: control file name relative to each directory from `dirlist`.
- `atm_in`: atmospheric input file name relative to each directory from `dirlist`.

Optional trailing arguments can be used to override control parameters
on the command line as `KEY VALUE` pairs.

Example from the repository:

```bash
./trac data/dirlist trac.ctl atm_split.tab
```

Example with command-line overrides:

```bash
./trac data/dirlist trac.ctl atm_split.tab \
  ATM_BASENAME atm_diff \
  GRID_BASENAME grid_diff
```

The effective configuration is defined by the `ctl_t` control
parameters. A full field-by-field reference is available in the
[Doxygen manual](https://slcs-jsc.github.io/mptrac/doxygen/structctl__t.html).
Practical examples, file syntax, and precedence rules are documented in
the [control parameters](../control-parameters.md) section of the user
manual.
