# MPTRAC regression tests

The tests in this directory exercise MPTRAC's command-line tools, data
formats, meteorological processing, particle transport, and selected physical
processes. Most suites generate results in a local `data/` directory and
compare them with checked-in references in `data.ref/`. Other suites check
diagnostics, round trips, hashes, or agreement between two executions.

These are regression and consistency checks for specific configurations.
Passing them does not by itself validate the scientific accuracy of MPTRAC
for other applications, meteorological data sets, or computing platforms.

## Test catalog

| Suite | Main purpose | Important checks and inputs | Default `make check` |
| --- | --- | --- | --- |
| [`atm_test`](atm_test/) | Atmospheric-particle initialization, conversion, selection, statistics, distance measures, and gridding | ASCII/binary/netCDF round trips; statistics and distance output compared with `data.ref/`. Some compiler-sensitive distance moments are intentionally disabled in `run.sh`. | Yes |
| [`cli_test`](cli_test/) | Consistent command-line interface across MPTRAC executables | No-argument errors and `-h`/`--help` behavior for the executables listed by `src/Makefile`; checks generated diagnostics rather than numerical references. | Yes |
| [`cms_test`](cms_test/) | Conversion and downsampling for Cmultiscale interoperability | NetCDF conversion and extracted map results; requires the external Cmultiscale setup referenced by `run.sh`. | No |
| [`compress_test`](compress_test/) | Binary and compressed meteorological-data conversion | Round trips through the compressed representation and SHA-256 checks against checked-in references. | Yes |
| [`coord_test`](coord_test/) | Transport using projected coordinates | Synthetic UTM meteorological input; time-resolved particle positions compared with `data.ref/`. | Yes |
| [`dd_test`](dd_test/) | MPI domain decomposition and particle exchange | Synthetic wind fields, subdomain migration, CPU/GPU comparisons, and optional performance benchmarks. See the [dedicated guide](dd_test/README.md). | No |
| [`depo_test`](depo_test/) | Aerosol and gas deposition, wet scavenging, radioactive decay, and ground inventories | Deposition and activity budgets, cloud-regime cases, and coupled radioactive processes compared with `data.ref/`. | Yes |
| [`dt_test`](dt_test/) | Model time steps, trajectory output intervals, and parameter validation | Short pressure-level trajectories, derived outputs, and rejection of invalid output intervals. | Yes |
| [`gpu_test`](gpu_test/) | GPU-enabled trajectory calculations and output consistency | Model-level and pressure-level transport; particle, grid, ensemble, profile, sample, station, budget, CSI, and VTK references. Requires an appropriate GPU build/environment. | No |
| [`grib_test`](grib_test/) | Reading GRIB meteorological data | Map and zonal-mean extraction from GRIB test inputs; requires the GRIB data and support expected by `run.sh`. | No |
| [`interoper_test`](interoper_test/) | CLaMS interoperability and selected transport/output edge cases | CLaMS position files, diabatic transport, missing-coordinate diagnostics, task isolation, ensemble IDs, sample VMR, and end-of-run outputs. | Yes |
| [`kpp_test`](kpp_test/) | Chemical integration through KPP | Trajectory calculations with the KPP configuration and comparison with particle references. Requires the corresponding KPP-enabled build. | No |
| [`met_test`](met_test/) | Meteorological analysis and sampling tools | Maps, profiles, zonal means, spectra, lapse rates, subgrid variability, and sampling using ERA-Interim and ERA5 input. | Yes |
| [`tools_test`](tools_test/) | Small numerical and physical utilities | Time conversion, sedimentation, and NAT-related calculations over parameter grids. | Yes |
| [`trac_test`](trac_test/) | Main serial/shared-memory transport regression | Model-level and pressure-level trajectories and the associated particle, gridded, ensemble, diagnostic, and VTK outputs. | Yes |
| [`tropo_test`](tropo_test/) | Tropopause diagnosis and derived products | Cold-point, dynamical, and first/second WMO tropopause definitions; sampling, climatology, and zonal means. | Yes |
| [`wrapper_test`](wrapper_test/) | C/Fortran interface consistency | Checks structure-variable mappings and compares the Fortran wrapper's transport output with the C implementation. Requires the wrapper build. | No |

The shared [`tests/data/`](data/) directory contains meteorological input
files, including ERA-Interim, ERA5 model-level, and projected ERA5 examples,
used by multiple suites. Individual test directories may additionally
contain configuration files, generated test inputs, or dedicated
`data.ref/` files.

## Running the tests

Build MPTRAC and run the default regression suite from the repository root:

```bash
cd src
make check
```

The default target runs `atm_test`, `cli_test`, `compress_test`,
`coord_test`, `depo_test`, `dt_test`, `interoper_test`, `met_test`,
`tools_test`, `trac_test`, and `tropo_test`.

To run a single test with its usual build dependencies, use its Makefile
target from `src/`:

```bash
cd src
make met_test
```

Alternatively, run a suite directly from its own directory after building
the required executables:

```bash
cd tests/met_test
./run.sh
```

The scripts generally assume their test directory as the working directory
and refer to executables under `../../src/` and libraries under
`../../libs/build/lib/`. Optional suites can be invoked using their
Makefile targets, such as `make grib_test`, `make kpp_test`,
`make gpu_test`, `make cms_test`, or `make wrapper_test`, provided their
additional dependencies and build configurations are available.

The domain-decomposition suite has additional launch modes, configuration
variables, and benchmark scripts. Consult [`dd_test/README.md`](dd_test/README.md)
before running it; its MPI/GPU and HPC workflows are not equivalent to
the ordinary single-process regression tests.

## References and reproducibility

Most suites compare generated text (`.tab`) or VTK files with the versioned
files in their own `data.ref/` directories, often using exact `diff`
comparisons. The compression and Cmultiscale tests additionally use
stored SHA-256 checksums. CLI and other focused assertions inspect program
exit status, diagnostics, or properties of generated output instead of
requiring an external baseline.

The individual `run.sh` files are the authoritative description of the
commands, input parameters, generated files, and comparison procedures.
Checked-in references preserve expected results for regression detection;
the repository does not necessarily document the original generation
environment or scientific provenance of every baseline.

When an intentional code change affects a reference result, inspect the
difference carefully, verify the physical or format-related reason, and
update only the affected reference files. Review baseline changes together
with the corresponding code changes. Exact comparisons can be sensitive to
compiler, library, floating-point, and platform differences; specialized
GPU and HPC tests should therefore be interpreted in their intended
environment.

## Regression testing and scientific validation

The test suite is part of MPTRAC's software quality assurance. Its limited
input cases and reference comparisons are not a substitute for independent
verification of numerical methods, comparisons with observations, or
application-specific scientific validation.
