# MPTRAC third-party libraries

This directory contains source archives for third-party libraries used by
MPTRAC and the [`build.sh`](build.sh) script for building them locally.
The script installs selected components under `libs/build/`; it does not
install them into the system directories.

## Bundled sources

| Archive | Version | Purpose |
| --- | --- | --- |
| `gsl-2.7.1.tar.bz2` | 2.7.1 | GNU Scientific Library, used by MPTRAC's numerical routines |
| `zlib-1.3.1.tar.bz2` | 1.3.1 | Compression support for the HDF5/netCDF stack |
| `szip-2.1.1.tar.bz2` | 2.1.1 | Additional HDF5 compression support |
| `hdf5-1.14.4-3.tar.bz2` | 1.14.4-3 | Hierarchical scientific data storage |
| `netcdf-c-4.9.2.tar.bz2` | 4.9.2 | netCDF input/output |
| `zfp-1.0.1.tar.bz2` | 1.0.1 | Floating-point data compression |
| `zstd-1.5.7.tar.bz2` | 1.5.7 | Zstandard compression |
| `lz4-1.9.4.tar.bz2` | 1.9.4 | LZ4 compression |
| `SZ3-3.2.1.tar.bz2` | 3.2.1 | Lossy scientific data compression |
| `eccodes-2.38.3.tar.bz2` | 2.38.3 | GRIB meteorological-data support |
| `KPP.tar.bz2` | Not encoded in archive name | Kinetic PreProcessor for chemical integration |
| `thrustsort-1.4.tar.bz2` | 1.4 | Optional CPU/GPU Thrust sorting wrappers |

The archives are checked into the repository. Their availability makes
the corresponding source files accessible without downloading them during
the build; a working compiler and any build-time tools are still required.
See `build.sh` for the exact configuration commands and build options.

## Building

Run the script **from this directory** so that its relative paths resolve
correctly:

```bash
cd libs
./build.sh -m
```

`-m` selects the mandatory library stack: GSL, zlib, SZIP, HDF5, and
netCDF-C. The script builds these in dependency order. A larger build is:

```bash
./build.sh -a
```

`-a` removes the existing `libs/build/` directory and builds all supported
components **except the optional Thrust wrappers**. This may also build KPP
chemistry and libraries not required by a particular MPTRAC configuration.

To select individual components, or limit compilation parallelism, use
options such as:

```bash
./build.sh -g -p 8
./build.sh -tu
./build.sh -h
```

The first command builds GSL with eight make jobs. The second builds both
CPU and GPU Thrust wrappers, requiring the corresponding Thrust headers
and, for the GPU variant, `nvcc`. The last command displays the complete
option list. Additional selections include `-e` for ecCodes, `-k` for
KPP, `-f` for ZFP, `-s` for Zstandard, `-x` for LZ4, and `-j`
for SZ3. The `-l` option enables an MPI-based parallel HDF5/netCDF build
and requires `mpicc`.

The script builds from the bundled source archives. For components that
depend on each other, build their prerequisites first or select the
mandatory stack. Compiler, MPI, and GPU toolchain compatibility depend
on the selected options and target platform.

## Build output

The local installation prefix is `libs/build/`, which contains directories
such as `bin/`, `include/`, `lib/`, and `src/`. The script adds
`libs/build/lib` to `LD_LIBRARY_PATH` for its own execution; users may
also need to configure library search paths when running the resulting
MPTRAC executables.

`./build.sh -c` removes the build directory without compiling libraries.
`./build.sh -a` also cleans it first. Both discard previously built
artifacts in `libs/build/`.

To build MPTRAC itself after installing the required dependencies, follow
the repository's main build documentation and [`src/Makefile`](../src/Makefile).
