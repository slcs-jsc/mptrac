# Code structure

MPTRAC is a Lagrangian particle transport model implemented in C.
This page is a guide to the main source files and interfaces, rather
than a replacement for the [Doxygen API reference](https://slcs-jsc.github.io/mptrac/doxygen/).

## Core library and data

- [`src/mptrac.h`](https://github.com/slcs-jsc/mptrac/blob/master/src/mptrac.h)
  declares the public API, constants, and data structures.
- [`src/mptrac.c`](https://github.com/slcs-jsc/mptrac/blob/master/src/mptrac.c)
  implements the core library routines.
- `ctl_t` holds simulation control parameters; `atm_t` describes the
  particle state; `met_t` stores meteorological fields; `clim_t`
  contains climatological data; and `cache_t` holds runtime and
  interpolation caches.

## Main workflow

The [`trac.c`](https://github.com/slcs-jsc/mptrac/blob/master/src/trac.c)
application is the main trajectory driver. The
[high-level interface](high-level-interface.md) describes the
`mptrac_*` entry points for allocation, initialization, input,
meteorological access, timestep execution, output, and cleanup.

A typical simulation reads its configuration, climatology, initial
particles, and meteorological fields; advances the particles through
the selected transport and physical processes; and writes the
requested outputs. Optional processes include diffusion, convection,
sedimentation, chemistry, deposition, and mixing.

## Tools and utilities

Standalone applications include `atm_init`, `atm_split`,
`atm_select`, and `atm_stat` for particle data, plus `met_map`,
`met_prof`, and `met_zm` for meteorological diagnostics.
For end-to-end usage, see the [User manual](quickstart.md).

## Parallel execution

MPTRAC supports independent MPI simulations and spatial domain
decomposition, with OpenMP CPU threading and optional OpenACC GPU
offloading. The [Parallelization](parallelization.md) chapter
documents the supported modes and domain-decomposition exchange paths.

## Where to look next

- [Building from source](building.md)
- [Testing and code quality](testing.md)
- [Profiling and benchmarking](profiling.md)
- [High-level interface](high-level-interface.md)
- [Doxygen API reference](https://slcs-jsc.github.io/mptrac/doxygen/)
