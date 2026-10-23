# Parallelization

## Introduction

MPTRAC combines MPI, OpenMP, and OpenACC. MPI supports two ways to
distribute work across ranks:

- **Independent simulations:** Each rank runs a separate simulation
  selected from a list of working directories. The ranks do not exchange
  particles. This is the established MPI ensemble mode.
- **Spatial domain decomposition (DD):** The ranks share one simulation.
  Each rank owns a geographical subdomain and exchanges particles with
  other ranks as they move between subdomains.

OpenMP distributes particle calculations among CPU threads within each
rank. OpenACC can run those calculations on a GPU. Either MPI strategy
can be combined with OpenMP and OpenACC, so multiple ranks may use the
CPUs or GPUs of one or more compute nodes.

## Parallelization strategy of MPTRAC

### Independent simulations

Build `trac` with `MPI=1` and use `DD=0` in the control files. The
directory list passed to `trac` contains one working directory per
simulation. MPI ranks take entries from this list in round-robin order;
a rank may process several simulations sequentially. Each simulation
has its own initial particles, control parameters, and meteorological
data. No particle exchange takes place between these simulations.

If the executable was also built with `DD=1`, set the runtime control
parameter `DD=0` explicitly for independent simulations. Otherwise,
DD is enabled by default when that executable starts with more than one
MPI rank.

### Spatial domain decomposition

Build `trac` with `MPI=1 DD=1` and set `DD=1` in the control files.
Divide one simulation into zonal and meridional subdomains using
`DD_SUBDOMAINS_ZONAL` and `DD_SUBDOMAINS_MERIDIONAL`. Start exactly one
MPI rank per subdomain, with one working-directory entry per rank in
the directory list. Each rank reads its local meteorological region
with halos and owns the particles currently in that region. Particles
crossing a subdomain boundary are sent to their new owner by MPI.
DD requires at least two subdomains and a latitude/longitude
meteorological grid; radioactive deposition is not supported.

`DD_EXCHANGE_MODE=0` (the default) uses the sort/push exchange path;
`DD_SORT_DT` controls periodic sorting and defaults to 1800 s.
`DD_EXCHANGE_MODE=1` classifies particles, packs outgoing particles, and
compacts retained particles without sorting. It also removes particles
with invalid subdomain indices. In mode 1, neither `DD_SORT_DT` nor the
general `SORT_DT` triggers sorting. Both modes use the same subdomain
assignment and MPI particle communication.

With OpenACC, mode 1 performs classification, packing, and local
compaction on the GPU. The active particle and cache arrays remain on
the device; only outgoing particle records are copied to the host for
the current MPI exchange, and received records are copied to the GPU.
The `iso_var` quantity and the three `uvwp` cache components follow particles
between ranks. `uvwp` is included in the MPI particle datatype and restored
to the cache after each exchange. Particle files remain unchanged; `uvwp`
is an internal state variable.

The Fortran wrapper supports both DD exchange modes when built with
`MPI=1 DD=1` and a matching MPI Fortran compiler (`MPIFC`). The DD branch
of `wrapper_test` compares `trac_fortran` with `trac` using mode 1. Build
from `src/` with the DD test dimensions:

```bash
make clean
make DD=1 MPI=1 MPIFC=mpifort \
  DEFINES='-DNP=32 -DEX_GLOB=364 -DEY_GLOB=186 -DEP=90 -DEX=130 -DEY=70 -DNQ=5' \
  wrapper_test
```

### OpenMP and OpenACC within each rank

OpenMP is enabled in MPTRAC builds and parallelizes particle
calculations within an MPI rank. Its threads share that rank's
meteorological data. Different MPI ranks have separate data, even
when they run on the same node. The OpenMP scalability of MPTRAC is
discussed by [Rößler et al. (2018)](https://doi.org/10.5194/gmd-11-575-2018).

Build with `GPU=1` to offload suitable calculations to GPUs through
OpenACC. For example, four ranks on a four-GPU node can run four
independent simulations or four subdomains of one DD simulation. The
choice of ranks, OpenMP threads, and GPUs depends on the problem size
and the target system.

## Supercomputers in Jülich

The JUWELS Cluster module has 2271 standard compute nodes, each with
48 CPU cores and 96 GiB of memory. Simultaneous multithreading provides
up to 96 hardware threads per node. The JUWELS Booster module has 936
compute nodes, each with 48 CPU cores, 512 GiB of memory, and four
NVIDIA A100 GPUs. See the [JUWELS hardware
configuration](https://apps.fz-juelich.de/jsc/hps/juwels/configuration.html)
for details.

Examples of job scripts and resource requests are available in the
[JUWELS batch-system documentation](https://apps.fz-juelich.de/jsc/hps/juwels/batchsystem.html).

Alternatively, MPTRAC has also been used on the [JURECA-DC
system](https://apps.fz-juelich.de/jsc/hps/jureca/quickintro.html).

## Further Reading

- [MPI Forum](https://www.mpi-forum.org/)
- [OpenMP](https://www.openmp.org/)
- [OpenACC](https://www.openacc.org/)
