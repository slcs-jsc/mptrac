# DD meteorological I/O scaling at 137 levels

`tests/dd_test/io_scaling.sh` extends the existing GPU-DD benchmark across
four synthetic NetCDF grids: 360 × 181 × 137, 720 × 361 × 137,
1440 × 721 × 137, and 2880 × 1441 × 137. All use `wind.c` with 137
equally spaced synthetic altitude levels from 0 to 60 km. These levels
are intended to increase data volume; they are not an ERA5 or ERA6
vertical coordinate. The default matrix uses four ranks and four GPUs on
one JUWELS Booster node, layouts `1x4`, `2x2`, and `4x1`, approximately
10 million identical global initial particle positions, exchange mode 1,
three repeats, and a six-hour run without diffusion. It reuses seven
hourly wind files across layouts within each resolution and leaves the
standard DD regression-test grid unchanged.

The script performs a dry run by default. It writes `plan.csv` with the
compiled dimensions and conservative memory and disk estimates, plus the
existing benchmark matrix for each resolution. Before execution, inspect
the file-system quota and the Slurm allocation as well as these estimates.
The finest wind generator needs about 41.7 GiB for its global `met_t`
and about 2.2 GiB for its largest write buffer on the host. The DD
executable allocates two *local* halo-extended `met_t` time slices per
rank. Its largest two-slice allocation is about 21.4 GiB per GPU;
the preflight reserves another 4 GiB for particle state, DD buffers,
and accelerator runtime. Seven finest-grid wind files request about
166 GiB before NetCDF metadata and file-system effects. The script
checks free space and limits against `DD_IO_GPU_MEMORY_GIB`,
`DD_IO_HOST_MEMORY_GIB`, and `DD_IO_AVAILABLE_DISK_GIB`; project quotas
must be checked separately.

Within an existing one-node, four-GPU Slurm allocation with the normal
MPTRAC NVHPC, MPI, GSL, and netCDF modules loaded, use:

```bash
export DD_IO_WORK_DIR=/path/to/scratch/dd_io_scaling_137
export DD_COMPILER=nvc DD_MPICC=mpicc DD_USE_GPU=1
./tests/dd_test/io_scaling.sh              # dry run and resource plan
DD_IO_RUN=1 ./tests/dd_test/io_scaling.sh  # run after checking the plan
```

Set `DD_IO_RESOLUTIONS=1` to preflight one resolution first; use a
separate work directory for exploratory runs. The full matrix comprises
36 simulations. By default, it compiles `trac` with `-DMET_IO_PROFILE`.
This optional DD instrumentation times `nc_open_par()`, the explicitly marked
core and periodic-halo float reads, zero-element collective reads, and
`nc_close()` without changing the NetCDF access mode or the existing MPTRAC
timers. It prints one `MET_IO` record per rank and wind file at `LOGLEV=2`.
Set `DD_IO_PROFILE=0` to use the older LD_PRELOAD interposer
(`netcdf_read_trace.c`) instead. The report parser accepts both formats;
archived traces do not contain open or MPI-redistribution timings. The
current interposer records open time but still infers halo reads from
the last hyperslab dimension.
One-second `nvidia-smi` samples record observed GPU memory use, and
`getrusage()` records peak host RSS. `io_scaling_report.py` writes
`io_summary.csv`, an English `report.md`, and plots (PNG with Matplotlib,
SVG otherwise). Missing or failed
runs remain visible in the CSV. Requested read throughput uses requested
hyperslab bytes, not measured file-system traffic.

## Timing semantics

`SELECT_TIMER()` ends the previously selected timer and starts a new one.
`TIMER_GROUP_INPUT` includes particle input and meteorological input across
the simulation. `TIMER_READ_MET_NC_LEVELS` starts before 3D reads and is not
ended until a later timer selection, so it can include `nc_close()` and other
work after the final 3D read. The `MET_IO` measurements are independent local
wall-clock intervals directly around individual NetCDF calls. Their sums are
rank-local accumulated call durations, not a synchronized elapsed phase.

`MET_NC_READ` is the sum of core, periodic-halo, and zero-element NetCDF
float-read call durations. The core and halo categories use the actual DD
reader call sites; they are not inferred from a slab width. A rank without
periodic halo data has zero halo calls and bytes, although the current 2D
reader can still make zero-element collective calls (`empty_calls`). The
requested byte count is the product of the hyperslab dimensions times
`sizeof(float)`; it does not measure HDF5 cache misses or filesystem traffic.
Coordinate reads, NetCDF metadata calls, `nc_var_par_access()`, buffer
allocation, and copying into meteorological arrays are outside these call
intervals. Open and close intervals include only their corresponding NetCDF
calls, not the profiler's logging or `getrusage()`.
They may include waiting inside NetCDF/HDF5; the profiler adds no MPI
barrier and does not separate storage time from library synchronization.

The reusable `MET_MPI_DISTRIBUTE` and `MET_MPI_PACK` hooks in
`src/mptrac.c` are reserved for a future reader. They use the same local
clock and are independent of `SELECT_TIMER()`; both counters remain zero for
the current reader. Time only the MPI operation (including time spent waiting
inside it) in the first category and CPU packing/unpacking in the second.
Place these intervals before the per-file `MET_IO` report and avoid
overlapping intervals within one category.

`io_rank_summary.csv` sums each metric over the seven default wind files
(or `DD_SIMULATION_HOURS + 1` for shorter runs) separately for each rank.
`io_summary.csv` reports the maximum and arithmetic mean of
those rank-local sums for each run. Rank-labeled MPTRAC total, input, and DD
timer means are also retained when available. The Markdown tables take the
median of the three run-level values for each grid/layout. Maxima of
different metrics can occur on different ranks and must not be added to
reconstruct `TIMER_GROUP_INPUT` or total simulation time. Older trace files
are still read; their unavailable metrics are left blank rather than estimated.
