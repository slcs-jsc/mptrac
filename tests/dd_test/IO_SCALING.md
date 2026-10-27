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
36 simulations. Component read times come from the built-in MPTRAC timers
(`TIMER_MET_*`).
Set `DD_MET_NC_READ_MODE=1` to benchmark the global horizontal slice reader;
the default `0` keeps the existing local hyperslab reader. Use a new
`DD_IO_WORK_DIR` for mode 1, since the script does not overwrite results.
The mode 1 reader allocates one global float32 horizontal field on each
rank assigned a level and reuses it over the file. At 0.25° this field is
about 4 MiB per rank; at 0.125° it is about 16 MiB. It also needs one
largest-destination packed slice and one local slice per rank. Existing
preflight estimates remain conservative for the much larger `met_t` arrays.
For direct field verification, run `bash tests/dd_test/met_reader_compare.sh`
on a CPU node with parallel NetCDF. It generates a small private synthetic
file, encodes global x/y/level coordinates and invalid samples, and compares
the complete local `met_t` structures bytewise for 1×1, 1×4, 2×2, and 4×1
with two halo points. It also checks a duplicate 0°/360° endpoint for 1×1
and 1×4. The 1×1 case exercises the reader geometry directly;
normal particle DD simulations require at least two ranks. A duplicate
0°/360° endpoint is currently accepted only without a zonal split: the
existing zonal DD geometry produces duplicate local longitudes and fails
its regular-spacing check, independent of read mode.
The DD met-input instrumentation uses MPTRAC's standard timers and requires
no extra compile flags. The timer output at `LOGLEV=1` includes
`TIMER_MET_READ` (core, halo, and global float reads), `TIMER_MET_PACK`
(CPU packing), `TIMER_MET_MPI` (MPI redistribution), and `TIMER_MET_CONVERT`
(conversion into the destination buffers), plus `TIMER_GET_MET`, which covers
`nc_open_par()`, grid reading, setup, and `nc_close()`.
Set `DD_IO_PROFILE=0` to additionally build the older LD_PRELOAD interposer
(`netcdf_read_trace.c`). The report parser still accepts `MET_IO` records
from archived logs of older instrumented builds, but new runs no longer emit
them. The
interposer records open time but still infers halo reads from
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

`MET_NC_READ` is the sum of core, periodic-halo, global-slice, and
zero-element NetCDF float-read call durations. The core and halo categories use the actual DD
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

`MET_MPI_DISTRIBUTE` measures blocking point-to-point transfers, including
waiting; `MET_MPI_PACK` measures CPU packing into local contiguous slices.
Both use a local clock independent of `SELECT_TIMER()` and remain zero for
mode 0. `file_elapsed_s` spans the full rank-local file read, preparation,
redistribution, and close. It has no added MPI barrier. Global read calls
are counted only on the assigned reader rank, and requested bytes are the
NetCDF hyperslab size rather than physical storage traffic.

`io_rank_summary.csv` sums each metric over the seven default wind files
(or `DD_SIMULATION_HOURS + 1` for shorter runs) separately for each rank.
`io_summary.csv` reports the maximum and arithmetic mean of
those rank-local sums for each run. Rank-labeled MPTRAC total, input, and DD
timer means are also retained when available. The Markdown tables take the
median of the three run-level values for each grid/layout. Maxima of
different metrics can occur on different ranks and must not be added to
reconstruct `TIMER_GROUP_INPUT` or total simulation time. Older trace files
are still read; their unavailable metrics are left blank rather than estimated.

The report's **DD and I/O timing breakdown** uses the maximum rank-local
sum for each component in a run, then the median across completed repeats.
`io_rank_summary.csv` retains the individual rank values and
`io_summary.csv` also records their mean. For detailed `TIMER_DD_*` entries,
the reporter uses the logged call mean times call count to reduce rounding
from the three-decimal printed total. `DD_GPU_CLASSIFY`, `DD_GPU_COMPACT`,
`DD_GPU_PACK`, `DD_GPU_SEND_D2H`, and `DD_COMMUNICATE_PARTICLES` supply the
particle rows; the communication timer includes host-side preparation.
`TIMER_GROUP_DD` is the separate particle-DD total. `MET_IO file_elapsed_s`
covers NetCDF file reading and MPI redistribution, but not later
meteorological processing. `TIMER_GROUP_INPUT` includes particle input as
well as meteorological input. These independent measurements are shown as
separate rows, never added to form a total. Archived traces without a
measurement display N/A.
