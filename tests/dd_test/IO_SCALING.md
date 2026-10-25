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
36 simulations. `io_scaling.sh` builds an optional NetCDF interposer
(`netcdf_read_trace.c`) and enables it only for benchmark processes. It
records per-rank core and periodic two-column halo hyperslab times,
requested bytes, shapes, `nc_close()` time, and peak host RSS without
changing the reader. One-second `nvidia-smi` samples record observed GPU
memory use. `io_scaling_report.py` writes `io_summary.csv`, an English
`report.md`, and plots (PNG with Matplotlib, SVG otherwise). Missing or failed
runs remain visible in the CSV. Requested read throughput uses requested
hyperslab bytes, not measured file-system traffic.
