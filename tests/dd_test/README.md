# DD test and benchmark guide

The scripts in this directory run DD regression tests and optional performance
measurements. Commands below assume the repository root as the working directory.

## DD performance benchmark

`tests/dd_test/benchmark.sh` runs a matrix of the sort/push exchange path
(`DD_EXCHANGE_MODE=0`) and the sort-free compaction path
(`DD_EXCHANGE_MODE=1`). It uses the existing DD setup in `run.sh` with
`DD_BENCHMARK=1` in the test script. This setting selects benchmark input and
output; MPTRAC diagnostics use the regular compile-time `LOGLEV`. The normal
regression mode and its reference files are unchanged. By default the
benchmark uses a 2 × 2 DD layout, approximately
100,000, 1 million, and 10 million particles, both exchange modes, a
six-hour smooth-wind simulation, and no diffusion. Set `DD_BENCH_COUNTS`
explicitly to include 100 million particles; that size is never run by
default. `DD_BENCH_MODES` and `DD_BENCH_REPEATS` select modes and repeats;
their execution order alternates on successive repeats.

`DD_WIND_ALPHA` controls the axis of the analytic horizontal wind field in
degrees. At `0`, the field is zonal (`u = U cos(latitude)`, `v = 0`). At
`90`, both components vary with position; it is predominantly meridional
in some regions, not a uniform north-south wind. The amplitude parameter
`WIND_U0`/`WIND_U1` is unchanged when rotating the field, although the
velocity sampled at a fixed particle position can change.

For a 1° wind grid, `WIND_PERT_K=18` to `36` gives a longitude wavelength
of 20° to 10°. For example, `WIND_PERT_U=3`, `WIND_PERT_V=3`, and
`WIND_PERT_W=0.03` add spatial variability for diffusion regression runs.

`DD_BENCH_LAYOUTS=1x4,2x2,4x1` compares layouts with the same MPI rank
count. For each particle target, the script chooses one global longitude ×
latitude particle grid divisible by every requested layout. The ordinary
`atm_init` path rounds each layout's local grid spacing to six decimal
places, so positions can differ slightly between layouts. Set
`DD_BENCH_CANONICAL_INIT=1` for an exact global cell-center grid when
comparing layouts; this leaves particle IDs and the total count unchanged.
Initial files are reused only within the same
layout. Each size, layout, repeat, and exchange mode has its own directory.
Meteorological input can be shared across layouts if the wind settings are
identical. Set both `DD_PARTICLES_LON_STRIDE` and
`DD_PARTICLES_LAT_STRIDE` only for a single count and a single layout.

Benchmark mode uses MPTRAC binary particle files (`ATM_TYPE=1`,
`ATM_TYPE_OUT=1`). `atm_init` writes `init.bin`; each rank writes binary
snapshots at the beginning and end of the default six-hour run. The output
interval is `DD_ATM_DT_OUT=86400` seconds by default and must include the
initial time and be divisible by the 600-second model step. Longer runs or
shorter output intervals can create additional snapshots. The script
checks expected disk use before compiling or running and limits it to
`DD_BENCH_MAX_DISK_GIB=64` by default. Initial files are hard-linked
between modes and repeats for the same particle count and layout. Meteorological
input is generated once and cached in the benchmark work directory, or
can be supplied with `DD_WIND_DIR`.
For a direction with no boundary crossings, such as purely zonal flow in
`1x4`, set `DD_BENCH_ALLOW_ZERO_MIGRATION=1` to retain the valid zero-send
case in the summary. The default still requires nonzero migration.

Particle counts come from the common global grid, divided among the MPI
ranks. `DD_NP_MAX` is the compiled particle capacity **per
rank**. Its benchmark default is the initial per-rank count plus a
configurable `DD_NP_MARGIN_PERCENT=25` margin and one extra slot. Increase
the margin if particles concentrate on particular ranks. The benchmark
checks that NP exceeds the initial count and that IDs fit in signed 32-bit
integers. At `NQ=5`, the atmospheric arrays, particle cache, and DD helper
arrays alone use approximately `152 × DD_NP_MAX` bytes per rank. Host and
GPU copies, meteorology, MPI buffers, and allocator overhead add to this;
the CSV memory estimate is therefore a lower bound. The preflight limit
`DD_BENCH_MAX_PARTICLE_MEMORY_GIB=24` applies to this lower bound and can be
overridden. Each binary snapshot uses about 72 bytes per particle, so 100 million particles require about
7.2 GB for each snapshot and another 7.2 GB for the initial files.

Run inside a JUWELS Booster Slurm allocation. The benchmark never submits
another job. `DD_BENCH_NODES`, `DD_BENCH_RANKS`,
`DD_BENCH_RANKS_PER_NODE`, and `DD_BENCH_GPUS_PER_RANK` specify the launch
topology; GPU benchmarks currently use one GPU per rank. The layout product
must equal the rank count, and ranks must be evenly distributed across
nodes. The launcher passes explicit `--nodes`, `--ntasks`,
`--ntasks-per-node`, and `--gpus-per-task` options to `srun`. MPTRAC records
the actual shared-memory node of every rank in `trac.log`; the report checks
the requested placement and writes `rank_mapping.csv`, including the
visible GPU assignment. If the actual placement differs, the case fails
validation. Slurm benchmark logs label each line with its rank, allowing
per-rank timer values to be inspected without ambiguity.
On JUWELS Booster, `srun` can expose the full allocated GPU list to every
rank even with GPU binding enabled. The benchmark launcher selects one
visible GPU by `SLURM_LOCALID` before starting `trac`, and the recorded
`CUDA_VISIBLE_DEVICES` values are checked for overlap within each node.

For a four-GPU, single-node allocation, replace the account and wall time
with suitable values:

```bash
#SBATCH --partition=booster
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --account=<account>
#SBATCH --time=02:00:00
module purge
module load nvidia-compilers ParaStationMPI CMake
```

Set a shared scratch work directory and launch from the repository root:

```bash
export DD_BENCH_WORK_DIR=/path/to/scratch/dd_benchmark
export DD_USE_GPU=1 DD_COMPILER=nvc DD_MPICC=mpicc
export DD_BENCH_COUNTS=100000,1000000,10000000
export DD_BENCH_MODES=0,1 DD_BENCH_REPEATS=3
./tests/dd_test/benchmark.sh
```

Run the three comparison phases in separate allocations and separate work
directories. Phase A compares domain geometry on one node:

```bash
export DD_BENCH_LAYOUTS=1x4,2x2,4x1 DD_BENCH_COUNTS=1000000,10000000
export DD_BENCH_MODES=0,1 DD_BENCH_REPEATS=3
export DD_BENCH_NODES=1 DD_BENCH_RANKS=4 DD_BENCH_RANKS_PER_NODE=4
./tests/dd_test/benchmark.sh
```

For phase B, run `2x2` with 10 million particles once in the allocation
above, then in a **two-node** allocation with four total tasks, two tasks
per node, and two GPUs on each node. Set a different work directory and:

```bash
#SBATCH --partition=booster
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
```

Within that allocation:

```bash
export DD_BENCH_LAYOUTS=2x2 DD_BENCH_COUNTS=10000000
export DD_BENCH_NODES=2 DD_BENCH_RANKS=4 DD_BENCH_RANKS_PER_NODE=2
export DD_BENCH_MODES=0,1 DD_BENCH_REPEATS=3
./tests/dd_test/benchmark.sh
```

For phase C, use 4, 8, or 16 tasks and GPUs in matching allocations.
Suitable layouts include `2x2` for four ranks, `2x4,4x2` for eight, and
`4x4,2x8,8x2` for sixteen. Set `DD_BENCH_RANKS` and
`DD_BENCH_LAYOUTS` accordingly; keep the total particle target fixed for
strong-scaling comparisons. The default layout is near square when an
explicit rank count is given without a layout.

For a single optional large run, set `DD_BENCH_COUNTS=100000000`,
`DD_BENCH_REPEATS=1`, and provide sufficient GPU memory, wall time, and
scratch space. Set `DD_BENCH_DRY_RUN=1` first to print capacity and disk
estimates without compiling or running. Start with a small matrix. The same
script accepts
`DD_DOMAINS_LON`, `DD_DOMAINS_LAT`, `DD_SIMULATION_HOURS`, `DD_SORT_DT`,
`DD_EXCHANGE_MODE` through `DD_BENCH_MODES`, `DD_DIFFUSION`,
`DD_TURB_MESOX`, `DD_TURB_MESOZ`, and the `DD_WIND_PERT_*` parameters.
`DD_THRUST=1` selects GPU Thrust sorting in mode 0; set it to 0 to compare
the GSL sorting backend. NVIDIA benchmark builds use `-Mlarge_arrays`
through `DD_BENCH_OPT` by default to support large NP values. A local CPU
smoke test can use `DD_USE_GPU=0`,
`DD_COMPILER=gcc`, and an existing `DD_WIND_DIR`. If `mpicc` in the
current `PATH` wraps NVIDIA's compiler, set `DD_MPICC` and `DD_MPIRUN`
to matching system MPI executables for the GCC build.

The script stores each run, full logs, and validated binary outputs in a
separate directory and appends one row to `summary.csv`. It checks every
particle ID in each snapshot with a compact bitset. `migrated_particles`
counts IDs whose final rank differs from the rank encoded by their initial
ID; it does not count repeat crossings or particles that return to their
original rank. `sent_particle_events`, `intra_node_particle_events`, and
`inter_node_particle_events` sum the MPI send counts over every exchange,
so repeated crossings count again. Byte columns use the MPI datatype's
transmitted size, which can be smaller than its in-memory extent.
`inter_node_event_fraction` is the fraction of sends crossing nodes. These
statistics use existing host-side send counts and do not copy particle
arrays from the GPU. Times are the maximum
reported across MPI ranks. `launcher_wall_time` separately includes launcher
startup and process shutdown.
The benchmark report requires `LOGLEV >= 2` for rank-to-node mapping and
communication totals. At level 1, MPTRAC reports the DD configuration,
initial local particle count, final `SIZE_NP`, and `SIZE_NP_PEAK`. Level 2
adds node placement, optional visible GPU assignment, and per-rank cumulative
intra-node and inter-node send counts and payload bytes. Level 3 adds one
send-count line per DD exchange. The node topology is determined once during
DD initialization; send statistics reuse the existing MPI counts.
`mpi_time` is MPTRAC's `DD_COMMUNICATE_PARTICLES` phase, which includes
host-side preparation in addition to MPI calls. `gpu_d2h_time` is the
separately timed outgoing transfer in mode 1. A complete GPU transfer time,
pure host-to-GPU time, and pure unpacking time are not isolated by the
existing timers and remain blank or are reported only as combined phases.
`TIMER_GROUP_INPUT` accumulates particle-file and meteorological reads
throughout the simulation. `TIMER_READ_MET_NC_LEVELS` includes the 3D
netCDF reads and parallel file close.
The optional NetCDF call profiler and its rank/repetition aggregation are
documented in [IO_SCALING.md](IO_SCALING.md).
Changing node placement also changes available CPUs and GPU affinity, so
differences in total runtime cannot be attributed solely to network traffic.

## Fortran wrapper DD test

The DD branch of `wrapper_test` compares `trac_fortran` with `trac` using
exchange mode 1. Build from `src/` with the DD test dimensions:

```bash
make clean
make DD=1 MPI=1 MPIFC=mpifort \
  DEFINES='-DNP=32 -DEX_GLOB=364 -DEY_GLOB=186 -DEP=90 -DEX=130 -DEY=70 -DNQ=5' \
  wrapper_test
```
