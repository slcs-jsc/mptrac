# High-Level interface

## Overview

The MPTRAC Lagrangian transport model provides a high-level interface for simulating atmospheric transport processes.

This high-level interface, implemented through the `mptrac_*` functions in [`mptrac.c`](https://github.com/slcs-jsc/mptrac/blob/master/src/mptrac.c) and [`mptrac.h`](https://github.com/slcs-jsc/mptrac/blob/master/src/mptrac.h), facilitates memory management, initialization, data input/output, and simulation execution.

This document describes the core functions of the interface and demonstrates their usage with an example workflow.

## Function categories

### 1. Memory Management

- **`void mptrac_alloc(...)`**: Allocates memory for control, cache, climatology, meteorology, and atmospheric data structures on the CPU and GPU.

- **`void mptrac_free(...)`**: Frees memory allocated for these structures.

### 2. Initialization and memory updates

- **`void mptrac_init(...)`**: Initializes the MPTRAC model with control, cache, climatology, and atmospheric data. It adjusts the time range of the simulation and initializes the random number generator.

- **`void mptrac_update_device(...)`**: Updates device memory with the latest data from CPU memory.

- **`void mptrac_update_host(...)`**: Updates host memory from GPU memory.

### 3. Data input

- **`void mptrac_read_ctl(...)`**: Reads control parameters from a file or the command line.

- **`void mptrac_read_clim(...)`**: Reads various climatology data from data files.

- **`int mptrac_read_met(...)`**: Reads meteorological data from a file.

- **`int mptrac_read_atm(...)`**: Reads air parcel data from a file.

### 4. Data output

- **`void mptrac_write_output(...)`**: Writes simulation results in various output types.

- **`void mptrac_write_met(...)`**: Writes meteorological data to a file.

- **`void mptrac_write_atm(...)`**: Writes air parcel data to a file.

### 5. Simulation execution

- **`void mptrac_run_timestep(...)`**: Executes a single timestep of the Lagrangian transport model. Simulates various processes such as advection, diffusion, convection, chemistry, etc.

### 6. Meteorological data handling

- **`void mptrac_get_met(...)`**: Retrieves meteorological data for a specific time.

## Example workflow

Below is an example of how to use the high-level MPTRAC interface in a typical simulation. Please see [`trac.c`](https://github.com/slcs-jsc/mptrac/blob/master/src/trac.c) for the full code.

### 1. Allocate memory

```
mptrac_alloc(&ctl, &cache, &clim, &met0, &met1, &atm);
```

### 2. Initialize the model

```
mptrac_read_ctl(filename, argc, argv, ctl);
mptrac_read_clim(ctl, clim);
mptrac_read_atm(filename, ctl, atm);
mptrac_init(ctl, cache, clim, atm, ntask);
```

### 3. Run the simulation

Within the time loop of the model, repeatedly call:

```
mptrac_get_met(ctl, clim, t, &met0, &met1);
mptrac_run_timestep(ctl, cache, clim, &met0, &met1, atm, t);
mptrac_write_output(dirname, ctl, met0, met1, atm, t);
```

### 4. Free memory

```
mptrac_free(ctl, cache, clim, met0, met1, atm);
```

## Notes

- Error handling is essential when reading input data or allocating memory.

- The example workflow shown here properly handles GPU offloading and data transfers between CPU and GPU memory. In other application, ensure proper synchronization between host and device memory when using GPUs.

## References

- [`mptrac.h`](https://github.com/slcs-jsc/mptrac/blob/master/src/mptrac.h): Header file with function declarations.

- [`trac.c`](https://github.com/slcs-jsc/mptrac/blob/master/src/trac.c): Example implementation of the MPTRAC interface.
