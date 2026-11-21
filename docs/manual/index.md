# Welcome to MPTRAC!

Massive-Parallel Trajectory Calculations (MPTRAC) is an open-source
Lagrangian particle transport and dispersion model. Driven by
meteorological reanalysis or forecast data, it tracks atmospheric
air parcels and represents transport, mixing, and selected physical
and chemical processes. MPTRAC is primarily designed for applications
in the free troposphere and stratosphere, including volcanic plumes,
wildfire emissions, and long-range transport.

![Lagrangian transport simulation of convective transport](img/convection_lowres.jpg)

## Key features

- **Trajectory simulations:** Forward and backward transport driven by
  meteorological wind and vertical velocity fields.
- **Dispersion and processes:** Diffusion, mixing, convection,
  sedimentation, chemistry, decay, and deposition.
- **Meteorological tools:** Preprocessing and diagnostics for
  reanalysis and forecast data.
- **Analysis and output:** Particle, grid, ensemble, profile, and
  sampling outputs for subsequent analysis and visualization.
- **Parallel computing:** Hybrid MPI, OpenMP, and OpenACC execution
  on CPUs and GPUs, including spatial domain decomposition.
- **Open source and integration:** C implementation under the GNU GPL v3,
  with a high-level API and Fortran wrapper.

![Geophysical modules and main software components of MPTRAC](img/mptrac_modules.png)

## Getting started

Use [Installation](installation.md) to set up MPTRAC and
[Quickstart](quickstart.md) to run a first trajectory simulation.
The [MPTRAC Web Runner](https://mptrac.jsc.fz-juelich.de)
provides a browser-based introduction without local installation.

## Explore the documentation

- **Scientific background:** [Model physics](model-physics.md) and
  [Diabatic transport](diabatic-transport.md).
- **Using the model:** [Applications and workflows](applications.md),
  [Control parameters](control-parameters.md), and
  [Meteorological input data](meteorological-input-data.md).
- **HPC and development:** [Parallelization](parallelization.md),
  [Code development](code-development.md), and
  [Testing and code quality](testing.md).

## Citation and contributions

Please consult [References](references.md) for scientific publications
and the [contribution guidelines](https://github.com/slcs-jsc/mptrac/blob/master/CONTRIBUTING.md)
if you would like to report an issue or contribute.

## Contact

Dr. Lars Hoffmann, <l.hoffmann@fz-juelich.de>  
Jülich Supercomputing Centre, Forschungszentrum Jülich, Germany
