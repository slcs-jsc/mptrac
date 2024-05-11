# Massive-Parallel Trajectory Calculations

Massive-Parallel Trajectory Calculations (MPTRAC) is a Lagrangian particle dispersion model for the analysis of atmospheric transport processes in the free troposphere and stratosphere.

![logo](https://github.com/slcs-jsc/mptrac/blob/master/docs/logo/MPTRAC_320px.png)

[![release (latest by date)](https://img.shields.io/github/v/release/slcs-jsc/mptrac)](https://github.com/slcs-jsc/mptrac/releases)
[![commits since latest release (by SemVer)](https://img.shields.io/github/commits-since/slcs-jsc/mptrac/latest)](https://github.com/slcs-jsc/mptrac/commits/master)
[![last commit](https://img.shields.io/github/last-commit/slcs-jsc/mptrac.svg)](https://github.com/slcs-jsc/mptrac/commits/master)
[![top language](https://img.shields.io/github/languages/top/slcs-jsc/mptrac.svg)](https://github.com/slcs-jsc/mptrac/tree/master/src)
[![code size in bytes](https://img.shields.io/github/languages/code-size/slcs-jsc/mptrac.svg)](https://github.com/slcs-jsc/mptrac/tree/master/src)
[![codacy](https://api.codacy.com/project/badge/Grade/a9de7b2239f843b884d2a4eb583726c9)](https://app.codacy.com/gh/slcs-jsc/mptrac?utm_source=github.com&utm_medium=referral&utm_content=slcs-jsc/mptrac&utm_campaign=Badge_Grade_Settings)
[![codecov](https://codecov.io/gh/slcs-jsc/mptrac/branch/master/graph/badge.svg?token=4X6IEHWUBJ)](https://codecov.io/gh/slcs-jsc/mptrac)
[![tests](https://img.shields.io/github/actions/workflow/status/slcs-jsc/mptrac/tests.yml?branch=master&label=tests)](https://github.com/slcs-jsc/mptrac/actions)
[![docs](https://img.shields.io/github/actions/workflow/status/slcs-jsc/mptrac/docs.yml?branch=master&label=docs)](https://slcs-jsc.github.io/mptrac)
[![license](https://img.shields.io/github/license/slcs-jsc/mptrac.svg)](https://github.com/slcs-jsc/mptrac/blob/master/COPYING)
[![doi](https://zenodo.org/badge/DOI/10.5281/zenodo.4400597.svg)](https://doi.org/10.5281/zenodo.4400597)

## Features

* MPTRAC calculates air parcel trajectories by solving the kinematic equation of motion using given horizontal wind and vertical velocity fields from global reanalyses or forecasts.
* Mesoscale diffusion and subgrid-scale wind fluctuations are simulated using the Langevin equation to add stochastic perturbations to the trajectories. A new inter-parcel exchange module represents mixing of air.
* Additional modules are implemented to simulate convection, sedimentation, exponential decay, gas and aqueous phase chemistry, and wet and dry deposition.
* Meteorological data pre-processing code provides estimates of the boundary layer, convective available potential energy, geopotential heights, potential vorticity, and tropopause data.
* Various output methods for particle, grid, ensemble, profile, sample, and station data. Gnuplot and ParaView interfaces for visualization.
* MPI-OpenMP-OpenACC hybrid parallelization and distinct code optimizations for efficient use from single workstations to HPC and GPU systems.
* Distributed open source under the terms of the GNU GPL.

## Getting started

### Prerequisites

This README file describes how to install MPTRAC on a Linux system.

The following software dependencies are required to compile MPTRAC:

* the [GNU make](https://www.gnu.org/software/make) build tool
* the C compiler of the [GNU Compiler Collection (GCC)](https://gcc.gnu.org)
* the [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl) for numerical calculations
* the [netCDF library](http://www.unidata.ucar.edu/software/netcdf) for file-I/O

The following optional software is required to enable additional features of MPTRAC:

* the distributed version control system [Git](https://git-scm.com/) to access the code repository
* the [HDF5 library](https://www.hdfgroup.org/solutions/hdf5) to enable the netCDF4 file format
* the [Zstandard library](https://facebook.github.io/zstd) and the [zfp library](https://computing.llnl.gov/projects/zfp) for compressed meteo data
* the [NVIDIA HPC Software Development Kit](https://developer.nvidia.com/hpc-sdk) for GPU support
* an MPI library such as [OpenMPI](https://www.open-mpi.org) or [ParaStationMPI](https://github.com/ParaStation/psmpi) for HPC support
* the graphing utility [gnuplot](http://www.gnuplot.info) for visualization

Some of the software is provided along with the MPTRAC repository, please see next section.

### Installation

Start by downloading the latest or one of the previous [MPTRAC releases](https://github.com/slcs-jsc/mptrac/releases). Unzip the release file:

    unzip mptrac-x.y.zip

Alternatively, you can get the development version of the software from the GitHub repository:

    git clone https://github.com/slcs-jsc/mptrac.git

Several libraries shipped along with MPTRAC can be compiled and installed by running a build script:

    cd [mptrac_directory]/libs
    ./build.sh -a

Then change to the source directory and edit the `Makefile` according to your needs:

    cd [mptrac_directory]/src
    emacs Makefile

In particular, you may want to check:

* Edit the `LIBDIR` and `INCDIR` paths to point to the directories where the GSL, netCDF, and other libraries are located on your system.

* By default, the MPTRAC binaries are linked statically, i.e., they can be copied and used on other machines. However, sometimes static compilation causes problems, e.g., in combination with dynamically compiled GSL and netCDF libraries or when using MPI or OpenACC. In this case, disable the `STATIC` flag and remember to set the `LD_LIBRARY_PATH` to include the paths to the shared libraries.

* To make use of the MPI parallelization of MPTRAC, the `MPI` flag must be enabled. Further steps will require an MPI library such as OpenMPI to be available on your system. To make use of the OpenACC parallelization, the `GPU` flag must be enabled. The NVIDIA HPC SDK is required to compile the GPU code. MPTRAC's OpenMP parallelization is always enabled.

Next, try compiling the code:

    make [-j]

To run the test cases to check the installation, use

    make check

This will run a series of tests sequentially. It will stop if any of the tests fail. Please check the log messages.

### Run the example

A simple example is provided to illustrate how to simulate the dispersion of volcanic ash from the eruption of the Puyehue-Cordón Caulle volcano, Chile, in June 2011.

The example can be found in the `projects/example/` subdirectory. The `projects/` subdirectory can also be used to store the results of your own simulation experiments with MPTRAC.

The example simulation is controlled by a shell script:

    cd mptrac/projects/example
    ./run.sh

See the `run.sh` script for how to invoke MPTRAC programs such as `atm_init` and `atm_split` to initialize the trajectory seeds and `trac` to compute the trajectories.

The script generates simulation output in the `examples/data` subdirectory. The corresponding reference data can be found in `examples/data.ref`.

A set of plots of the simulation output at different time steps after the eruption, generated by means of the `gnuplot` plotting tool, can be found in `examples/plots`. The plots should look similar to the output provided in `examples/plots.ref`.

This is an example showing the particle positions and grid output on 6th and 8th of June 2011:
<p align="center"><img src="projects/example/plots.ref/atm_2011_06_06_00_00.tab.png" width="45%"/> &emsp; <img src="projects/example/plots.ref/grid_2011_06_06_00_00.tab.png" width="45%"/></p>
<p align="center"><img src="projects/example/plots.ref/atm_2011_06_08_00_00.tab.png" width="45%"/> &emsp; <img src="projects/example/plots.ref/grid_2011_06_08_00_00.tab.png" width="45%"/></p>

## Further information

These are the main scientific publications that provide information about MPTRAC:

* Hoffmann, L., Baumeister, P. F., Cai, Z., Clemens, J., Griessbach, S., Günther, G., Heng, Y., Liu, M., Haghighi Mood, K., Stein, O., Thomas, N., Vogel, B., Wu, X., and Zou, L.: Massive-Parallel Trajectory Calculations version 2.2 (MPTRAC-2.2): Lagrangian transport simulations on graphics processing units (GPUs), Geosci. Model Dev., 15, 2731–2762, https://doi.org/10.5194/gmd-15-2731-2022, 2022.

* Hoffmann, L., T. Rößler, S. Griessbach, Y. Heng, and O. Stein, Lagrangian transport simulations of volcanic sulfur dioxide emissions: Impact of meteorological data products, J. Geophys. Res. Atmos., 121, 4651-4673, https://doi.org/10.1002/2015JD023749, 2016. 

Additional references are collected on the [references web page](https://slcs-jsc.github.io/mptrac/references/).

More detailed information for users of MPTRAC is provided in the [user manual](https://slcs-jsc.github.io/mptrac).

Information for developers of MPTRAC can be found in the [doxygen manual](https://slcs-jsc.github.io/mptrac/doxygen).

## Contributing

We are interested in supporting operational and research applications with MPTRAC.

You can submit bug reports or feature requests on the [issue tracker](https://github.com/slcs-jsc/mptrac/issues).

Proposed code changes and fixes can be submitted as [pull requests](https://github.com/slcs-jsc/mptrac/pulls).

Please do not hesitate to contact us if you have any questions or need assistance.

## License

MPTRAC is distributed under the terms of the [GNU General Public License v3.0](https://github.com/slcs-jsc/mptrac/blob/master/COPYING).

Please see the [citation file](https://github.com/slcs-jsc/mptrac/blob/master/CITATION.cff) for more information about citing the MPTRAC model in scientific publications.

## Contact

Dr. Lars Hoffmann

Jülich Supercomputing Centre, Forschungszentrum Jülich

e-mail: l.hoffmann@fz-juelich.de
