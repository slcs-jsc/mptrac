# Massive-Parallel Trajectory Calculations

Massive-Parallel Trajectory Calculations (MPTRAC) is a Lagrangian particle dispersion model for the analysis of atmospheric transport processes in the free troposphere and stratosphere.

[![release (latest by date)](https://img.shields.io/github/v/release/slcs-jsc/mptrac)](https://github.com/slcs-jsc/mptrac/releases)
[![commits since latest release (by SemVer)](https://img.shields.io/github/commits-since/slcs-jsc/mptrac/latest)](https://github.com/slcs-jsc/mptrac/commits/master)
[![last commit](https://img.shields.io/github/last-commit/slcs-jsc/mptrac.svg)](https://github.com/slcs-jsc/mptrac/commits/master)
[![code size in bytes](https://img.shields.io/github/languages/code-size/slcs-jsc/mptrac.svg)](https://github.com/slcs-jsc/mptrac/tree/master/src)
[![top language](https://img.shields.io/github/languages/top/slcs-jsc/mptrac.svg)](https://github.com/slcs-jsc/mptrac/tree/master/src)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/a9de7b2239f843b884d2a4eb583726c9)](https://app.codacy.com/gh/slcs-jsc/mptrac?utm_source=github.com&utm_medium=referral&utm_content=slcs-jsc/mptrac&utm_campaign=Badge_Grade_Settings)
[![codecov](https://codecov.io/gh/slcs-jsc/mptrac/branch/master/graph/badge.svg?token=4X6IEHWUBJ)](https://codecov.io/gh/slcs-jsc/mptrac)
[![test case](https://img.shields.io/github/workflow/status/slcs-jsc/mptrac/testcase?label=test%20case)](https://github.com/slcs-jsc/mptrac/actions)
[![doxygen](https://img.shields.io/github/workflow/status/slcs-jsc/mptrac/doxygen?label=doxygen)](https://slcs-jsc.github.io/mptrac)
[![license](https://img.shields.io/github/license/slcs-jsc/mptrac.svg)](https://github.com/slcs-jsc/mptrac/blob/master/COPYING)
[![doi](https://zenodo.org/badge/DOI/10.5281/zenodo.4400597.svg)](https://doi.org/10.5281/zenodo.4400597)

## Features

* MPTRAC calculates air parcel trajectories by solving the kinematic equation of motion using given horizontal wind and vertical velocity fields.
* Mesoscale diffusion and sub-grid scale wind fluctuations are simulated using the Langevin equation to add stochastic perturbations to the trajectories.
* Additional modules are implemented to simulate convection, sedimentation, radioactive decay, hydroxyl chemistry, dry deposition, and wet deposition.
* Various output methods for particle, ensemble, gridded, sample, and station data. Gnuplot interface for direct visualization.
* MPTRAC features an MPI/OpenMP/OpenACC hybrid parallelization for efficient use on HPC and GPU systems.

## Getting started

### Prerequisites

This documentation describes the installation of MPTRAC on a Linux system. Some standard tools (e.g., gcc, git, make) are needed for compilation and installation. The [GNU Scientific Library](https://www.gnu.org/software/gsl) is required for numerical calculations and the [Unidata netCDF library](http://www.unidata.ucar.edu/software/netcdf) is needed for file-I/O. The graphing utility [gnuplot](http://www.gnuplot.info) is recommend for visualization.

Commands to install the necessary software on Ubuntu Linux:

    sudo apt-get update
    sudo apt-get install gcc git make gnuplot libgsl-dev libnetcdf-dev libhdf5-dev

Additional software being recommended for code development:

    sudo apt-get install cppcheck doxygen graphviz indent lcov

### Installation

Start by downloading the MPTRAC source code from the git repository:

    git clone https://github.com/slcs-jsc/mptrac.git

To update an existing installation, please use:

    git pull https://github.com/slcs-jsc/mptrac.git

In case the GSL and netCDF libraries are missing on your system, it is recommended to install the versions of the libraries provided along with MPTRAC by running the library build script:

    cd mptrac/libs
    ./build.sh

Next, change to the source directory and edit the Makefile according to your needs.

    cd mptrac/src
    emacs Makefile

In particular, you might want to check:

* Edit the `LIBDIR` and `INCDIR` paths to point to the directories where the GSL and netCDF libraries are located on your system.

* By default, the MPTRAC binaries will be linked statically, i.e., they can be copied and used on other machines. However, sometimes static compilations causes problems, e.g., in combination with dynamically compiled netCDF and GSL libraries or when using MPI and OpenACC. In this case, disable the `STATIC` flag in the Makefile and remember to set the `LD_LIBRARY_PATH` to include the libraries.

* To make use of the MPI parallelization of MPTRAC, the `MPI` flag needs to be enabled in the Makefile. Further steps of the installation will require an MPI library such as [OpenMPI](https://www.open-mpi.org) or [MPICH](https://www.mpich.org).

* To make use of the OpenACC parallelization, the `GPU` flag needs to be enabled. The [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk) will be required to compile the GPU code.

Try to compile the code:

    make

To run a test case of the installation, please use:

    make check

After compilation, the MPTRAC binaries are located in the directory `mptrac/src/`. To install them to `/usr/local/bin` (or another directory specified by `DESTDIR`), please use:

    make install

### Run the example

An example is provided, illustrating how to simulate the dispersion of volcanic ash from the eruption of the Puyehue-Cordón Caulle volcano, Chile, in June 2011.

It is recommended that you create a project directory for testing the example and also to store the results of other experiments:

    mkdir -p mptrac/projects
    cp -a mptrac/example mptrac/projects
    
This shows how to run the example:

    cd mptrac/projects/example
    ./run.sh

Please see the example script `run.sh` on how to invoke MPTRAC programs such as `atm_init` and `atm_split` to initialize trajectory seeds and `trac` to calculate the trajectories.

The script generates a number of plots of the simulation output at different time steps after the eruption by means of `gnuplot`. These plots should look similar to the output already provided in the repository.

This is an example showing the particle position and grid output on 6th and 8th of June 2011:
<p align="center"><img src="example/plots.ref/atm_2011_06_06_00_00.tab.png" width="45%"/> &emsp; <img src="example/plots.ref/grid_2011_06_06_00_00.tab.png" width="45%"/></p>
<p align="center"><img src="example/plots.ref/atm_2011_06_08_00_00.tab.png" width="45%"/> &emsp; <img src="example/plots.ref/grid_2011_06_08_00_00.tab.png" width="45%"/></p>

## Further information

More detailed information for new users and developers of MPTRAC is collected in the [GitHub wiki](https://github.com/slcs-jsc/mptrac/wiki).

These are the main references for citing the MPTRAC model in scientific publications:

* Hoffmann, L., Baumeister, P. F., Cai, Z., Clemens, J., Griessbach, S., Günther, G., Heng, Y., Liu, M., Haghighi Mood, K., Stein, O., Thomas, N., Vogel, B., Wu, X., and Zou, L.: Massive-Parallel Trajectory Calculations version 2.2 (MPTRAC-2.2): Lagrangian transport simulations on graphics processing units (GPUs), Geosci. Model Dev., 15, 2731–2762, https://doi.org/10.5194/gmd-15-2731-2022, 2022.

* Hoffmann, L., T. Rößler, S. Griessbach, Y. Heng, and O. Stein, Lagrangian transport simulations of volcanic sulfur dioxide emissions: Impact of meteorological data products, J. Geophys. Res. Atmos., 121, 4651-4673, https://doi.org/10.1002/2015JD023749, 2016. 

* You can cite the source code of MPTRAC by using the DOI https://doi.org/10.5281/zenodo.4400597. This DOI represents all versions, and will always resolve to the latest one. Specific DOIs for each release of MPTRAC can be found on the zenodo web site.

Please see the [citation file](https://github.com/slcs-jsc/mptrac/blob/master/CITATION.cff) for further information.

## Contributing

We are interested in sharing MPTRAC for operational and research applications. Please do not hesitate to contact us, if you have any further questions or need support.

## License

MPTRAC is distributed under the [GNU General Public License v3.0](https://github.com/slcs-jsc/mptrac/blob/master/COPYING).

## Contact

Dr. Lars Hoffmann

Jülich Supercomputing Centre, Forschungszentrum Jülich

e-mail: l.hoffmann@fz-juelich.de

website: https://www.fz-juelich.de/ias/jsc/slcs
