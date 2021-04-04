# Massive-Parallel Trajectory Calculations (MPTRAC)

Massive-Parallel Trajectory Calculations (MPTRAC) is a Lagrangian particle dispersion model for the analysis of atmospheric transport processes in the troposphere and stratosphere.

[![release (latest by date)](https://img.shields.io/github/v/release/slcs-jsc/mptrac)](https://github.com/slcs-jsc/mptrac/releases)
[![commits since latest release (by SemVer)](https://img.shields.io/github/commits-since/slcs-jsc/mptrac/latest)](https://github.com/slcs-jsc/mptrac/commits/master)
[![last commit](https://img.shields.io/github/last-commit/slcs-jsc/mptrac.svg)](https://github.com/slcs-jsc/mptrac/commits/master)
[![code size in bytes](https://img.shields.io/github/languages/code-size/slcs-jsc/mptrac.svg)](https://github.com/slcs-jsc/mptrac/tree/master/src)
[![top language](https://img.shields.io/github/languages/top/slcs-jsc/mptrac.svg)](https://github.com/slcs-jsc/mptrac/tree/master/src)
[![test case](https://img.shields.io/github/workflow/status/slcs-jsc/mptrac/testcase?label=test%20case)](https://github.com/slcs-jsc/mptrac/actions)
[![doxygen](https://img.shields.io/github/workflow/status/slcs-jsc/mptrac/doxygen?label=doxygen)](https://slcs-jsc.github.io/mptrac)
[![license](https://img.shields.io/github/license/slcs-jsc/mptrac.svg)](https://github.com/slcs-jsc/mptrac/blob/master/COPYING)
[![doi](https://zenodo.org/badge/DOI/10.5281/zenodo.4400597.svg)](https://doi.org/10.5281/zenodo.4400597)


![test](https://img.shields.io/github/workflow/status/slcs-jsc/mptrac/doxygen?label=doxygen)

![test2](https://img.shields.io/github/workflow/status/slcs-jsc/mptrac/testcase?label=test%20case)



## Features

* MPTRAC calculates particle trajectories by solving the kinematic equation of motion using given wind fields.
* Mesoscale diffusion and sub-grid scale wind fluctuations are simulated using the Langevin equation to add stochastic perturbations to the trajectories.
* Additional modules are implemented to simulate convection, sedimentation, radioactive decay, hydroxyl chemistry, dry deposition, and wet deposition.
* Various output methods for particle, ensemble, gridded, sample, and station data. Gnuplot interface for direct visualization.
* MPTRAC features an MPI/OpenMP/OpenACC hybrid parallelization for efficient use on HPC and GPU systems.

## Getting started

### Prerequisites

This documentation describes the installation of MPTRAC on a Linux system. A number of standard tools (gcc, git, make) and software libraries are needed to install MPTRAC. The [GNU Scientific Library](https://www.gnu.org/software/gsl) is required for numerical calculations and the [Unidata netCDF library](http://www.unidata.ucar.edu/software/netcdf), the [HDF5 library](https://www.hdfgroup.org/solutions/hdf5), and [zlib](http://www.zlib.net/) are needed for file-I/O. Copies of these libraries are provided in the MPTRAC git repository.

Start by downloading the MPTRAC source code from the git repository:

    git clone https://github.com/slcs-jsc/mptrac.git

To update an existing installation use:

    git pull https://github.com/slcs-jsc/mptrac.git

### Installation

First, compile the GSL, netCDF, HDF5, and zlib libraries required by MPTRAC by running the build script:

    cd mptrac/libs
    ./build.sh <nc2|nc4>

Please select `nc2`, if you want to use meteorological data files in netCDF classic format, or select `nc4`, if you want to use both, netCDF classic and netCDF-4 data files. The HDF5 and zlib libraries are needed only for netCDF-4. Sometimes, the compilation of netCDF-4 may fail, and netCDF classic may be useful as a fall-back option in that case.

Next, change to the source directory and edit the Makefile according to your needs.

    cd mptrac/src
    emacs Makefile

In particular, comment or uncomment the `NC4` flag, depending on whether you want to use netCDF classic or netCDF-4 data files. You may also want to edit the LIBDIR and INCDIR paths to point to the directories where the libraries are located on your system. By default, LIBDIR and INCDIR will point to `../libs/build/lib` and `../libs/build/include`, respectively.

To make use of the MPI parallelization of MPTRAC, the MPI flag needs to be uncommented in the Makefile. Further steps of the installation will require an MPI library to be installed or loaded as a module. To make use of the OpenACC parallelization, the GPU flag needs to be uncommented. The PGI Compiler Suite will be required to compile the GPU code. The OpenMP parallelization of MPTRAC is always enabled.

Load any software modules that might also be needed on your target platform, and try to compile the code:

    make [-j4]

The argument `-j` is optional. It can be used to specify the number of parallel threads to speed up compilation.

After compilation, the MPTRAC binaries are located in the mptrac/src/ directory.

By default, the binaries will be linked statically, i.e., they can be copied and run on other machines. However, sometimes static compilations causes problems, in particular in combination with MPI and OpenACC, as static versions of some libraries might be missing. In this case, remove the `-static` flag from the CFLAGS in the Makefile and compile again. To run dynamically linked binaries, the LD_LIBRARY_PATH needs to be set to include the mptrac/libs/build/lib directory.

By default we apply rather strict compiler warnings to catch problems. Also, all warning messages will be turned into errors and no binaries will be produced. This behavior is enforced by the flag `-Werror`. It should not be removed from the Makefile, unless you know what you are doing.

### Run the example

An example is provided, illustrating how to simulate the dispersion of volcanic ash from the eruption of the Puyehue-Cordón Caulle volcano, Chile, in June 2011.

It is recommended that you create a project directory for testing the example and to store the results also of other experiments:

    mkdir -p mptrac/projects
    cp -a mptrac/example mptrac/projects
    
This shows how to run the example:

    cd mptrac/projects/example
    ./run.sh

At first call, the run script will download meteorological input data from a data server. This step may take a while as the input data comprise several hundred MByte in size. The input data are saved for later runs and need to be downloaded only once.

Please see the example script (run.sh) on how to invoke MPTRAC programs such as `atm_init` and `atm_split` to initialize trajectory seeds and `trac` to calculate the trajectories.

The script generates a number of plots of the simulation output at different times after the eruption by means of `gnuplot`. These plots should look similar to the output already provided in the repository.

This is an example showing the particle position and grid output on 6th and 8th of June 2011:
<p align="center"><img src="example/plots.ref/atm_2011_06_06_00_00.tab.png" width="45%"/> &emsp; <img src="example/plots.ref/grid_2011_06_06_00_00.tab.png" width="45%"/></p>
<p align="center"><img src="example/plots.ref/atm_2011_06_08_00_00.tab.png" width="45%"/> &emsp; <img src="example/plots.ref/grid_2011_06_08_00_00.tab.png" width="45%"/></p>

## Further information

More detailed information for new users and developers of MPTRAC is collected here:

* [GitHub wiki](https://github.com/slcs-jsc/mptrac/wiki)

* [reference manual](https://slcs-jsc.github.io/mptrac)

These are the main references for citing the MPTRAC model in scientific publications:

* Hoffmann, L., T. Rößler, S. Griessbach, Y. Heng, and O. Stein, Lagrangian transport simulations of volcanic sulfur dioxide emissions: Impact of meteorological data products, J. Geophys. Res. Atmos., 121, 4651-4673, https://doi.org/10.1002/2015JD023749, 2016. 

* You can cite the source code of MPTRAC by using the DOI https://doi.org/10.5281/zenodo.4400597. This DOI represents all versions, and will always resolve to the latest one. Specific DOIs for each release of MPTRAC can be found on the zenodo web site.

A [list of selected papers](https://github.com/slcs-jsc/mptrac/wiki/References) is also available from the GitHub wiki.

## Contributing

We are interested in sharing MPTRAC for operational or research applications.

Please do not hesitate to contact us, if you have any further questions or need support.

## License

MPTRAC is distributed under the [GNU General Public License v3.0](https://github.com/slcs-jsc/mptrac/blob/master/COPYING).

## Contact

Lars Hoffmann  

Jülich Supercomputing Centre, Forschungszentrum Jülich

e-mail: l.hoffmann@fz-juelich.de

website: https://www.fz-juelich.de/ias/jsc/slcs
