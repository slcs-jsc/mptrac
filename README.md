# Massive-Parallel Trajectory Calculations (MPTRAC)

Massive-Parallel Trajectory Calculations (MPTRAC) is a Lagrangian particle dispersion model for the analysis of atmospheric transport processes in the free troposphere and stratosphere.

![GitHub tag (latest SemVer)](https://img.shields.io/github/tag/slcs-jsc/mptrac.svg)
![GitHub top language](https://img.shields.io/github/languages/top/slcs-jsc/mptrac.svg)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/slcs-jsc/mptrac.svg)
![GitHub last commit](https://img.shields.io/github/last-commit/slcs-jsc/mptrac.svg)
![GitHub](https://img.shields.io/github/license/slcs-jsc/mptrac.svg)

## Features

* MPTRAC calculates air parcel trajectories by solving the kinematic equation of motion using given wind fields.
* Mesoscale diffusion and sub-grid scale wind fluctuations are simulated with a Markov chain model.
* Additional modules are implemented to simulate the sedimentation of air parcels and the decay of particle mass.
* MPTRAC features an MPI/OpenMP hybrid parallelization for efficient use on supercomputers.

Further information can be found at:
http://www.fz-juelich.de/ias/jsc/mptrac

## Installation

This documentation describes the installation of MPTRAC on a Linux system.
A number of standard tools such as the GNU Compiler Collection (gcc)
and 'make' are required to install MPTRAC.

Start by downloading the source code from the github repository:

    git clone https://github.com/slcs-jsc/mptrac

Change to the directory mptrac/ which holds source codes,
libraries, documentation, etc:

    cd mptrac

The GNU Scientific Library (https://www.gnu.org/software/gsl)
is required for numerical calculations and the Unidata netCDF library
(http://www.unidata.ucar.edu/software/netcdf) is needed for file-I/O.
Copies of these libraries can be found in the repository, if they are
not available on your system. A script is provided to build the libraries:

    cd lib
    ./build.sh

Next, change to the source directory and edit the Makefile according to
your needs. In particular, check the paths to the libraries
(INCDIR and LIBDIR). Then try to compile the code:

    cd ../src
    emacs Makefile
    make

The binaries will be linked statically, i.e., they can be copied to other
machines. Sometimes static compilations causes problems, in particular in
combination with MPI. In this case remove the '-static' flag from the
CFLAGS in the Makefile and compile again.

By default we use rather strict compiler warnings.
All warning messages will be turned into errors and no binaries will be
produced. This behavior is enforced by the flag '-Werror'.

The binaries will remain in the src/ directory.

## Getting started

This script illustrates how to use MPTRAC:

    cd ../example
    ./run.sh

The example illustrates how to simulate the dispersion of volcanic ash from the
eruption of the Puyehue-Cordón Caulle volcano, Chile, in June 2011.

The script generates a number of plots of the simulation output
at different times after the eruption by means of 'gnuplot'.
These plots should look similar to the output already
provided in the repository.

More details on the control parameters, data structures, and algorithms
can be found in the MPTRAC reference manual:

     evince ../doc/refman.pdf

## Contact

We are interested in sharing MPTRAC for research applications.

Please do not hesitate to contact us if you have any further questions:

Dr. Lars Hoffmann  
Forschungszentrum Jülich  
Jülich Supercomputing Centre  
52425 Jülich  
Germany  

e-mail: l.hoffmann@fz-juelich.de

## License

MPTRAC is distributed under the GNU GPL v3.
