# Installation

Here we describe the installation of MPTRAC on a Linux system. The installation comprises several steps, i.e. to meet the prerequisites, downloading the source code from the GitHub repository, the compilation of the source code, and conducting checks to confirm the integrity of the installation.

## Prerequisites

The following software dependencies are mandatory to compile and install MPTRAC:

* the distributed version control system [Git](https://git-scm.com/)
* the C compiler of the [GNU Compiler Collection](https://gcc.gnu.org)
* the [GNU Make](https://www.gnu.org/software/make) build tool
* the [GNU Scientific Library](https://www.gnu.org/software/gsl) for numerical calculations
* the [netCDF library](http://www.unidata.ucar.edu/software/netcdf) for file-I/O

Optionally, the following software is needed to enable further capabilities of the model:

* the graphing utility [gnuplot](http://www.gnuplot.info) for visualization
* the [HDF5 library](https://www.hdfgroup.org/solutions/hdf5) to support netCDF4
* the [Zstandard library](https://facebook.github.io/zstd) and the [zfp library](https://computing.llnl.gov/projects/zfp) for compressed meteo data
* the [NVIDIA HPC Software Development Kit](https://developer.nvidia.com/hpc-sdk) for GPU support
* an MPI library such as [OpenMPI](https://www.open-mpi.org) or [ParaStation MPI](https://github.com/ParaStation/psmpi) for HPC support

Some of the software is provided along with the MPTRAC repository, please see next section.

## Downloading the source code

Start by downloading the MPTRAC source code from the GitHub repository:

    git clone https://github.com/slcs-jsc/mptrac.git

To update an existing installation, please use:

    git pull https://github.com/slcs-jsc/mptrac.git

## Building the libraries

Several libraries provided along with MPTRAC can be compiled and installed by running a build script:

    cd mptrac/libs
    ./build.sh

As the build process is rather time-consuming, it might be more efficient to apply versions of the libraries readily available on your system.

## Compilation of the source code

Change to the source directory and edit the Makefile according to your needs.

    cd mptrac/src
    emacs Makefile

In particular, you might want to edit the `LIBDIR` and `INCDIR` paths to point to the directories where GSL, netCDF, and other libraries are located on your system.

By default, the MPTRAC binaries will be linked statically, i.e., they can be copied and used on other machines. However, sometimes static compilations causes problems, e.g., in combination with dynamically compiled GSL and netCDF libraries or when using MPI and OpenACC. In this case, disable the `STATIC` flag and remember to set the `LD_LIBRARY_PATH` to include the libraries.

To make use of the MPI parallelization of MPTRAC, the `MPI` flag needs to be enabled. Further steps will require an MPI library such as OpenMPI to be available. To make use of the OpenACC parallelization, the `GPU` flag needs to be enabled. The NVIDIA HPC SDK is required to compile the GPU code. The OpenMP parallelization of MPTRAC is always enabled.

Finally, use `make` in the source directory to compile the code. Please carefully check the logs of the compilation process for any warnings or errors.

## Run the test cases

Please run the test cases in order to check the installation:

    make check

The test cases have been designed to cover a wide range of the functionality of the MPTRAC model. Please carefully inspect whether the test cases are conducted successfully before applying the model in any other application.
