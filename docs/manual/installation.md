# Installation

This guide covers the installation of MPTRAC on a Linux system. The
process involves several steps: meeting prerequisites, downloading the
source code, building the required libraries, compiling the source
code, and verifying the installation.

## Prerequisites

To build and run MPTRAC, you will need some basic tools and libraries,
including [Git](https://git-scm.com/),
[GNU Make](https://www.gnu.org/software/make),
[GCC](https://gcc.gnu.org/),
[GSL](https://www.gnu.org/software/gsl),
[HDF5](https://www.hdfgroup.org/solutions/hdf5), and
[netCDF](https://www.unidata.ucar.edu/software/netcdf).

For additional features such as HPC and GPU support, optional
dependencies like [OpenMPI](https://www.open-mpi.org) and the
[NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk) are required.

Some of the required dependencies are included with the MPTRAC
repository. See the next section for more details.

For a complete list of dependencies, including specific versions and
installation instructions, refer to the
[dependencies file](https://github.com/slcs-jsc/mptrac/blob/master/DEPENDENCIES.md).

## Downloading the source code

Get the latest or a previous version from the
[MPTRAC releases](https://github.com/slcs-jsc/mptrac/releases)
page. After downloading, extract the release file:

    unzip mptrac-x.y.zip

Alternatively, to get the latest development version, clone the GitHub
repository:

    git clone https://github.com/slcs-jsc/mptrac.git

To update an existing installation, navigate to the directory and pull
the latest changes:

    git pull https://github.com/slcs-jsc/mptrac.git

## Building the libraries

MPTRAC includes several libraries that can be compiled and installed
using a build script:

    cd [mptrac_directory]/libs
    ./build.sh -a

The build process of the libraries can take considerable time. If your
system already has compatible versions of the libraries, consider
using those instead.

## Configure the Makefile

Navigate to the source directory and adjust the `Makefile` as needed:

    cd [mptrac_directory]/src
    emacs Makefile

Pay special attention to the following settings:

* Edit the `LIBDIR` and `INCDIR` paths to point to the directories
  where the GSL, netCDF, and other required libraries are located on
  your system.

* By default, the MPTRAC binaries are linked dynamically. Ensure that
  the `LD_LIBRARY_PATH` is properly configured to include the paths to
  the shared libraries. If you prefer static linking, you can enable
  it by setting the `STATIC` flag, which allows you to copy and use
  the binaries on other machines. However, in some cases, either
  static or dynamic linking may not be feasible or could cause
  specific issues.

* To enable MPI parallelization in MPTRAC, you must set the `MPI`
  flag. Additionally, an MPI library, such as OpenMPI, must be
  installed on your system. To utilize OpenACC parallelization, enable
  the `GPU` flag, and ensure the NVIDIA HPC SDK is installed to
  compile the GPU code. OpenMP parallelization is always enabled.

* Some options in the Makefile are labeled as experimental. These
  features are still under development and may not be fully functional
  or tested. Use them at your own risk.

## Compile and test the installation

Once the Makefile is configured, compile the code using:

    make [-j]

To verify the installation, run the test suite:

    make check

These tests cover a wide range of MPTRAC functionalities. Make sure
all tests pass successfully before using the model in other
applications. If any test fails, check the log messages for further
details.
