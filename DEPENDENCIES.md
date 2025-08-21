# MPTRAC Dependencies

This file lists the software dependencies required to build and run MPTRAC. It includes both **mandatory** and **optional** dependencies. Ensure that the correct versions are installed for optimal performance.

## Mandatory Dependencies

These dependencies are required to compile MPTRAC:

| Dependency | Version | Description |
|------------|---------|-------------|
| [Git](https://git-scm.com/) | v2.43.0 | Version control and code repository access |
| [GNU Make](https://www.gnu.org/software/make) | v4.3 | Build tool |
| [GNU Compiler Collection (GCC)](https://gcc.gnu.org) | v13.3.0 | C, C++ and Fortran compilers |
| [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl) | v2.7.1 | Numerical calculations |
| [netCDF library](http://www.unidata.ucar.edu/software/netcdf) | v4.9.2 | Portable file I/O |
| [HDF5 library](https://www.hdfgroup.org/solutions/hdf5) | v1.14.4-3 | Enable netCDF4 file format |

## Optional Dependencies

These dependencies enable additional features of MPTRAC. They are not required for basic compilation but enhance functionality like GPU support, parallel computing, and data compression.

| Dependency | Version | Description |
|------------|---------|-------------|
| [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk) | v25.1 | GPU support |
| [CMake](https://cmake.org/) | v3.28.3 | Build tool (required for ecCodes support) |
| [Clang (LLVM)](https://clang.llvm.org) | v18.1.8 | Alternative C/C++ compiler |
| [Intel Compilers (ICC and IFORT)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/compilers.html) | v2024.2.0  | High-performance C, C++, and Fortran compilers from Intel |
| [Open MPI](https://www.open-mpi.org) | v5.1.0 | MPI library for HPC support |
| [ParaStation MPI](https://github.com/ParaStation/psmpi) | v5.11.0 | MPI library for HPC support |
| [ecCodes](https://confluence.ecmwf.int/display/ECC/ecCodes+Home) | v2.38.3 | Encoding and decoding GRIB messages |
| [Kinetic PreProcessor (KPP)](https://github.com/KineticPreProcessor/KPP) | v3.2.0 | Preprocessor for generating chemical mechanism code |
| [Thrust](https://developer.nvidia.com/thrust) | v1.14.0 | C++ template library for parallel algorithms and data structures |
| [Zstandard](https://facebook.github.io/zstd) | v1.5.5 | Compression library |
| [zfp](https://computing.llnl.gov/projects/zfp) | v1.0.1 | Compression library |
| [SZ3](https://github.com/szcompressor/SZ3) | v3.2.1 | Compression library |
| [gnuplot](http://www.gnuplot.info) | v5.4.5 | Graphing utility for visualization |
| [ParaView](https://www.paraview.org) | v5.11.0 | Visualization tool for scientific data |

## Installing Dependencies

For detailed installation instructions on each of these dependencies, refer to the official documentation linked above. Below are a few common ways to install them:

- **Ubuntu/Debian:** Use `apt-get` for system-wide installs:
  ```bash
  sudo apt-get install git make gcc gfortran libgsl-dev libnetcdf-dev libhdf5-dev zlib1g-dev
  ```

- **CentOS/RHEL:** Use `yum` to install most dependencies:
  ```bash
  sudo yum install git make gcc gsl-devel netcdf-devel hdf5-devel zlib-devel
  ```
  
- **MacOS:** Use `brew` (Homebrew):
  ```bash
  brew install git make gcc gsl netcdf hdf5 zlib
  ```
  
## Notes on Optional Dependencies

  1. **NVIDIA HPC SDK**: Required for GPU support. Follow the official NVIDIA installation guide for your system.

  2. **OpenMPI**: For high-performance computing (HPC) and parallel computing. Install via your package manager or build from source.

  3. **KPP (Kinetic PreProcessor)**: If you're working with chemical mechanisms, KPP is needed for preprocessing. Follow the KPP documentation for installation.
