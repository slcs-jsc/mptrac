# Massive-Parallel Trajectory Calculations (MPTRAC)

Massive-Parallel Trajectory Calculations (MPTRAC) is a Lagrangian particle dispersion model for the analysis of atmospheric transport processes in the free troposphere and stratosphere.

[GitHub tag (latest SemVer)](https://img.shields.io/github/tag/slcs-jsc/jurassic.svg)
[GitHub last commit](https://img.shields.io/github/last-commit/slcs-jsc/jurassic.svg)
[GitHub code size in bytes](https://img.shields.io/github/languages/code-size/slcs-jsc/jurassic.svg)
[GitHub top language](https://img.shields.io/github/languages/top/slcs-jsc/jurassic.svg)
[GitHub](https://img.shields.io/github/license/slcs-jsc/jurassic.svg)

## Features

* MPTRAC calculates air parcel trajectories by solving the kinematic equation of motion using given wind fields.
* Mesoscale diffusion and sub-grid scale wind fluctuations are simulated using a Markov chain model.
* Additional modules are implemented to simulate the sedimentation of particles and the decay of air parcel mass.
* Various output methods for particle, ensemble, gridded, and station data. Gnuplot interface for direct visualization.
* MPTRAC features an MPI/OpenMP and MPI/GPU hybrid parallelization for efficient use on supercomputers.

## Getting started

### Prerequisites

This documentation describes the installation of MPTRAC on a Linux system. A number of standard tools (gcc, git, make) and software libraries are needed to install MPTRAC. The [GNU Scientific Library](https://www.gnu.org/software/gsl) is required for numerical calculations and the [Unidata netCDF library](http://www.unidata.ucar.edu/software/netcdf) is needed for file-I/O. Copies of these libraries can be found in the git repository.

Start by downloading the source code from the git repository:

    git clone https://github.com/slcs-jsc/mptrac.git

To update an existing installation use:

    git pull https://github.com/slcs-jsc/mptrac.git

### Installation

First, compile the netCDF and GSL libraries needed for MPTRAC by using the build script:

    cd mptrac/lib
    ./build.sh

Next, change to the source directory, edit the Makefile according to your needs, and try to compile the code:

    cd mptrac/src
    emacs Makefile
    make

The binaries will be linked statically, i.e., they can be copied to other machines. Sometimes static compilations causes problems, in particular in combination with MPI. In this case remove the '-static' flag from the CFLAGS in the Makefile and compile again.

By default we use rather strict compiler warnings. All warning messages will be turned into errors and no binaries will be produced. This behavior is enforced by the flag '-Werror'.

The binaries will remain in the mptrac/src/ directory.

### Try the example

An example is provided, illustrating how to simulate the dispersion of volcanic ash from the eruption of the Puyehue-Cordón Caulle volcano, Chile, in June 2011.

It is recommended that you create a project directory for testing the example and also to store other experiments:

    mkdir -p mptrac/projects

This shows how to run the example:

    cp -a mptrac/example mptrac/projects
    cd mptrac/projects/example
    ./run.sh

Please see the example script (run.sh) on how to invoke programs such as atm_init and atm_split to initialize trajectory seeds and trac to calculate the trajectories.

The script generates a number of plots of the simulation output at different times after the eruption by means of 'gnuplot'. These plots should look similar to the output already provided in the repository.

This is an example showing the particle position output for 7 June 2011:
<p align="center"><img src="example/plots.org/atm_diff_2011_06_07_00_00.tab.png" width="60%"/></p>

## Further information

This is the main reference for citing the MPTRAC model in scientific publications:

* Hoffmann, L., T. Rößler, S. Griessbach, Y. Heng, and O. Stein, Lagrangian transport simulations of volcanic sulfur dioxide emissions: Impact of meteorological data products, J. Geophys. Res. Atmos., 121, 4651-4673, https://doi.org/10.1002/2015JD023749, 2016. 

This is a list of selected papers in which MPTRAC was applied:

* Hoffmann, L., Günther, G., Li, D., Stein, O., Wu, X., Griessbach, S., Heng, Y., Konopka, P., Müller, R., Vogel, B., and Wright, J. S.: From ERA-Interim to ERA5: the considerable impact of ECMWF's next-generation reanalysis on Lagrangian transport simulations, Atmos. Chem. Phys., 19, 3097-3124, https://doi.org/10.5194/acp-19-3097-2019, 2019.

* Wu, X., Griessbach, S., and Hoffmann, L.: Long-range transport of volcanic aerosol from the 2010 Merapi tropical eruption to Antarctica, Atmos. Chem. Phys., 18, 15859-15877, https://doi.org/10.5194/acp-18-15859-2018, 2018.

* Rößler, T., Stein, O., Heng, Y., Baumeister, P., and Hoffmann, L.: Trajectory errors of different numerical integration schemes diagnosed with the MPTRAC advection module driven by ECMWF operational analyses, Geosci. Model Dev., 11, 575-592, https://doi.org/10.5194/gmd-11-575-2018, 2018.

* Wu, X., Griessbach, S., and Hoffmann, L.: Equatorward dispersion of a high-latitude volcanic plume and its relation to the Asian summer monsoon: a case study of the Sarychev eruption in 2009, Atmos. Chem. Phys., 17, 13439-13455, https://doi.org/10.5194/acp-17-13439-2017, 2017.

* Hoffmann, L., Hertzog, A., Rößler, T., Stein, O., and Wu, X.: Intercomparison of meteorological analyses and trajectories in the Antarctic lower stratosphere with Concordiasi superpressure balloon observations, Atmos. Chem. Phys., 17, 8045-8061, https://doi.org/10.5194/acp-17-8045-2017, 2017.

* Heng, Y., Hoffmann, L., Griessbach, S., Rößler, T., and Stein, O.: Inverse transport modeling of volcanic sulfur dioxide emissions using large-scale simulations, Geosci. Model Dev., 9, 1627-1645, https://doi.org/10.5194/gmd-9-1627-2016, 2016.

More details on the data structures and algorithms can be found in the [MPTRAC reference manual](doc/refman.pdf).

## Contributing

We are interested in sharing MPTRAC for operational or research applications.

Please do not hesitate to contact us, if you have any further questions or need support.

## License

MPTRAC is distributed under the GNU GPL v3.

## Contact

Dr. Lars Hoffmann  

Jülich Supercomputing Centre, Forschungszentrum Jülich

e-mail: l.hoffmann@fz-juelich.de

website: https://www.fz-juelich.de/ias/jsc/slcs
