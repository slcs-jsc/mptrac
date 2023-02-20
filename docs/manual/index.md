# Welcome to MPTRAC!

Massive-Parallel Trajectory Calculations (MPTRAC) is a Lagrangian particle dispersion model for the analysis of atmospheric transport processes in the free troposphere and stratosphere.

## Features

These are some of the most important features of the MPTRAC model:

* MPTRAC calculates air parcel trajectories by solving the kinematic equation of motion using given horizontal wind and vertical velocity fields.

* Mesoscale diffusion and sub-grid scale wind fluctuations are simulated using the Langevin equation to add stochastic perturbations to the trajectories.

* Additional modules are implemented to simulate convection, sedimentation, exponential decay, gas and aqueous phase chemistry, wet and dry deposition.

* Meteo data pre-processing code to calculate boundary layer, convective available potential energy, geopotential heights, potential vorticity, and tropopause data.

* Various output methods for particle, ensemble, gridded, sample, and station data. Gnuplot interface for direct visualization.

* MPTRAC features an MPI-OpenMP-OpenACC hybrid parallelization for efficient use on HPC and GPU systems.

* Distributed open source under the terms and conditions of the GNU GPL.

## Contact

We are interested in sharing MPTRAC for research applications. Please do not hesitate to contact us, if you have any questions or need support.

Dr. Lars Hoffmann, <l.hoffmann@fz-juelich.de>

Jülich Supercomputing Centre, Forschungszentrum Jülich
