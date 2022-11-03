# Welcome to MPTRAC!

Massive-Parallel Trajectory Calculations (MPTRAC) is a Lagrangian particle dispersion model for the analysis of atmospheric transport processes in the free troposphere and stratosphere.

This manual provides information for users and developers of the model.

## Features

These are some of the most important features of the MPTRAC model:

* MPTRAC calculates air parcel trajectories by solving the kinematic equation of motion using given horizontal wind and vertical velocity fields.

* Mesoscale diffusion and sub-grid scale wind fluctuations are simulated using the Langevin equation to add stochastic perturbations to the trajectories.

* Additional modules are implemented to simulate convection, sedimentation, radioactive decay, hydroxyl chemistry, dry deposition, and wet deposition.

* Various output methods for particle, ensemble, gridded, sample, and station data. Gnuplot interface for direct visualization.

* MPTRAC features an MPI-OpenMP-OpenACC hybrid parallelization for efficient use on HPC and GPU systems.

* MPTRAC is distributed under the [GNU General Public License v3.0](https://github.com/slcs-jsc/mptrac/blob/master/COPYING).

## Contact

We are interested in sharing MPTRAC for operational and research applications. Please do not hesitate to contact us, if you have any questions or need support.

Dr. Lars Hoffmann

Jülich Supercomputing Centre, Forschungszentrum Jülich

e-mail: l.hoffmann@fz-juelich.de

website: https://www.fz-juelich.de/ias/jsc/slcs
