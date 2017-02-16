# Massive-Parallel Trajectory Calculations

Massive-Parallel Trajectory Calculations (MPTRAC) is a new
Lagrangian particle dispersion model for the free troposphere
and stratosphere. MPTRAC has been developed to support
the analysis of atmospheric transport processes.

## Algorithms

The primary task of MPTRAC is to calculate air parcel trajectories by
solving the kinematic equation of motion. Turbulent diffusion and
sub-grid scale wind fluctuations are simulated with Markov chain models.
Additional modules are implemented to simulate the sedimentation
of air parcels and the decay of particle mass.

## Implementation

MPTRAC is written in the C programming language. It uses the
GNU Scientific Library (https://www.gnu.org/software/gsl) for
numerical operations. The Unidata netCDF library
(http://www.unidata.ucar.edu/software/netcdf) is used for
file-I/O. MPTRAC can be used on a desktop PC, but it also features
an MPI/OpenMP hybrid parallelization for efficient use on supercomputers.

## Contact

We are interested in sharing MPTRAC for research applications.

Please do not hesitate to contact us if you have further questions:

Dr. Lars Hoffmann  
Forschungszentrum Jülich  
Jülich Supercomputing Centre  
52425 Jülich  
Germany

e-mail: l.hoffmann@fz-juelich.de

## License

MPTRAC is distributed under the GNU GPL v3.
