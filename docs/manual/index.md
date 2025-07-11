# Welcome to MPTRAC!

Massive-Parallel Trajectory Calculations (MPTRAC) is a Lagrangian
particle dispersion model to analyze atmospheric transport processes
in the free troposphere and stratosphere. Leveraging high-performance
computing techniques, MPTRAC efficiently handles large-scale
trajectory simulations, making it a powerful tool for both research
and operational applications.

![Lagrangian transport simulation of convective transport](img/convection_lowres.jpg)

## Features

MPTRAC is a powerful tool for atmospheric Lagrangian transport
simulations, offering a wide range of features to enhance accuracy,
performance, and usability:

- **Advanced Trajectory Calculations**: MPTRAC calculates air parcel
    trajectories by solving the kinematic equation of motion using
    horizontal wind and vertical velocity fields from global
    reanalyses or forecast datasets, enabling precise tracking of
    atmospheric transport processes in the free troposphere and
    stratosphere.

- **Stochastic Perturbation and Mixing**: Mesoscale diffusion and
    subgrid-scale wind fluctuations are simulated using the Langevin
    equation, introducing stochastic perturbations to trajectories. An
    inter-parcel exchange module represents mixing of air between
    neighboring particles, capturing realistic atmospheric dispersion.

- **Comprehensive Process Modeling**: MPTRAC includes modules to
    simulate convection, sedimentation, exponential decay, gas and
    aqueous phase chemistry, and wet and dry deposition, allowing for
    accurate modeling of physical and chemical transformations.

- **Meteorological Data Pre-Processing**: The model pre-processes
    meteorological data to estimate variables such as boundary layer
    height, convective available potential energy (CAPE), geopotential
    heights, potential vorticity, and tropopause data, ensuring
    seamless integration with diverse datasets.

- **Flexible Output and Visualization**: MPTRAC supports various
    output formats for particle trajectories, gridded fields, ensemble
    statistics, vertical profiles, point samples, and station
    data. Visualization interfaces with Gnuplot and ParaView make it
    easy to analyze complex data.

- **High-Performance Computing**: The model employs hybrid
    parallelization using MPI, OpenMP, and OpenACC, allowing efficient
    utilization of resources from single workstations to HPC clusters
    and GPU-based systems.

- **Web-Based Accessibility**: The new MPTRAC Web Runner provides an
    intuitive, browser-based interface for running trajectory
    simulations without local installation, making the tool more
    accessible for educational, research, and operational users.

- **Open Source and Community Driven**: MPTRAC is distributed as
    open-source software under the GNU General Public License (GPL),
    promoting collaborative development and ensuring transparency.

![Geophysical modules and main software components of MPTRAC](img/clusters.png)

## Get Involved

We encourage collaboration and welcome contributions from the research
community. Feel free to reach out if you have questions, suggestions,
or need support with MPTRAC.

## Contact

Dr. Lars Hoffmann

Jülich Supercomputing Centre, Forschungszentrum Jülich, Germany

e-mail: <l.hoffmann@fz-juelich.de>
