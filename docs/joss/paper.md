---
title: 'MPTRAC: A high-performance Lagrangian transport model for atmospheric air parcel dispersion'

tags:
  - Atmospheric modeling
  - Lagrangian transport
  - Particle dispersion
  - Troposphere
  - Stratosphere
  - High-performance computing (HPC)
  - Graphical processing unit (GPU)
  - Numerical simulation
  - Advection and diffusion
  - Meteorology
  - Environmental monitoring
  - Aerosol transport
  - Pollutant tracking

authors:
  - name: Lars Hoffmann
    affiliation: "1, 2"
    orcid: 0000-0003-3773-4377
  - name: Jan Clemens
    affiliation: "2, 3"
    orcid: 0000-0002-2422-9454
  - name: Sabine Griessbach
    affiliation: "1, 2"
    orcid: 0000-0003-3792-3573
  - name: Kaveh Haghighi Mood
    affiliation: "1"
    orcid: 0000-0002-8578-4961
  - name: Yi Heng
    affiliation: "4"
    orcid: 0000-0001-6894-6976
  - name: Farahnaz Khosrawi
    affiliation: "1, 2"
    orcid: 0000-0002-0261-7253
  - name: Mingzhao Liu
    affiliation: "4"
    orcid: 0000-0001-7314-9568
  - name: Yen-Sen Lu
    affiliation: "1, 2"
    orcid: 0000-0002-3255-824X
  - name: Catrin Meyer
    affiliation: "1, 2"
    orcid: 0000-0002-9271-6174
  - name: Nils Nobre Wittwer
    affiliation: "1, 2"
    orcid: 0009-0001-3175-0889
  - name: Xue Wu
    affiliation: "5, 6"
    orcid: 0000-0002-0427-782X
  - name: Ling Zou
    affiliation: "1, 2"
    orcid: 0000-0001-6563-9815

affiliations:
 - name: Jülich Supercomputing Centre, Forschungszentrum Jülich, Jülich, Germany
   index: 1
 - name: Centre for Advanced Simulation and Analytics, Forschungszentrum Jülich, Jülich, Germany
   index: 2
 - name: Institute of Climate and Energy Systems, Forschungszentrum Jülich, Jülich, Germany
   index: 3
 - name: School of Computer Science and Engineering, Sun Yat-sen University, Guangzhou, China
   index: 4
 - name: Key Laboratory of Middle Atmosphere and Global Environment Observation, Institute of Atmospheric Physics, Chinese Academy of Sciences, Beijing, China
   index: 5
 - name: University of Chinese Academy of Sciences, Beijing, China
   index: 6

date: 2 April 2025
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

The Massive-Parallel Trajectory Calculations (MPTRAC) Lagrangian transport model [@hoffmann16; @hoffmann22] simulates the movement, dispersion, and chemical transformation of aerosol particles and trace gases in the atmosphere. Unlike grid-based models, MPTRAC tracks individual air parcels, offering high-resolution insights into transport, mixing, chemistry, and deposition. Trajectory calculations are driven by meteorological fields from data assimilation systems and forecast models. MPTRAC has been applied in various contexts, especially for tracking volcanic plumes in the troposphere and stratosphere. Optimized for high-performance computing (HPC) and graphics processing unit (GPU) systems, MPTRAC efficiently supports large-scale, high-resolution simulations, making it a powerful tool in atmospheric research.

# Statement of need

Accurate atmospheric transport modeling is essential for applications such as air quality assessments, climate studies, and public health protection. Traditional Eulerian models often struggle with fine-scale resolution, limiting their ability to capture complex transport and mixing processes. Lagrangian models, such as MPTRAC and others [@mckenna02; @lin03; @stohl05; @jones07; @stein15; @pisso19], address this limitation by tracking individual air parcels, offering fine-grained representations of atmospheric processes.

MPTRAC supports both research and operational applications requiring high-precision transport simulations. It models long-range transport and chemical transformations, enabling pollution studies, emission tracking, and environmental impact forecasting. MPTRAC runs efficiently on HPC and GPU systems, supporting fast, large-scale simulations that facilitate real-time decision-making and scientific analysis.

# Features

MPTRAC computes air parcel trajectories using wind and vertical velocity fields from global reanalysis or forecast data [@hersbach20; @gelaro17]. It accounts for eddy diffusion and subgrid-scale winds via the Langevin equation. An inter-parcel exchange module represents air mixing. Additionally, MPTRAC simulates convection, sedimentation, radioactive decay, gas and aqueous phase chemistry, and wet and dry deposition. Meteorological pre-processing routines compute boundary layer heights, CAPE, geopotential heights, and tropopause data. \autoref{fig:clusters} shows MPTRAC's geophysical and chemical modules and core software components.

![Geophysical and chemical modules and main software components of the MPTRAC Lagrangian transport model. Image adapted from @hoffmann24.\label{fig:clusters}](clusters.png){ width=95% }

Optimized for computational efficiency, MPTRAC features an MPI-OpenMP-OpenACC hybrid parallelization for scalable deployment on workstations, HPC systems, and GPU platforms [@liu20; @hoffmann22; @hoffmann24]. It supports multiple output formats for diverse data types, enabling visualization and analysis using Gnuplot, ParaView, or Python. MPTRAC is open-source, distributed under the GNU GPL v3 license. Code and documentation are available in a [GitHub repository](https://github.com/slcs-jsc/mptrac), with software releases archived on [Zenodo](https://doi.org/10.5281/zenodo.4400597).

# Applications

MPTRAC simulates the transport and dispersion of atmospheric constituents from both natural and anthropogenic sources. For volcanic eruptions, it helps estimate emissions and track the spread of volcanic ash and sulfate aerosols, which impact air traffic, climate, and ecosystems [@heng16; @wu17; @wu18; @cai22; @mishra22]. Similarly, MPTRAC models the dispersion of carbon dioxide from wildfires, including vertical ascent and widespread transport into the upper atmosphere [@liao24]. Further studies have examined the long-range transport of aerosol particles and trace gases in the upper troposphere and lower stratosphere [@smoydzin22; @wu23; @clemens24].

\autoref{fig:convection} illustrates MPTRAC's use in studying convective transport of air from the boundary layer into the free troposphere [@hoffmann23]. The simulation tracks air parcels lifted by updrafts from tropical storms and mid-latitude systems, revealing how boundary layer air and pollutants spread aloft, improving understanding of atmospheric circulation and pollutant distribution at regional and global scales.

![Lagrangian transport simulation of convective transport from the boundary layer into the free troposphere. One million trajectories are initialized near the surface on 1 July 2017, 00:00 UTC, and tracked over 10 days using ERA5 reanalysis data. Color coding indicates the geopotential height of the air parcels.\label{fig:convection}](convection.png){ width=95% }

# Evolution and Future Directions

MPTRAC has evolved since 2013, originally designed for HPC applications. Initially using OpenMP for multi-core CPUs, it later incorporated MPI for large-scale ensemble simulations. In 2019, OpenACC offloading boosted performance by porting all geophysical modules to NVIDIA GPUs. In 2023, we explored GPU offloading using OpenMP for AMD GPUs, broadening portability. MPTRAC is now ready for simulations on the Exascale [JUPITER](https://www.fz-juelich.de/en/ias/jsc/jupiter) system.

Recent technical efforts have improved documentation and usability. Written in C, MPTRAC includes a Fortran wrapper and high-level API, enabling seamless integration of Lagrangian transport simulations into other models. Continuous testing via GitHub Actions and multiple HPC systems, including the JUPITER Exascale Development Instrument ([JEDI](https://www.fz-juelich.de/en/ias/jsc/systems/supercomputers/jedi)), ensures robust performance and reliability.

# Acknowledgements

The development of MPTRAC is supported by the SDL Climate Science at JSC and the JL-ExaESM of the Helmholtz Association. We thank JSC for providing compute and storage resources. Funding was provided by BMBF (project ADAPTEX, FKZ 16ME0670) and DFG (project AeroTrac, HO 5102/1-1). AI tools supported language editing.

# References
