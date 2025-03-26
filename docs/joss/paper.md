---
title: 'MPTRAC: A high-performance Lagrangian transport model for atmospheric air parcel dispersion'
tags:
  - Atmospheric modeling
  - Lagrangian transport
  - Particle dispersion
  - Troposphere
  - Stratosphere
  - Aerosol transport
  - High-performance computing (HPC)
  - GPU acceleration
  - Numerical simulation
  - Advection and diffusion
  - Environmental monitoring
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
    orcid: TODO
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
    orcid: TODO
  - name: Xue Wu
    affiliation: "5, 6"
    orcid: TODO
  - name: Ling Zou
    affiliation: "1, 2"
    orcid: 0000-0001-6563-9815
affiliations:
 - name: Jülich Supercomputing Centre, Forschungszentrum Jülich, Jülich, Germany
   index: 1
 - name: Centre for Advanced Simulation and Analytics, Forschungszentrum Jülich, Jülich, Germany
   index: 2
 - name: Institute of Climate and Energy, Forschungszentrum Jülich, Jülich, Germany
   index: 3
 - name: School of Computer Science and Engineering, Sun Yat-sen University, Guangzhou, China
   index: 4
 - name: Key Laboratory of Middle Atmosphere and Global Environment Observation, Institute of Atmospheric Physics, Chinese Academy of Sciences, Beijing, China
   index: 5
 - name: University of Chinese Academy of Sciences, Beijing, China
   index: 6

date: 26 March 2025
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

The Massive-Parallel Trajectory Calculations (MPTRAC) Lagrangian transport model simulates the movement, dispersion, and chemical transformation of air parcels carrying pollutants, aerosols, and trace gases. Unlike Eulerian models, which use a fixed grid, MPTRAC tracks trajectories of individual air parcels, offering a high-resolution view of atmospheric transport and chemistry. Optimized for high-performance computing, MPTRAC combines advanced numerical techniques with computational efficiency, making it a powerful tool for studying dispersion, long-range transport, deposition, and trace gas reactions under diverse meteorological conditions.

# Statement of need

Accurate pollutant transport modeling is crucial for assessing air quality, climate impacts, and public health. While traditional Eulerian models often lack fine-scale detail, Lagrangian models track individual air parcels for a more precise representation of dispersion, deposition, and atmospheric interactions. MPTRAC has been successfully applied to studying pollutant transport in the free troposphere and stratosphere, proving its value for atmospheric research, environmental monitoring, and emergency response. Unlike many Lagrangian models, MPTRAC is optimized for high-performance computing (HPC), enabling efficient large-scale, high-resolution simulations.

# Scientific background

Lagrangian transport models track individual air parcels to simulate the long-range transport, dispersion, and chemical transformation of aerosols and trace gases [@mckenna02; @lin03; @stohl05; @jones07; @stein15; @pisso19]. Unlike Eulerian models, which use a fixed grid, Lagrangian models provide greater flexibility, making them ideal for studying dispersion, mixing, and chemical reactions. They solve the equations of motion for air parcels, accounting for advection, diffusion, settling, and deposition, and can simulate atmospheric chemistry. Driven by meteorological data from weather prediction systems, they are used for applications such as volcanic ash tracking, wildfire smoke dispersion, and monitoring nuclear or industrial accidents. Advances in turbulence modeling, chemistry, and computational efficiency continue to improve their accuracy and expand their use in environmental monitoring and emergency response.

# Features

MPTRAC [@hoffmann16; @hoffmann22] is a Lagrangian particle dispersion model for analyzing atmospheric transport in the free troposphere and stratosphere. It computes air parcel trajectories using wind and vertical velocity fields from global reanalysis or forecast data, such as ECMWF's ERA5 [@hersbach20] or NASA's MERRA-2 [@gelaro17], while accounting for mesoscale diffusion and wind fluctuations with the Langevin equation. An inter-parcel exchange module represents air mixing. The model also simulates convection, sedimentation, decay, gas and aqueous phase chemistry, and wet/dry deposition, with meteorological pre-processing for boundary layer height, convective available potential energy, geopotential heights, and tropopause data.

Optimized for computational efficiency, MPTRAC features an MPI-OpenMP-OpenACC hybrid parallelization for scalable deployment on workstations, HPC systems, and GPU platforms [@liu20; @hoffmann22; @hoffmann24]. It supports multiple output formats, including particle, grid, ensemble, profile, sample, and station data, with visualization via Gnuplot and ParaView. MPTRAC is open-source software distributed under the GNU GPL v3 license. \autoref{fig:clusters} illustrates MPTRAC's geophysical modules and core software components.

![Geophysical modules and main software components of the MPTRAC model. Image adapted from @hoffmann24.\label{fig:clusters}](clusters.png){ width=80% }

# Applications

MPTRAC simulates the transport and dispersion of aerosols and trace gases from both natural and anthropogenic sources. For volcanic eruptions, it tracks plumes carried by wind, helping estimate the emissions and predict the spread of volcanic ash and sulfate aerosols, which impact air travel, climate, and ecosystems [@heng16; @wu17; @wu18; @cai22]. Similarly, MPTRAC models the dispersion of carbon dioxide and smoke from wildfires, including their ascent into the upper atmosphere and eventual removal via deposition or precipitation [@liao24]. Other studies with MPTRAC addressed the transport of aerosols and trace gases in the upper troposphere and lower stratosphere region [@smoydzin22; @wu23; @clemens24]. These applications are vital for air quality monitoring, hazard assessment, and environmental impact studies.

An example, shown in \autoref{fig:convection}, illustrates MPTRAC’s use in studying the transport of air from the planetary boundary layer (PBL) into the free troposphere. MPTRAC tracks how air parcels lifted by updrafts or convective currents move into the more stable free troposphere, where pollutants or aerosols can spread over large areas. This helps provide insights into atmospheric circulation and pollutant distribution on regional to global scales.

![Lagrangian transport simulation of convective transport from the planetary boundary layer (PBL) into the free troposphere. One million trajectories are initialized in the PBL at 00:00 UTC on 1 July 2017 and tracked over 10 days using ERA5 reanalysis data.\label{fig:convection}](convection.png)

# Evolution and Future Directions

MPTRAC development began in 2013, designed from the ground up for HPC applications. Initially using OpenMP for multi-core CPUs, it later incorporated MPI for large-scale ensemble simulations. In 2019, OpenACC offloading enabled execution on NVIDIA GPUs, significantly boosting performance. Over time, all geophysical modules were ported to GPUs, optimizing efficiency and minimizing data transfers. Recent efforts have improved documentation, enhanced usability, and integrated continuous testing across GitHub Actions and multiple HPC systems, including the JUPITER Exascale machine development platform.

Future development will expand MPTRAC’s applications, particularly in the planetary boundary layer (PBL), where turbulence and surface interactions affect transport. Planned enhancements include terrain-following coordinates for better airflow representation over complex topography and advanced turbulence parametrizations for more accurate small-scale transport modeling. These improvements will strengthen MPTRAC’s role in studying air quality, pollutant dispersion, and other environmental challenges.

# Acknowledgements

The development of MPTRAC has been continuously supported by the Simulation and Data Laboratory Climate Science at the Jülich Supercomputing Centre (JSC) and the Joint Lab Exascale Earth System Modeling of the Helmholtz Association (HGF). We acknowledge JSC for providing essential compute time and storage resources. Financial support for scientific applications and software development has been provided by the DFG project AeroTrac (HO 5102/1-1) and the BMBF project ADAPTEX (FKZ 16ME0670). AI-assisted tools, including ChatGPT and DeepL, were used for manuscript preparation and language editing.

# References
