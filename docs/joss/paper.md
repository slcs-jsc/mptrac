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
    affiliation: "5"
    orcid: TBD
  - name: Ling Zou
    affiliation: "1, 2"
    orcid: 0000-0001-6563-9815
affiliations:
 - name: Jülich Supercomputing Centre, Forschungszentrum Jülich, Germany
   index: 1
 - name: Centre of Advanced Simulation and Analytics, Forschungszentrum Jülich, Germany
   index: 2
 - name: Institute of Climate and Energy, Forschungszentrum Jülich, Germany
   index: 3
 - name: Institute of..., Sun Yat-sen University, China
   index: 4
 - name: IAP Beijing, China
   index: 5
date: 25 March 2025
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

The Massive-Parallel Trajectory Calculations (MPTRAC) Lagrangian transport model simulates the movement, dispersion, and chemical transformation of air parcels carrying pollutants, aerosols, and trace gases. Unlike Eulerian models, which operate on a fixed spatial grid, MPTRAC tracks individual air parcels along their trajectories, providing a high-resolution view of atmospheric transport and chemical reactions. It enables detailed studies of dispersion, long-range transport, deposition, and trace gas chemistry. Optimized for high-performance computing, MPTRAC combines advanced numerical techniques with computational efficiency, making it a versatile tool for simulating atmospheric transport and chemistry under diverse meteorological conditions.

# Statement of need

Accurate pollutant transport modeling is essential for assessing air quality, climate impacts, and public health. While traditional Eulerian models often lack fine-scale detail, Lagrangian models track individual air parcels, offering a more precise representation of dispersion, deposition, and atmospheric interactions. MPTRAC has been successfully used to study pollutant transport from volcanic eruptions and wildfires, proving its value in environmental monitoring, climate research, and emergency response. Unlike many existing Lagrangian models, MPTRAC is specifically optimized for high-performance computing (HPC), enabling efficient large-scale, high-resolution simulations.

# Scientific background

Lagrangian transport models track individual air parcels to simulate the long-range transport, dispersion, and chemical transformation of aerosols and trace gases [@mckenna02; @lin03; @stohl05; @jones07; @stein15; @pisso19]. Unlike Eulerian models, which use a fixed grid, Lagrangian models provide greater flexibility, making them ideal for studying dispersion, mixing, and chemical reactions. They solve the equations of motion for air parcels, accounting for advection, diffusion, settling, and deposition, and can simulate atmospheric chemistry. Driven by meteorological data from weather prediction systems, they are used for applications such as volcanic ash tracking, wildfire smoke dispersion, and monitoring nuclear or industrial accidents. Advances in turbulence modeling, chemistry, and computational efficiency continue to improve their accuracy and expand their use in environmental monitoring and emergency response.

# Features

MPTRAC [@hoffmann16; @hoffmann22] is a Lagrangian particle dispersion model designed to analyze atmospheric transport in the free troposphere and stratosphere. It computes air parcel trajectories using wind and vertical velocity fields from global reanalysis or forecast data, such as ECMWF ERA5 [@hersbach20]. The model accounts for mesoscale diffusion and wind fluctuations with the Langevin equation and includes an inter-parcel exchange module for air mixing.

MPTRAC also simulates convection, sedimentation, decay, gas and aqueous phase chemistry, and wet/dry deposition. It provides meteorological data processing for boundary layers, convective potential energy, geopotential heights, and tropopause data. The model supports multiple output methods, including particle, grid, ensemble, profile, sample, and station data, with visualization via Gnuplot and ParaView.

Optimized for efficiency, MPTRAC uses MPI-OpenMP-OpenACC hybrid parallelization, making it suitable for deployment on workstations, HPC, and GPU platforms. It is available as open-source software under the GNU GPL license.

\autoref{fig:clusters} illustrates MPTRAC's main software components.

![Main software components of the MPTRAC model. Image adapted from TODO.\label{fig:clusters}](clusters.png){ width=50% }

# Applications

The MPTRAC model has important applications in simulating the transport and dispersion of aerosols and trace gases from natural or anthropogenic sources. During a volcanic eruption, ash plumes are ejected high into the atmosphere, and the model can track the movement of these particles as they are carried by wind currents across vast distances. MPTRAC helps in predicting the spread of volcanic ash and sulfate aerosol, which can have significant impacts on air travel, climate, and ecosystems. Similarly, the model has been used to simulate the transport of carbon dioxide from wildfires. As wildfires release large amounts of smoke and aerosols into the atmosphere, MPTRAC can simulate how these particles are dispersed, including their ascent into the upper atmosphere and eventual removal through deposition or precipitation. These applications are crucial for air quality monitoring, hazard assessment, and environmental impact studies related to both volcanic activity and wildfire events.

An example application of MPTRAC, illustrated in \autoref{fig:convection}, is in studying the transport of air from the planetary boundary layer (PBL) into the free troposphere. The PBL, the lowest layer of the atmosphere, is characterized by turbulent mixing and direct interaction with the Earth's surface, which often leads to the accumulation of pollutants or aerosols. Using MPTRAC, we can track how air parcels, once lifted by strong updrafts or convective currents, move into the more stable free troposphere. Once in the free troposphere, the air parcels are subject to less turbulent mixing, and the pollutants or aerosols they carry can be dispersed over large areas within a matter of days. MPTRAC simulates how these air parcels are transported, spread, and ultimately removed from the atmosphere through processes like deposition or precipitation, providing crucial insights into atmospheric circulation and pollutant distribution on regional to global scales.

![Lagrangian transport simulation showing convective transport from the planetary boundary layer (PBL) into the free troposphere. Trajectories are initialized in the PBL at 00:00 UTC on 1 July 2017 and are tracked over a 10-day period. The simulation is driven by ERA5 meteorological reanalysis data.\label{fig:convection}](convection.png)

Evolution and Future Directions of MPTRAC

MPTRAC development began in 2013, designed from the ground up for high-performance computing (HPC). Initially using OpenMP for multi-core CPUs, it later incorporated MPI for large-scale ensemble simulations. In 2019, OpenACC offloading enabled execution on NVIDIA GPUs, significantly boosting performance. Over time, all geophysical modules were ported to GPUs, optimizing efficiency and minimizing data transfers. Recent efforts have improved documentation, enhanced usability, and integrated continuous testing across GitHub Actions and multiple HPC systems, including the JUPITER Exascale platform.

Future development will expand MPTRAC’s applications, particularly in the planetary boundary layer (PBL), where turbulence and surface interactions affect transport. Planned enhancements include terrain-following coordinates for better airflow representation over complex topography and advanced turbulence parametrizations for more accurate small-scale transport modeling. These improvements will strengthen MPTRAC’s role in studying air quality, pollutant dispersion, and other environmental challenges.

# Acknowledgements

The development of MPTRAC has been continuously supported by the Simulation and Data Laboratory Climate Science at the Jülich Supercomputing Centre (JSC) and the Joint Lab Exascale Earth System Modeling of the Helmholtz Association (HGF). We acknowledge JSC for providing essential compute time and storage resources. Financial support for scientific applications and software development has been provided by the DFG project AeroTrac (HO 5102/1-1) and the BMBF project ADAPTEX (FKZ 16ME0670). AI-assisted tools, including ChatGPT and DeepL, were used for manuscript preparation and language editing.

# References
