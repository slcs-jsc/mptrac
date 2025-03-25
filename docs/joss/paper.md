---
title: 'MPTRAC: A high-performance Lagrangian transport model for atmospheric particle dispersion'
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
 - name: Jlich Supercomputing Centre, Forschungszentrum Jlich, Germany
   index: 1
 - name: Centre of Advanced Simulation and Analytics, Forschungszentrum Jlich, Germany
   index: 2
 - name: Institute of Climate and Energy, Forschungszentrum Jlich, Germany
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

The MPTRAC Lagrangian transport model is a numerical tool used to simulate the transport of particles, such as pollutants or aerosols, in the atmosphere. It operates on the Lagrangian framework, which follows individual particles as they move through the atmosphere, instead of tracking the air mass itself (as in Eulerian models). This model is particularly useful for studying the dispersion, transport, and deposition of substances, allowing researchers to analyze how pollutants travel over time and space. It helps in understanding atmospheric processes such as air quality, climate modeling, and the spread of hazardous materials. The model is often employed in environmental studies, weather forecasting, and assessing the impact of pollutants on different regions.

# Statement of need

As global environmental concerns grow, particularly regarding air quality and the dispersion of pollutants, there is an increasing need for advanced tools to accurately model the transport and behavior of particles in the atmosphere. Understanding how pollutants, such as aerosols, greenhouse gases, and hazardous materials, move through air currents is crucial for assessing their impact on climate, human health, and ecosystems. Current models often rely on Eulerian frameworks that track the movement of air masses, but they fall short in capturing the fine details of particle transport.

The MPTRAC Lagrangian transport model addresses this gap by simulating the movement of individual particles, offering a more precise representation of how pollutants disperse, deposit, and interact with various atmospheric conditions over time. This model is vital for improving air quality predictions, informing climate change studies, and assisting in emergency response planning for chemical spills or other atmospheric hazards. There is a clear need for a reliable and efficient tool like MPTRAC to enhance environmental monitoring, policy-making, and public health protection efforts.

# Background

The MPTRAC Lagrangian transport model simulates the movement of particles in the atmosphere by following their individual trajectories over time, using a Lagrangian approach. Unlike Eulerian models, which focus on fixed locations, this method tracks the particles as they are transported by wind, allowing for a more detailed understanding of their dispersion. The primary forces acting on the particles include atmospheric advection, where the wind moves particles horizontally and vertically, and turbulent diffusion, which accounts for random small-scale air movements that spread particles across larger areas.

The model also considers gravitational settling, where heavier particles fall to the ground more quickly, and both dry and wet deposition processes. Dry deposition involves the direct settling of particles on surfaces, while wet deposition occurs when particles are removed by precipitation. Additionally, MPTRAC incorporates atmospheric fields like wind patterns, temperature, and pressure, as well as the boundary layer's influence, which affects how particles behave near the surface. By solving a set of differential equations that account for these physical processes, the model provides accurate predictions of pollutant and aerosol transport, making it a valuable tool for environmental monitoring and air quality forecasting.

# Features

The MPTRAC (Massive-Parallel Trajectory Calculations) model is a Lagrangian particle dispersion tool designed to analyze atmospheric transport processes in the free troposphere and stratosphere. It computes air parcel trajectories by solving the kinematic equation of motion, utilizing horizontal wind and vertical velocity fields from global reanalyses or forecasts. To account for mesoscale diffusion and subgrid-scale wind fluctuations, MPTRAC employs the Langevin equation, introducing stochastic perturbations to the trajectories. Additionally, it features an inter-parcel exchange module to represent air mixing.​
GitHub+1GitHub+1

The model includes modules for simulating convection, sedimentation, exponential decay, gas and aqueous phase chemistry, and both wet and dry deposition processes. Its meteorological data pre-processing capabilities provide estimates for boundary layers, convective available potential energy, geopotential heights, potential vorticity, and tropopause data. MPTRAC offers various output methods, including particle, grid, ensemble, profile, sample, and station data, and supports visualization through Gnuplot and ParaView interfaces. Designed for efficiency, it utilizes MPI-OpenMP-OpenACC hybrid parallelization and specific code optimizations, making it suitable for deployment on systems ranging from single workstations to high-performance computing (HPC) and GPU platforms. The model is distributed as open-source software under the GNU GPL license.

Figure \autoref{fig:clusters} shows the main software components of MPTRAC.

![Main software coponents of the MPTRAC model. Image adapted from TODO.\label{fig:clusters}](clusters.png)

# Applications

The MPTRAC model has important applications in simulating the transport and dispersion of particles from natural events like volcanic eruptions and wildfires. During a volcanic eruption, ash plumes are ejected high into the atmosphere, and the model can track the movement of these particles as they are carried by wind currents across vast distances. MPTRAC helps in predicting the spread of volcanic ash, which can have significant impacts on air travel, climate, and ecosystems. Similarly, the model is used to simulate the transport of smoke and particulate matter from wildfires. As wildfires release large amounts of smoke and aerosols into the atmosphere, MPTRAC can simulate how these particles are dispersed, including their ascent into the upper atmosphere and eventual removal through deposition or precipitation. These applications are crucial for air quality monitoring, hazard assessment, and environmental impact studies related to both volcanic activity and wildfire events.

An example application of the MPTRAC Lagrangian transport model is in studying the transport of air parcels from the planetary boundary layer (PBL) into the free troposphere. The PBL, the lowest layer of the atmosphere, is characterized by turbulent mixing and direct interaction with the Earth's surface, which often leads to the accumulation of pollutants or aerosols. Using MPTRAC, researchers can track how these air parcels, once lifted by strong updrafts or convective currents, move into the more stable free troposphere. Once in the free troposphere, the air parcels are subject to less turbulent mixing, and the pollutants or aerosols they carry can be dispersed over large areas within a matter of days. The model simulates how these particles are transported, spread, and ultimately removed from the atmosphere through processes like deposition or precipitation, providing crucial insights into atmospheric circulation and pollutant distribution on regional to global scales.

# Development history and outlook

Development of the MPTRAC model began in 2013, with the code built from scratch specifically for high-performance computing (HPC) applications. Initially, the model utilized OpenMP parallelization, allowing it to efficiently run on multi-core CPU systems. In 2019, OpenACC offloading was introduced to enable execution on NVIDIA GPUs, significantly enhancing performance. Over time, all geophysical modules were ported to GPUs, ensuring that trajectory calculations are now fully offloaded, minimizing host-device data transfers and optimizing efficiency. In recent years, substantial efforts have been made to improve software documentation for both users and developers, as well as to implement continuous integration and deployment (CI/CD) methods for thorough code testing. Today, MPTRAC undergoes routine testing in GitHub Actions and is validated on three different HPC systems at JSC, including the development platform for the upcoming JUPITER Exascale machine, ensuring its reliability and scalability for cutting-edge atmospheric transport simulations.

While the HPC capabilities of MPTRAC are already highly advanced, future development will focus on expanding its applicability to a broader range of atmospheric transport simulations. A key objective is to enhance the realism of simulations in the planetary boundary layer (PBL), where turbulence and surface interactions play a crucial role in particle transport. To achieve this, we plan to implement the direct use of terrain-following coordinates from meteorological input data, improving the model’s representation of airflow over complex topography. Additionally, more advanced parametrizations of turbulence and diffusion in the PBL will be introduced, allowing for a more accurate depiction of small-scale transport processes. These improvements will extend MPTRAC’s capabilities beyond the free troposphere and stratosphere, making it an even more powerful tool for studying near-surface air quality, pollutant dispersion, and other environmental applications.

# Acknowledgements

The development of MPTRAC is continuously supported by the Simulation and Data Laboratory Climate Science, located at the Jülich Supercomputing Centre (JSC), Germany, and the Joint Lab Exascale Earth System Modeling of the Helmholtz Association (HGF). We gratefully acknowledge JSC for providing compute time and storage resources, which have been essential for the advancement of MPTRAC. Financial support for scientific applications and software development has been provided by the DFG project AeroTrac (grant #TODO) and the BMBF project ADAPTEX (grant #TODO). We sincerely thank TODO for their valuable comments and suggestions on this manuscript. Additionally, we acknowledge the use of AI-assisted tools, including ChatGPT and DeepL, for manuscript preparation and language editing.

# References
