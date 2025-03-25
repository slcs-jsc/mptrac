---
title: 'MPTRAC: A Lagrangian transport model'
tags:
  - C
  - atmosphere
  - troposphere
  - stratosphere
  - transport
  - particle methods
  - GPU
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
  - name: Farahnaz Khosrawi
    affiliation: "1, 2"
    orcid: 0000-0002-0261-7253
  - name: Mingzhao Liu
    affiliation: "4"
    orcid: 0000-0001-7314-9568
  - name: Catrin Meyer
    affiliation: "1, 2"
    orcid: 0000-0002-9271-6174
  - name: Yen-Sen Lu
    affiliation: "1, 2"
    orcid: 0000-0002-3255-824X
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

Figure \autoref{fig:clusters} shows the sturcture and main features of the software.


Figures can be included like this:
![Caption for example figure.\label{fig:clusters}](clusters.png) { width=60% }
and referenced from text using .


# Applications

The MPTRAC model has important applications in simulating the transport and dispersion of particles from natural events like volcanic eruptions and wildfires. During a volcanic eruption, ash plumes are ejected high into the atmosphere, and the model can track the movement of these particles as they are carried by wind currents across vast distances. MPTRAC helps in predicting the spread of volcanic ash, which can have significant impacts on air travel, climate, and ecosystems. Similarly, the model is used to simulate the transport of smoke and particulate matter from wildfires. As wildfires release large amounts of smoke and aerosols into the atmosphere, MPTRAC can simulate how these particles are dispersed, including their ascent into the upper atmosphere and eventual removal through deposition or precipitation. These applications are crucial for air quality monitoring, hazard assessment, and environmental impact studies related to both volcanic activity and wildfire events.

# Example application

An example application of the MPTRAC Lagrangian transport model is in studying the transport of air parcels from the planetary boundary layer (PBL) into the free troposphere. The PBL, the lowest layer of the atmosphere, is characterized by turbulent mixing and direct interaction with the Earth's surface, which often leads to the accumulation of pollutants or aerosols. Using MPTRAC, researchers can track how these air parcels, once lifted by strong updrafts or convective currents, move into the more stable free troposphere. Once in the free troposphere, the air parcels are subject to less turbulent mixing, and the pollutants or aerosols they carry can be dispersed over large areas within a matter of days. The model simulates how these particles are transported, spread, and ultimately removed from the atmosphere through processes like deposition or precipitation, providing crucial insights into atmospheric circulation and pollutant distribution on regional to global scales.


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Acknowledgements

We acknowledge contributions from ...

# References
