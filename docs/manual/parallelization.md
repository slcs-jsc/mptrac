# Parallelization

## Introduction

MPTRAC features an MPI-OpenMP-OpenACC hybrid parallelization for efficient use on supercomputers and GPUs. An application built with the hybrid model of parallel programming can run on a compute cluster using both OpenMP and MPI, such that OpenMP is used for parallelism within a (multi-core) node while MPI is used for parallelism between nodes. In combination with MPI, OpenACC can be used to conduct multi-GPU simulations.

MPI (Message Passing Interface) is a communication protocol for programming parallel computers. Both point-to-point and collective communication are supported. MPI is a message-passing application programmer interface, together with protocol and semantic specifications for how its features must behave in any implementation. MPI's goals are high performance, scalability, and portability. MPI remains the dominant model used in high-performance computing today.

The application programming interface OpenMP (Open Multi-Processing) supports multi-platform shared-memory multiprocessing programming in C, C++, and Fortran, on many platforms, instruction-set architectures and operating systems. It consists of a set of compiler directives, library routines, and environment variables that influence run-time behavior. OpenMP uses a portable, scalable model that gives programmers a simple and flexible interface for developing parallel applications for platforms ranging from the standard desktop computer to the supercomputer.

OpenACC (for open accelerators) is a programming standard for parallel computing developed by Cray, CAPS, Nvidia and PGI. The standard is designed to simplify parallel programming of heterogeneous CPU/GPU systems. As in OpenMP, the programmer can annotate C, C++ and Fortran source code to identify the areas that should be accelerated using compiler directives and additional functions.

## Parallelization strategy of MPTRAC

In MPTRAC, the OpenMP parallelization is always enabled to exploit the multi-core compute capabilities of modern computing systems. The OpenMP parallelization is used to distribute the trajectory calculations for the air parcels of a single simulation over the cores of a single compute node. The calculations share the same meteorological data, i.e., only a single copy of this potentially large data set needs to be kept in memory. As the trajectory calculations typically need most of the total runtime of the simulations but can be conducted independently of each other, this parallelization scales very effectively over a large range of compute cores. The OpenMP scalability of MPTRAC has been discussed in the study of [Rößler et al. [2018]](https://doi.org/10.5194/gmd-11-575-201).

The MPI parallelization of MPTRAC can be used to distribute a set of independent simulations (i.e., an ensemble of simulations) over a range of compute nodes. The simulations will run independently of each other, with their own individual input data and control parameters. A list of directories providing the input data for the different simulations needs to be provided as a control parameter to the trac tool. The trac tool will distribute this list of simulations over the different compute nodes. On each compute node, the individual simulations will be parallelized with OpenMP to make use of the different compute cores of each node. Also, mixed MPI/OpenMPI setups are possible. For example, on a single compute node consisting of 48 compute cores, 4 independent simulations can be conducted by means of MPI, where each simulation employes 12 compute cores of the node by means of OpenMP.

If the compute nodes are equipped with GPUS, the MPTRAC can make use of these by means of then OpenACC parallelization. In particular, the compute-intensive parts of the trajectory calculations will be offloaded to the GPUs. As modern computing systems typically provide more than GPU per compute node, the MPI parallelization of MPTRAC needs to be used to exploit all the GPUs of a node for the calculations. For example, for compute nodes equipped with 4 GPUs, at least 4 individual simulations need to be performed to make use of all GPUs. At the same time, the OpenMP parallelization should be used to make efficient use of the CPUs. If the compute nodes have 48 cores, then 12 cores per simulation should be employed by means of OpenMP.

Various combinations of an MPI/OpenMP/OpenACC hybrid parallelization can be configured. The best parallelization setup will depend on the requirements for the individual simulations and the target computing system.

## Supercomputers in Jülich

The JUWELS cluster module is composed of more than 2000 standard compute nodes. Each compute node consists of 48 compute cores and is equipped with 96 GByte of memory. Although there are 48 compute cores physically, a larger number of cores, e.g., 96 cores, can be used by means of the hyper-threading approach in the OpenMP parallelization. Hyper-treading can help to use the compute cores more effectively.

The JUWELS booster module is composed of 936 compute nodes with 48 cores and 512 GByte per node. Each compute node is equipped with 4 NVIDIA A100 Tensor Core-GPUs, which can be used by means of the OpenACC parallelization.

Examples of how to set up the job scripts for different types of MPI/OpenMP/OpenACC simulations can be found in [JUWELS Quick Introduction](https://apps.fz-juelich.de/jsc/hps/juwels/quickintro.html).

Alternatively, MPTRAC has also been used on the [JURECA-DC system](https://apps.fz-juelich.de/jsc/hps/jureca/quickintro.html).

## Further Reading

- [MPI Forum](https://www.mpi-forum.org/)
- [OpenMP](https://www.openmp.org/)
- [OpenACC](https://www.openacc.org/)
