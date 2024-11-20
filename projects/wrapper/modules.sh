#!bin/bash

module --force purge
ml Stages/2024
ml GCC/12.3.0 ParaStationMPI/5.9.2-1

export LANG=C
export LC_ALL=C

# try other compilers...
#ml Intel/2023.2.1  ParaStationMPI/5.9.2-1
#ml NVHPC/23.7-CUDA-12  OpenMPI/4.1.5
#ml NVHPC/23.7-CUDA-12 ParaStationMPI/5.9.2-1

ml GCCcore/.12.3.0
ml GSL/2.7
ml netCDF/4.9.2 #(automatically loads HDF5)

