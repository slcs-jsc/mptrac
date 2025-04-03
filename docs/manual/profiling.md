# Profiling

## Introduction

In software engineering, profiling ("program profiling", "software profiling") is a form of dynamic program analysis that measures, for example, the space (memory) or time complexity of a program, the usage of particular instructions, or the frequency and duration of function calls. Most commonly, profiling information serves to aid program optimization.

Profiling is achieved by instrumenting either the program source code or its binary executable form using a tool called a profiler (or code profiler). Profilers may use a number of different techniques, such as event-based, statistical, instrumented, and simulation methods.

Profilers use a wide variety of techniques to collect data, including hardware interrupts, code instrumentation, instruction set simulation, operating system hooks, and performance counters. Profilers are used in the performance engineering process.

## Timers of the trac tool

A set of timers has been implemented directly in the trac tool. Next to the runtime for different parts of the code, also estimates of the memory needs and the problem size are reported.

Example of output reported by the trac tool:

```
SIZE_NP = 10000
SIZE_MPI_TASKS = 1
SIZE_OMP_THREADS = 4
MEMORY_ATM = 1449.58 MByte
MEMORY_CACHE = 648.499 MByte
MEMORY_CLIM = 15.3182 MByte
MEMORY_METEO = 7408.74 MByte
TIMER_ALLOC = 0.000 s    (min= 6.9767e-05 s, mean= 6.9767e-05 s, max= 6.9767e-05 s, n= 1)
TIMER_READ_CTL = 0.007 s    (min= 0.00717609 s, mean= 0.00717609 s, max= 0.00717609 s, n= 1)
TIMER_READ_CLIM = 0.025 s    (min= 0.0254529 s, mean= 0.0254529 s, max= 0.0254529 s, n= 1)
TIMER_READ_ATM = 0.028 s    (min= 0.0275215 s, mean= 0.0275215 s, max= 0.0275215 s, n= 1)
TIMER_MODULE_TIMESTEPS_INIT = 0.002 s    (min= 0.002085 s, mean= 0.002085 s, max= 0.002085 s, n= 1)
TIMER_GET_MET = 0.026 s    (min= 1.546e-06 s, mean= 1.76945e-05 s, max= 0.00997372 s, n= 1444)
TIMER_READ_MET_GRID = 0.003 s    (min= 0.000771007 s, mean= 0.000873196 s, max= 0.000986955 s, n= 4)
TIMER_READ_MET_SURFACE = 0.219 s    (min= 0.0519605 s, mean= 0.0547709 s, max= 0.0580699 s, n= 4)
TIMER_READ_MET_LEVELS = 5.490 s    (min= 1.30604 s, mean= 1.37255 s, max= 1.39654 s, n= 4)
TIMER_READ_MET_EXTRAPOLATE = 0.253 s    (min= 0.0239204 s, mean= 0.06324 s, max= 0.103068 s, n= 4)
TIMER_READ_MET_POLAR_WINDS = 0.006 s    (min= 0.00125896 s, mean= 0.00153137 s, max= 0.00172609 s, n= 4)
TIMER_READ_MET_PERIODIC = 0.002 s    (min= 0.000341735 s, mean= 0.000458583 s, max= 0.000582423 s, n= 4)
TIMER_READ_MET_GEOPOT = 1.594 s    (min= 0.394668 s, mean= 0.398482 s, max= 0.402873 s, n= 4)
TIMER_READ_MET_PV = 0.238 s    (min= 0.0565131 s, mean= 0.0594941 s, max= 0.0644101 s, n= 4)
TIMER_READ_MET_PBL = 0.122 s    (min= 0.0269698 s, mean= 0.0306072 s, max= 0.0331155 s, n= 4)
TIMER_READ_MET_TROPO = 0.737 s    (min= 0.183587 s, mean= 0.184353 s, max= 0.18602 s, n= 4)
TIMER_READ_MET_CLOUD = 0.085 s    (min= 0.0204375 s, mean= 0.0213186 s, max= 0.0221769 s, n= 4)
TIMER_READ_MET_CAPE = 2.116 s    (min= 0.517294 s, mean= 0.528923 s, max= 0.549389 s, n= 4)
TIMER_READ_MET_OZONE = 0.122 s    (min= 0.020632 s, mean= 0.0305517 s, max= 0.041007 s, n= 4)
TIMER_MODULE_CHEM_INIT = 0.000 s    (min= 0.000247314 s, mean= 0.000247314 s, max= 0.000247314 s, n= 1)
TIMER_MODULE_TIMESTEPS = 0.302 s    (min= 2.6811e-05 s, mean= 0.000209307 s, max= 0.000455049 s, n= 1441)
TIMER_MODULE_POSITION = 0.629 s    (min= 5.9016e-05 s, mean= 0.000218254 s, max= 0.000583772 s, n= 2882)
TIMER_MODULE_ADVECT = 1.685 s    (min= 0.000162452 s, mean= 0.00116961 s, max= 0.00311676 s, n= 1441)
TIMER_MODULE_DIFF_TURB = 1.742 s    (min= 0.000707676 s, mean= 0.00120864 s, max= 0.00243357 s, n= 1441)
TIMER_MODULE_DIFF_MESO = 1.544 s    (min= 0.000488041 s, mean= 0.00107166 s, max= 0.00232684 s, n= 1441)
TIMER_MODULE_CONVECTION = 0.962 s    (min= 0.000295857 s, mean= 0.000667748 s, max= 0.00142307 s, n= 1441)
TIMER_MODULE_METEO = 0.014 s    (min= 0.00228027 s, mean= 0.0035394 s, max= 0.00568389 s, n= 4)
TIMER_MODULE_BOUND_COND = 0.916 s    (min= 4.1641e-05 s, mean= 0.000317876 s, max= 0.000974128 s, n= 2882)
TIMER_MODULE_DECAY = 0.427 s    (min= 0.000157237 s, mean= 0.000296458 s, max= 0.000618609 s, n= 1441)
TIMER_MODULE_MIXING = 2.408 s    (min= 0.0298906 s, mean= 0.0329851 s, max= 0.0553491 s, n= 73)
TIMER_MODULE_CHEM_GRID = 1.744 s    (min= 0.000420581 s, mean= 0.00121059 s, max= 0.00318095 s, n= 1441)
TIMER_MODULE_OH_CHEM = 1.873 s    (min= 0.000114869 s, mean= 0.00129959 s, max= 0.00288309 s, n= 1441)
TIMER_MODULE_H2O2_CHEM = 0.776 s    (min= 5.8315e-05 s, mean= 0.00053884 s, max= 0.00215969 s, n= 1441)
TIMER_MODULE_TRACER_CHEM = 2.582 s    (min= 1.9956e-05 s, mean= 0.00179167 s, max= 0.0040358 s, n= 1441)
TIMER_MODULE_WET_DEPO = 0.532 s    (min= 1.9342e-05 s, mean= 0.000369044 s, max= 0.000838083 s, n= 1441)
TIMER_MODULE_DRY_DEPO = 0.453 s    (min= 2.7147e-05 s, mean= 0.000314659 s, max= 0.000882621 s, n= 1441)
TIMER_WRITE_ATM = 0.192 s    (min= 0.0456917 s, mean= 0.0480724 s, max= 0.049666 s, n= 4)
TIMER_WRITE_GRID = 0.328 s    (min= 0.0776646 s, mean= 0.0819251 s, max= 0.0879495 s, n= 4)
TIMER_WRITE_CSI = 1.025 s    (min= 0.000456732 s, mean= 0.000711579 s, max= 0.00215784 s, n= 1441)
TIMER_WRITE_ENS = 0.003 s    (min= 0.000652666 s, mean= 0.000704509 s, max= 0.000817902 s, n= 4)
TIMER_WRITE_PROF = 5.827 s    (min= 0.000380689 s, mean= 0.0040439 s, max= 0.0294752 s, n= 1441)
TIMER_WRITE_SAMPLE = 0.024 s    (min= 9.81003e-07 s, mean= 1.65024e-05 s, max= 0.0191332 s, n= 1441)
TIMER_WRITE_STATION = 0.620 s    (min= 0.000347288 s, mean= 0.000430191 s, max= 0.00285026 s, n= 1441)
TIMER_WRITE_VTK = 0.152 s    (min= 0.0374231 s, mean= 0.0378852 s, max= 0.0385151 s, n= 4)
TIMER_FREE = 0.098 s    (min= 0.0981576 s, mean= 0.0981576 s, max= 0.0981576 s, n= 1)
TIMER_GROUP_MEMORY = 0.098 s
TIMER_GROUP_INPUT = 5.798 s
TIMER_GROUP_PHYSICS = 18.593 s
TIMER_GROUP_METPROC = 5.276 s
TIMER_GROUP_OUTPUT = 8.171 s
TIMER_TOTAL = 37.936 s
```

SIZE_NP refers to the number of particles used in the simulation. SIZE_MPI_TASKS refers to the number of MPI tasks. SIZE_OMP_THREADS refers to the number of OpenMP threads per MPI task. SIZE_ACC_DEVICES refers to the number of GPU devices per MPI task.

The memory needs have been derived by analyzing the data structures used in the code and are not based on measurements.

Some timers are zero or not shown because the corresponding modules have not been used in this run.

A simple optimization of memory needs (and runtime!) is to adjust the constants EX, EY, and EZ defined in mptrac.h to match the grid dimensions of the meteorological data set and the constant NP to match the number of particles used in the simulation.

## Profiling of CPU runs

CPU profiling of MPTRAC is enabled by means of gprof. Set "PROF = 1" in the Makefile, recompile the code, and run a test case.

Runtime information will be collected in the file "gmon.out" in the working directory. Use "gprof [binary]" to see the runtime profile.

Example output for the trac tool:

```
Flat profile:
Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total           
 time   seconds   seconds    calls  ms/call  ms/call  name    
 27.09      2.85     2.85 21071464     0.00     0.00  intpol_met_space_3d
 15.49      4.48     1.63       96    16.98    17.61  write_output
  8.27      5.35     0.87     7770     0.11     0.68  frame_dummy
  7.13      6.10     0.75 11636847     0.00     0.00  intpol_met_time_3d
  5.04      6.63     0.53                             gomp_team_barrier_wait_end
  4.56      7.11     0.48                             cspline_eval
  3.33      7.46     0.35  2357631     0.00     0.00  clim_tropo
  3.33      7.81     0.35                             __cos_avx
  3.04      8.13     0.32                             __printf_fp_l
  2.57      8.40     0.27                             gsl_ran_gaussian_ziggurat
  2.28      8.64     0.24  3977062     0.00     0.00  intpol_met_space_2d
  2.28      8.88     0.24        3    80.00    80.00  spline
  2.09      9.10     0.22                             __mpn_divrem
  1.71      9.28     0.18                             gomp_barrier_wait_end
  1.62      9.45     0.17  7475552     0.00     0.00  locate_irr
  1.05      9.56     0.11  2658144     0.00     0.00  intpol_met_time_2d
  ...
```

## Profiling of GPU runs

GPU profiling is enabled by means of NVIDIA Nsight Systems. To see NVTX markers in the timeline, set "NVTX = 1" in the Makefile before the compilation.

Please find an example here on how to enable Nsight systems profilings in a job script:

```
    #### Nsight Sytems notes                                                                                                                                                                                            
    # For MPI profiling add "mpi" to the trace list:                                                                                                                                                                    
    #     --trace=nvtx,osrt,openacc,mpi                                                                                                                                                                                     
    # to start profiling after x sec add "--delay x":                                                                                                                                                                   
    #     nsys profile -o wodiff --delay 5 ...      starts after 5 sec                                                                                                                                               
    # To limit the profiling to x sec add "--duration x":                                                                                                                                                               
    #     nsys profile -o wodiff --delay 7 ...      limits the profiling to 7 sec                                                                                                                                    

    # Calculate trajectories...                                                                                                                                                                     
    nsys profile -f true -o wodiff --trace=nvtx,osrt,openacc --stats=true \
        $trac/trac data/dirlist trac.ctl atm_split.tab meteo/ei ATM_BASENAME atm | tee data/log_trac.txt
```

After the execution, you can visualize .qdrep files in nsys-ui. Take a look at "Nsight Sytems notes" for other useful options.

## Further reading

- [gprof manual](http://sourceware.org/binutils/docs/gprof)

- [Nsight Systems documentation](https://developer.nvidia.com/nsight-systems)
