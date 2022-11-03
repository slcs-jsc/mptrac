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
SIZE_ACC_DEVICES = 0
MEMORY_ATM = 1220.7 MByte
MEMORY_CACHE = 1885.26 MByte
MEMORY_METEO = 6812.07 MByte
MEMORY_DYNAMIC = 4019.6 MByte
MEMORY_STATIC = 1307.61 MByte
TIMER_INIT = 1.134 s
TIMER_INPUT = 0.001 s
TIMER_OUTPUT = 0.852 s
TIMER_ADVECT = 0.705 s
TIMER_DECAY = 0.141 s
TIMER_DIFFMESO = 0.275 s
TIMER_DIFFTURB = 0.322 s
TIMER_ISOSURF = 0.000 s
TIMER_METEO = 1.140 s
TIMER_POSITION = 0.050 s
TIMER_SEDI = 0.000 s
TIMER_OHCHEM = 0.000 s
TIMER_WETDEPO = 0.000 s
TIMER_TOTAL = 4.639 s
```

SIZE_NP refers to the number of particles used in the simulation. SIZE_MPI_TASKS refers to the number of MPI tasks. SIZE_OMP_THREADS refers to the number of OpenMP threads per MPI task. SIZE_ACC_DEVICES refers to the number of GPU devices per MPI task.

The memory needs have been derived by analyzing the data structures used in the code and are not based on measurements.

Some timers are zero because the corresponding modules have not been used in this run.

A simple optimization of memory needs (and runtime!) is to adjust the constants EX, EY, and EZ defined in libtrac.h to match the grid dimensions of the meteorological data set and the constant NP to match the number of particles used in the simulation.

## Profiling of CPU runs

CPU profiling of MPTRAC is enabled by means of gprof. Uncomment "PROF = 1" in the Makefile, recompile the code, and run a test case.

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

GPU profiling is enabled by means of NVIDIA Nsight Systems. To see NVTX markers in the timeline, uncomment "USE_NVTX = 1" in the Makefile before the compilation.

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

https://en.wikipedia.org/wiki/Profiling_(computer_programming)

https://developer.nvidia.com/nsight-systems

http://sourceware.org/binutils/docs/gprof
