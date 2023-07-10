# atm_dist

This applocation calculates the transport deviations of trajectories.

atm_dist can be called by:
```
# calling atm_dist
$ ./atm_dist  <ctl> <dist.tab> <param> <atm1a> <atm1b> [<atm2a> <atm2b> ...]
```
The required arguments to run atm_dist are as follows:
* ctl: is the control (.ctl) file of the trajectory simulation  
* dist.tab: the output file where the deviations should be written into  
* param: the deviation paramters to be calculated (e.g. abs_dev etc. see below)
* atm1a and atm1b: the trajectory files that should be compared

Statistical paramters that can be calculated are:

* absdev
* mean
* stddv
* min
* max
* median
* mad

The here used statistics are the gnu statistic library. More informati on the [gnu statistics library](https://www.gnu.org/software/gsl/doc/html/statistics.html) can be found on their webpage. 

An example command line for running atm_dist looks like as follows:

```
./atm_dist trac.ctl dist_absdev.tab abs_dev traj_sim1_2017_01_01_00.tab traj_sim2_2017_01_01_00.tab
```

Note: The command line needs to be adjusted to the respective directories were the files (or executables in case of atm_dist) are located. 
Currently, transport deviations with atm_dist can (for the time being) only be caluculated for the ASCCI output files (.tab files) 

Additional paramters that can be used when running atm_dist are:

DIST_LAT0 -89 DIST_LAT1 89  - The latitude range to be considered    
DIST_Z0 8 DIST_Z1 16 - The altitude range to be considered   
DISTZCORE - can be used for filtering outliers  
