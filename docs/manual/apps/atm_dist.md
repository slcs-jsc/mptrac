# atm_dist

This application calculates the transport deviations of trajectories.

atm_dist can be called by:
```
# calling atm_dist
$ ./atm_dist  <ctl> <dist.tab> <param> <atm1a> <atm1b> [<atm2a> <atm2b> ...]
```
The required arguments to run atm_dist are as follows:
* ctl: is the control (.ctl) file of the trajectory simulation  
* dist.tab: the output file where the deviations should be written into  
* param: the deviation parameters to be calculated, for example `absdev`.
* atm1a and atm1b: the trajectory files that should be compared

Statistical parameters that can be calculated are:

* mean
* stddev
* min
* max
* skew
* kurt
* absdev
* median
* mad

The statistics are calculated with the GNU Scientific Library. More information on the [GNU Scientific Library statistics functions](https://www.gnu.org/software/gsl/doc/html/statistics.html) can be found on their webpage.

An example command line for running atm_dist looks like as follows:

```
./atm_dist trac.ctl dist_absdev.tab absdev traj_sim1_2017_01_01_00.tab traj_sim2_2017_01_01_00.tab
```

Note: The command line needs to be adjusted to the respective directories were the files (or executables in case of atm_dist) are located. 

Additional parameters that can be used when running atm_dist are:

DIST_LAT0 -89 DIST_LAT1 89  - The latitude range to be considered    
DIST_Z0 8 DIST_Z1 16 - The altitude range to be considered   
DIST_ZSCORE - can be used for filtering outliers
