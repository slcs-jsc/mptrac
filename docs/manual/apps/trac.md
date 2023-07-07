# trac

The major app for calculating forward or backward trajectories is trac. It is called with at least four (4) arguments:

```
# calling trac
$ ./trac  <dirlist> <ctl> <atm_in> <metbase>
```

The required arguments are:

* dirlist: A file containing directories to be processed. Each directory has to have an own control parameter file and a starting point file.

* ctl: In the control parameter file the configuration parameters can be set.

* atm_in: The starting point file contains a list of starting points for the trajectory calculation.

* metbase: Here, the path to the meteorological data files and their basename shall be given.

Example from the repository:

```
# minimal required input to run trac
$ ./trac data/dirlist trac.ctl atm_split.tab meteo/ei
```

In addition it is possible to append control parameters, giving the flag name first, followed by the parameter, all separated by a space. Example:

```
#
$ ./trac data/dirlist trac.ctl atm_split.tab meteo/ei ATM_BASENAME atm_diff GRID_BASENAME grid_diff
```

The following control parameters are effective in trac:

Todo

add flaglist here …. 
