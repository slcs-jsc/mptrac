# tropo_zm

This application extracts the zonal mean of a tropopause data set.

```
# Calling tropo_zm
./tropo_zm <ctl> <zm.tab> <var> <tropo.nc>
```
The required parameters are:
* ctrl: The file containing the control parameters.
* zm.tab: The ouput file in which the zonal means should be stored.
* var: The variables???? Which can be chosen?
* tropo.nc: The tropopause data set of which the zonal means should be calculated.

Add an example how to use tropo_zm. Which input/output formats are supported? How does the ouput of tropo_zm look like?
Is "var" really used? From the code it looks like that generally all (or at least several) parameters are written to the output file.
