# tropo_zm

This application extracts the zonal mean of a tropopause data set.

```
# Calling tropo_zm
./tropo_zm <ctl> <zm.tab> <var> <tropo.nc>
```
The required parameters are:
* ctrl: The file containing the control parameters.
* zm.tab: The ouput file in which the zonal means should be stored.
* var: The tropopause variables
* tropo.nc: The tropopause data set of which the zonal means should be calculated.

| var     | description                     |
|:------- |:------------------------------- |
| clp     | cold point tropopause           |
| dyn     | dynamical tropopause            |
| wmo_1st | WMO first lapse rate tropopause |
| wmo_2nd | WMO secon lapse rate tropopause |

Example:
```
for var in clp dyn wmo_1st wmo_2nd ; do
    ../src/tropo_zm - zm_$var.tab $var tropo_2011_06_05.nc
done
```
Here:
* ../src/tropo_zm: calling tropo_zm
* no control file is given
* zm_$var.tab: output of the zonal mean tropopause with the same latitude as the input data
* tropo_2011_06_05.nc: input, the tropopause dataset produced by "tropo"

output:
* $1 = time [s]
* $2 = latitude [deg]
* $3 = tropopause height (mean) [km]
* $4 = tropopause pressure (mean) [hPa]
* $5 = tropopause temperature (mean) [K]
* $6 = tropopause water vapor (mean) [ppv]
* $7 = tropopause ozone (mean) [ppv]
* $8 = tropopause height (sigma) [km]
* $9 = tropopause pressure (sigma) [hPa]
* $10 = tropopause temperature (sigma) [K]
* $11 = tropopause water vapor (sigma) [ppv]
* $12 = tropopause ozone (sigma) [ppv]
* $13 = number of data points
* $14 = occurrence frequency [%]
