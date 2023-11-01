# tropo

This application creates a tropopause data set from meteorological data.

```
# calling tropo
$ ./tropo  <ctl> <tropo.nc> <met0> [<met1> ...]
```
The required arguments are:
* ctl: The file containing the control parameters.
* tropo.nc: The output file in which the climatology should be stored.
* met0 \[\<met1\> ...\]: The meteorological data from which the data set should be created.

The following paramerters can be defined in <.ctl> file:
| parameter  | description                                    | default | Options                                            |
|:---------- |:---------------------------------------------- | ------- | -------------------------------------------------- |
| DT_MET     | time step of meteorological data [s].          | 21600   |                                                    |
| MET_TROPO  | tropopause definition                          | 3       | 0:none, 1:clim, 2:cold point, 3:WMO_1ST, 4:WMO_2nd, 5:dynamical |
| MET_TROPO_PV | dynamical tropopause threshold (PV)          | 3.5     |
| TROPO_LON0 | minimum longitude                              | -180    |                                                    |
| TROPO_LON1 | maximum longitude                              | 180     |                                                    |
| TROPO_DLON | distance between longitudes                    | -999    |                                                    |
| TROPO_LAT0 | minimum latitude                               | -90     |                                                    |
| TROPO_LAT1 | maximum latitude                               | 90      |                                                    |
| TROPO_DLAT | distance between latitudes                     | -999    |                                                    |
| TROPO_H2O  | extract water vapor mixing ratio at tropopause | 1       |                                                    |
| TROPO_O3   | extract ozone volume at tropopause             | 1       |                                                    |
|            |                                                |         |                                                    |

Example:
```
# define the control file
cat > tropo.ctl <<EOF
DT_MET = 21600  
MET_TROPO = 2
TROPO_DLON = 6
TROPO_DLAT = 3
EOF
# running the tropo application
../src/tropo tropo.ctl tropo_2011_06_05.nc ../ei_2011_06_05_00.nc ../ei_2011_06_06_00.nc
```
here:
* ../src/tropo: calling tropo
* tropo.ctl: using control file
* tropo_2011_06_05.nc: output of the tropopause 
* ../ei_...: inputs of the meteorological data (data directory should be indicated)
* At JSC, Era_interim data can be found: /p/fastdata/slmet/slmet111/met_data/ecmwf/era_interim/pressure_0.75deg_v2/nc/2011/ei_...nc
