# atm_init

The application atm_init creates an atmospheric data file with initial air parcel positions. The code will generate air parcels in a simple manner with user-defined settings.

```
# calling atm_init
$ ./atm_init  <ctl> <atm_out>
```
The following specific configuration parameters can be used to determine the release position:

| parameter | purpose | default | 
|:-----------|:---------|---------:|
| INIT_T0 | Release time range start [s] | 0 |
| INIT_T1 | Release time range end [s] | 0 |
| INIT_DT | Release time interval [s] | 1 |
| INIT_Z0 | Release altitude range start [km] | 0 |
| INIT_Z1 | Release altitude range end [km] | 0 |
| INIT_DZ | Release altitude interval [km] | 1 |
| INIT_LON0 | Release longitude range start | 0 |
| INIT_LON1 | Release longitude range end | 0 |
| INIT_DLON | Release longitude interval | 1 |
| INIT_LAT0 | Release latitude range start | 0 |
| INIT_LAT1 | Release latitude range end | 0 |
| INIT_DLAT | Release latitude interval | 1 |

Followings are some optional parameters for controlling the release of air parcels:

* INIT_REP: This parameter control the number of air parcels in each release position. The default value is 1.
* INIT_MASS: Total release mass [kg]. The default value is 0.
* INIT_VMR: Volume mixing ratio of each air parcel [ppv]. The default value is 0.
* INIT_ST/INIT_SZ/INIT_SLON/INIT_SLAT: If this parameter is set to non-zero value, the air parcels are released with a Gaussian distribution in time/altitude/longitude/latitude. The parameter value represents the full width at half maximum. 
* INIT_UT/INIT_UZ/INIT_ULON/INIT_ULAT: If this parameter is set to non-zero value, the air parcels are released with a uniform distribution in time/altitude/longitude/latitude. The parameter value represents release range.
* INIT_EVENLY: If this parameter is set to 1, the number of air parcels released will be weighted by the cosine of latitude so that the air parcels are evenly distributed globally. 

**Note:** Instead of using the control file the configuration parameters can also be appended to the function call. 