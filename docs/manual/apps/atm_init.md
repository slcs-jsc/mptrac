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
| INIT_RELZ0 | Release start position relative to the local PBL depth `z/zi` | -999 |
| INIT_RELZ1 | Release end position relative to the local PBL depth `z/zi` | -999 |
| INIT_DRELZ | Release interval in relative PBL coordinates | 1 |
| INIT_RANDOM_RELZ | Uniform random sampling between `INIT_RELZ0` and `INIT_RELZ1` | 0 |
| INIT_LON0 | Release longitude range start | 0 |
| INIT_LON1 | Release longitude range end | 0 |
| INIT_DLON | Release longitude interval | 1 |
| INIT_LAT0 | Release latitude range start | 0 |
| INIT_LAT1 | Release latitude range end | 0 |
| INIT_DLAT | Release latitude interval | 1 |

For the standard longitude/latitude grid these horizontal coordinates
are given in degrees. For the special-purpose UTM-grid mode
(`MET_COORD_TYPE = 1`), the same parameters refer to UTM easting and
northing in metres. The horizontal spread parameters (`INIT_SLON`,
`INIT_SLAT`, `INIT_ULON`, `INIT_ULAT`, and `INIT_SX`) use the same
horizontal coordinate units as the respective met grid: degrees for longitude/latitude grids and
metres for UTM grids.

The following optional parameters can be used for controlling the release of air parcels:

* INIT_REP: This parameter controls the number of air parcels in each release position. The default value is 1.
* INIT_MASS: Total release mass [kg]. The default value is 0.
* INIT_VMR: Volume mixing ratio of each air parcel [ppv]. The default value is 0.
* INIT_ST/INIT_SZ/INIT_SLON/INIT_SLAT: If this parameter is set to non-zero value, the air parcels are released with a Gaussian distribution in time/altitude/longitude/latitude. The parameter value represents the full width at half maximum. 
* INIT_UT/INIT_UZ/INIT_ULON/INIT_ULAT: If this parameter is set to non-zero value, the air parcels are released with a uniform distribution in time/altitude/longitude/latitude. The parameter value represents release range.
* INIT_EVENLY: If this parameter is set to 1, the number of air parcels released will be weighted by the cosine of latitude so that the air parcels are evenly distributed globally.
* INIT_COL_MASS: If this parameter is set to 1, air parcels are sampled uniformly in surface pressure, proportional to the column mass above the model top. The surface pressure is taken from the meteo data. Requires meteo data and a fixed release time (INIT_T0 == INIT_T1, INIT_ST = INIT_UT = 0). Default is 0.
* INIT_WELL_MIXED: If this parameter is set to 1, air parcels are distributed uniformly in pressure between the model top and the surface pressure at their location, resulting in a well-mixed tracer distribution. Requires INIT_EVENLY=1, INIT_COL_MASS=1, and a fixed release time. Default is 0.
* INIT_NP: Target number of air parcels when using INIT_WELL_MIXED sampling. The simulation will retry the grid sampling until this number is reached. Must be between 0 and NP. Requires INIT_WELL_MIXED=1. Default is 0 (grid-based sampling, no retry loop).
* INIT_IDX_OFFSET: Offset value added to air parcel indices. This parameter allows customization of the starting index for air parcels, which can be useful for domain decomposition or when combining multiple initialization runs. The default value is 0. 

Relative-PBL releases use the local surface pressure and local PBL pressure from the meteorological data to convert `z/zi` to pressure. Therefore `INIT_RELZ0` and `INIT_RELZ1` require meteorological input and a fixed release time (`INIT_T0 == INIT_T1`, `INIT_ST = INIT_UT = 0`). They cannot be combined with `INIT_Z0` and `INIT_Z1`.

If `INIT_RANDOM_RELZ = 0`, the code constructs a regular release grid in relative PBL coordinates using `INIT_RELZ0`, `INIT_RELZ1`, and `INIT_DRELZ`. If `INIT_RANDOM_RELZ = 1`, the vertical position is sampled uniformly at random between `INIT_RELZ0` and `INIT_RELZ1` for each parcel.

`INIT_REP` and `INIT_NP` serve different purposes: `INIT_REP` repeats every release position in the configured grid, while `INIT_NP` is only used together with `INIT_WELL_MIXED = 1` to keep resampling until a target total number of parcels has been accepted.

**Note:** Instead of using the control file the configuration parameters can also be appended to the function call. 
