# met_map 

Extracts map from meteorological data.

```
# calling met_map
$ ./met_map  <ctl> <map.tab> <met0> [<met1> ...]
```

The required parameters are as follows:

* ctl: The file containing the control parameters.
* met_map: The output file in which the map should be stored.
* met0 \[\<met1\> ...]: The meteorological data from which the climatology should be created.

Set optional `MAP_*` and `MET_*` parameters in the control file for
reliable operation. Appending them after the meteorological input files
is not recommended for this tool.

Example control-file settings for extracting a 1° x 1° map from 0.3° x
0.3° ERA5 data:

```
MAP_DLON = 1
MAP_DLAT = 1
```

Example call:

```
./met_map map.ctl map_era5_2017_01_08_17_1_1.tab era5_2017_01_08_17.nc
```

The optional control parameters required are listed as follows:

| parameter | purpose | default | 
|:-----------|:---------|---------:|
| MAP_Z0 | The altitude to extract the map data [km]| 10 |
| MAP_LON0 | Longitude range start | -180 |
| MAP_LON1 | Longitude range end | 180 |
| MAP_DLON | Longitude interval | -999 |
| MAP_LAT0 | Latitude range start | -90 |
| MAP_LAT1 | Latitude range end | 90 |
| MAP_DLAT | Latitude interval | -999 |
| MAP_THETA | Use theta level instead of altitude | -999 |

To extract a map from a compressed meteorological data file such as ZFP,
set the meteorological file type in the control file:

```
MET_TYPE = 3
```

Example call:

```
./met_map map.ctl map_era5_2017_01_08_17.tab era5_2017_01_08_17.zfp
```

Currently, to enable the usage of zfp and zstd format MPTRAC needs to be compiled using ZSTD=1 and ZFP=1:

```
# Run make enabling the usage of the ZFP and ZSTD compression
$ make ZSTD=1 ZFP=1
```
