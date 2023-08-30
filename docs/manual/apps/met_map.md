# met_map 

Extracts map from meteorological data.

```
# calling met_map
$ ./met_map  <ctl> <map.tab> <met0> [<met1> ...]
```

The required parameters are as follows:

* ctl: The file containing the control parameters. The control parameters can also be appended at the end of the function call, giving the flag name first, followed by the parameter, all separated by a space.
* met_map: The output file in which the map should be stored.
* met0 \[\<met1\> ...]: The meteorological data from which the climatology should be created.

Calling met_map and creating out of the 0.3° x 0.3° ER5 data a data set on a 1° x 1° would look like as follows:

```
# Calling met_map
$./met_map - map_era5_2017_01_08_17_1_1.tab era5_2017_01_08_17.nc MAP_DLON 1 MAP_DLAT 1
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

To extract a map from a compressed data file (as e.g. a ZFP file) one would need to give also the meteorological data file format as control parameter. Calling met_map would then look as follows:

```
# Calling met_map for a zfp file
$./met_map - map_era5_2017_01_08_17.tab era5_2017_01_08_17.zfp MET_TYPE 3
```

Currently, to enable the usage of zfp and zstd format MPTRAC needs to be compiled using ZSTD=1 and ZFP=1:

```
# Run make enabling the usage of the ZFP and ZSTD compression
$ make ZSTD=1 ZFP=1
````
