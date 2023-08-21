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


