# met_sample

Sample meteorological data at given geolocations 

```
# calling met_sample
$ ./met_sample  <ctl> <sample.tab> <atm_in>
```
The required input parameters are:
* ctl: The control file, defining the atmospheric grid.
* sample.tab: Output ascii file with re-sampled meteorological data to the locations provided in atm_in.
* atm_in: Input file of atm-type. Expected input parameters are a list of data sets containing time, altitude, longitude, and latitude.

Add an example for  calling/applying met_sample.
```
# Calling met_sample
$./met_sample sample.ctl new_file.tab volcanoes.tab 
```
The following configuration parameters are effective in met_sample:

| parameter | purpose | default | 
|:-----------|:---------|---------:|
| SAMPLE_GEOPOT     | sample geopotential height | 0 |
| SAMPLE_GRID_TIME  | sample time      | 0 |
| SAMPLE_GRID_Z     | sample altitude  | 0 |
| SAMPLE_GRID_LON   | sample longitude | 0 |
| SAMPLE_GRID_LAT   | sample latitude  | 0 |
