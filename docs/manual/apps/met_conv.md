# met_conv

Conversion of meteorological data

```
# calling met_conv
./met_conv <ctl> <met_in> <met_in_type> <met_out> <met_out_type>
```
```
# An example command line for running atm_conv looks like
./met_conv trac.ctl era5_2017_01_08_17.nc 0 era5_2017_01_08_17.zstd 4
```
Here the meteorological data file in netcdf format is converted (compressed) to the zstd format.

The following input/output files are supported by met_conv:

|Type In/Out| Format |
|:--------- |:-----:|
|0          | netcdf |
|1          | binary |
|2          | pck    |
|3          | zfp    |
|4          | zstd |

whereby netcdf and binary are standard data formats. Pck (layer packing, https://gmd.copernicus.org/articles/10/413/2017), zfp (https://computing.llnl.gov/projects/zfp) and zst (zstandard, https://github.com/facebook/zstd) are compression formats. 

Currently to enable the usage of zfp and zstd format MPTRAC needs to be compiled using ZSTD=1 and ZFP=1:

```
# Run make enabling the usage of the ZFP and ZSTD compression
$ make ZSTD=1 ZFP=1
````

