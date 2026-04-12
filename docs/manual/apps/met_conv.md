# met_conv

Conversion of meteorological data

```
# calling met_conv
./met_conv <ctl> <met_in> <met_in_type> <met_out> <met_out_type>
```
The required parameters are as follows:

* ctl: The control file.
* met_in: The name of the meteorological file that should be converted.
* met_in_type: The format type of the meteorological input file.
* met_out: The name of the meteorological file in the new data format. 
* met_out_type: The type into the meteorological data file should be converted.

```
# An example command line for running met_conv looks like
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
|5          | cms    |
|7          | sz3    |

whereby netCDF and binary are standard data formats. Pck (layer packing, https://gmd.copernicus.org/articles/10/413/2017), zfp (https://computing.llnl.gov/projects/zfp), zstd (zstandard, https://github.com/facebook/zstd), cms, and sz3 are compression formats.

The output formats store the pressure-level meteorological fields used by
standard trajectory simulations. The generic MPTRAC netCDF writer also
includes derived 2-D diagnostics, but it does not write native model-level
fields such as `PL`, `UL`, `VL`, `WL`, `ZETA`, and `ZETA_DOT`. Keep the
original netCDF meteo files for diabatic or zeta-coordinate workflows.

The optional compression formats need to be enabled at compile time. For example:

```
# Enable selected compression backends.
$ make ZFP=1 ZSTD=1 SZ3=1 CMS=1
```
