# tropo_sample

Sample the tropopause height from tropopause data file <tropo.nc> created by e.g. tropo at a given time, latitude, and longitude provided in the atmosphere data file <atm_in>. The tropopause data files may contain snapshots at a certain time or climatologies. The results in terms of tropopause geometric height, pressure, temperature, and water vapour volume mixing ratio are written to the tropopause sample output file <sample.tab>.

```
# calling tropo_sample
$ tropo_sample  <ctl> <sample.tab> <tropo.nc> <var> <atm_in>
```

The type of tropopause information to be extracted is selected with the
\<var\> argument. This argument is the variable prefix used in the
NetCDF file created by [tropo](tropo.md). The following prefixes are
available:

| parameter | description |
| :--------- | :----------- |
| `clp` | cold point tropopause |
| `wmo_1st` | WMO first lapse rate tropopause |
| `wmo_2nd` | WMO second lapse rate tropopause |
| `dyn` | dynamical tropopause |

For example, using `dyn` samples the variables `dyn_z`, `dyn_p`,
`dyn_t`, and, if available in the file, `dyn_q` and `dyn_o3`.

The following configuration parameters are effective in tropo_sample:

| parameter |  purpose | default |
| :--------- | :-------- | --------: |
| TROPO_SAMPLE_METHOD | Define interpolation method. Default=1: linear interpolation. Everything else: nearest neighbor. | 1 |
