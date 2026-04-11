# atm_stat

Calculate air parcel statistics.

```
# calling atm_stat
$ atm_stat  <ctl> <stat.tab> <param> <atm1> [<atm2> ...]
```

Statistical parameters that can be calculated are:

* mean
* stddev
* min
* max
* skew
* kurt
* absdev
* median
* mad

The statistics are calculated with the GNU Scientific Library. More information on the [GNU Scientific Library statistics functions](https://www.gnu.org/software/gsl/doc/html/statistics.html) can be found on their webpage.
