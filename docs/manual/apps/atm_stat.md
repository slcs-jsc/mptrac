# atm_stat

Calculate air parcel statistics.

```
# calling atm_stat
$ atm_stat  <ctl> <stat.tab> <param> <atm1> [<atm2> ...]
````

Statistical paramters that can be calculated are:

* mean
* stddv
* min
* max
* skew
* kurt
* absdev
* median
* mad

The here used statistics are the gnu statistic library. More information on the [gnu statistics library](https://www.gnu.org/software/gsl/doc/html/statistics.html) can be found 
on their webpage. 


#JS - Comment - First User Perspective
#Open Questions:
#It is not quite clear to me, what is the difference between atm_dist and atm_stat?
#control parameter: absdev or abs_dev?


