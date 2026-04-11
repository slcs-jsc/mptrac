# tnat

This application calculates the PSC formation temperatures (in Kelvin).

```
# Calling tnat
$./tnat <p> <h2o> <hno3>
```

Required parameters are:
* p: The atmospheric pressure in hPa.
* h2o: The volume mixing ratio (ppv) of water vapour.
* hno3: The volume mixing ratio (ppv) of nitric acid.

Assuming a pressure of 50 hPa, a H<sub>2</sub>O volume mixing ratio of
5e-6 ppv, and a HNO<sub>3</sub> volume mixing ratio of 15e-9 ppv:
```
$./tnat 50 5e-6 15e-9
```
one receives the following output:

```
$ p = 50 hPa
$ q_H2O = 5e-06 ppv
$ q_HNO3 = 1.5e-08 ppv
$ T_dew = 184.543 K
$ T_ice = 188.559 K
$ T_NAT = 196.312 K
```
