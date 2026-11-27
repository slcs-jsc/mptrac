# atm_split

Split air parcels into a larger number of parcels.

The calling sequence is:

```
# calling atm_split
$ ./atm_split  <ctl> <atm_in> <atm_out>
```
The atmospheric input data file <atm_in> contains a list of data points with the structure:

```
time altitude longitude latitude quantity_1 quantity_2 ...
```

Based on the locations provided in the input atmospheric data file <atm_in> atm_split generates a specified total number of data points by generating new data points in a region defined by the configuration parameters and stores the result in a new atmospheric data file with the file name provided in <atm_out>.

The following atm_split specific configuration parameters can be used:

| parameter | purpose | default | 
|:-----------|:---------|---------:|
| SPLIT_N | total number of data points to generate |
| SPLIT_M | total mass assigned to the new particles [kg]; negative values use the input total | -999 |
| SPLIT_UM | absolute width of a uniform child-mass perturbation [kg] | 0 |
| SPLIT_DT | delta time [s] | 0 |
| SPLIT_T0 | time window start [s] | 0 |
| SPLIT_T1 | time window end [s] | 0 |
| SPLIT_DZ | vertical radius [km] | 0 |
| SPLIT_Z0 | altitude range start [km] | 0 |
| SPLIT_Z1 | altitude range end [km] | 0 |
| SPLIT_DX | horizontal radius [km] | 0 |
| SPLIT_LON0 | longitude range start | 0 |
| SPLIT_LON1 | longitude range end; LON0<LON1 |	0 |
| SPLIT_LAT0 | latitude range start | 0 |
| SPLIT_LAT1 | latitude range end; LAT0<LAT1 | 0 |
| SPLIT_KERNEL | kernel file name | – |

If a range is provided the data points are distributed uniformly. If a distance is provided (DT, DZ, DX), the data points are distributed with a Gaussian random variate where e.g. ‘distance’/2.3548 = FWHM.

With `N = SPLIT_N`, each child receives `m_child = [M + U * (V - 0.5)] / N`, where `M` is the input total mass or `SPLIT_M`, `U = SPLIT_UM`, and each `V` is drawn independently and uniformly from [0, 1]. `SPLIT_UM` is an absolute uniform perturbation, not a relative or Gaussian uncertainty. At `SPLIT_UM = 0`, the selected total mass is preserved up to rounding; for positive `SPLIT_UM`, it is preserved only in expectation, and sufficiently large perturbations can produce negative child masses. Other configured quantities are copied from each selected parent except where explicitly reassigned; sums of other extensive quantities, such as accumulated mass-loss diagnostics, need not be preserved.

**Note:** Instead of using the control file the configuration parameters can also be appended to the function call. 
