# atm_select

Extract subsets of air parcels from atmospheric data files.

The calling sequence is:

```
# calling atm_select
$ atm_select  <ctl> <atm_select> <atm1> [<atm2> ...]
```
The atmospheric data file(s) <atm1> [<atm2> …] contain(s) a list of data points with the structure:

```
time altitude longitude latitude quantity_1 quantity_2 ...
```
Following the selection criteria provided in the control file <ctl> atm_select extracts a subset of the data points and writes them into a new atmospheric data file with the file name provided in <atm_select>. The following atm_select specific configuration parameters can be used:

| parameter | purpose  | default |
| :--- | :----------- | -----------:
| SELECT_STRIDE | select every i-th data point |  1 |
| SELECT_IP0    | set index range start  | 0 | 
| SELECT_IP1    | set index range end  | 0 |
| SELECT_T0 | set time range start [s] | 0 |
| SELECT_T1 | set time range end [s] | 0 |
| SELECT_Z0 | set altitude range start [km] | 0 |
| SELECT_Z1 | set altitude range end [km] | 0 |
| SELECT_LON0 | set longitude range start [deg] | 0 |
| SELECT_LON1 | set longitude range end [deg] | 0 |
| SELECT_LAT0 | set latitude range start [deg] | 0 |
| SELECT_LAT1 | set latitude range end [deg] | 0 |
| SELECT_R0 | set radius range start [km] | 0 |
| SELECT_R1 | set radius range end [km] | 0 |
| SELECT_RLON | set station longitude [deg] | 0 |
| SELECT_RLAT | set station latitude [deg] | 0 |

The parameters should be given as an indexed list in the control file (example).

**Note:** The parameters can also be appended to the function call, but they will cause error messages, as a loop runs over all input parameters >= 3 and tries to read files named as the control parameters. Despite the annoying error messages the app will create a correct output file.
