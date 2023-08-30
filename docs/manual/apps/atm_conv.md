# atm_conv

Converts file format of air parcel data files.

```
# calling atm_conv
$ ./atm_conv  <ctl> <atm_in> <atm_in_type> <atm_out> <atm_out_type>
```
The required paramters are as follows:

* ctl: The control file.
* atm_in: The name of the air parcel data file that should be converted.
* atm_in_type: The format type of air parcel data input file.
* atm_out: The name of the air parcel in the new data format. 
* atm_out_type: The type into the air parcel data file should be converted.

Types of atmospheric data:

|Type In/Out| Format |
|:--------- |:-----:|
| 0         | ASCII |    
| 1         | binary| 
| 2         | netcdf|
| 3         | clams |  

Correct?? Needs to be checked.
