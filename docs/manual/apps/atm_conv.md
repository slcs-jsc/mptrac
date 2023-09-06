# atm_conv

``` atm_conv ``` converts atmosphere files from one format into the other. Possible input files are ASCII, binary, netcdf and netcdf files following format guidelines for CLaMS. The data can be converted to these files as well. The usage is as follows:

```
# calling atm_conv
$ ./atm_conv  <ctl> <atm_in> <atm_type_in> <atm_out> <atm_type_out>
```

* ```ctl``` : The control file is used to set the atmospheric types and files if not set in the command line. It is also needed if it is empty. 
* ```atm_in```: Input atmospheric file
* ```atm_type_in```: Format of input atmospheric file
* ```atm_out```: Output atmospheric file
* ```atm_type_in```: Format of output atmospheric file

The type of data can be chosen from the options below:

| ATM_TYPE | output format |
|---|---|
| 0  | ASCII (default) |
| 1  | binary  |
| 2  | netcdf  |
| 3  | netcdf (CLaMS: trajectory and position file) |
| 4  | netcdf (CLaMS: position file)  |

Note:  If trajectory files similar to those of CLaMS are required (```atm_type=3```), the stop time, as seen below, must be set in the control file. However, to obtain position files as for CLaMS the option ```atm_type=4```, which does not rely on these setting, is recommended instead.

```                                           
cat > atm_conv.ctl <<EOF
T_STOP = 360547200.00          
EOF   
```

Example: 
> ```atm_conv atm_conv.ctl input.tab 0 output.nc 2```
> 
> , converts the input.tab ASCII file into the new > netcdf file output.nc.

