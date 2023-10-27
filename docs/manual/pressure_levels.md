## Predefined pressure levels


#JS - Comment - First User Perspective
#Open Questions:
#For what purpose/application is it usefull/necessary to define pressure levels?
#Are the pressure levels read from met file per default, anyway?
#Where does the pressure level information come from?
#The bottom pressure level in the ECMWFs model level definitions link is 1013.25 and not 1044.45?




With the parameter ``` PRESS_LEVEL_DEF ``` vertical pressure levels can be selected from a number of definitions. 
The definitions follow the [ECMWFs model level definitions](https://confluence.ecmwf.int/display/UDOC/Model+level+definitions). 
Per default, the value is set to -1 and hence ignored.

| PRESS_LEVEL_DEF |  Name  | bottom  | top | nbr. of levels    |
| --------------- | --- | --- | --- | --- |
| 0               |   L137  |  1044.45 hPa|0.02 hPa| 137    |
| 1               |   L91  |1044.45 hPa|0.02 hPa|  92   |
| 2               |   L60  |1044.45 hPa|0.01 hPa|   61  |
