## Predefined pressure level sets

With the parameter `MET_PRESS_LEVEL_DEF` vertical pressure levels can
be selected from a number of definitions.  The definitions follow the
[ECMWF's model level definitions](https://confluence.ecmwf.int/display/UDOC/Model+level+definitions).
Additional levels up to ~1045 hPa have been added to avoid
extrapolation errors near the surface.  Per default, the value is set
to -1 and hence pre-defined pressure level sets ignored.  The
recommended pressure level set for ERA5 is number 6.

| MET_PRESS_LEVEL_DEF |  Name  | bottom      | top      | number of levels |
| ------------------- | ------ | ----------- | -------- | ---------------- |
| 0                   |  L137  | 1044.45 hPa | 0.02 hPa | 138              |
| 1                   |   L91  | 1044.45 hPa | 0.02 hPa |  92              |
| 2                   |   L60  | 1044.45 hPa | 0.01 hPa |  60              |
| 3                   |  L137  | 1044.45 hPa | 0.02 hPa | 147              |
| 4                   |   L91  | 1044.45 hPa | 0.02 hPa | 101              |
| 5                   |   L60  | 1044.45 hPa | 0.01 hPa |  62              |
| 6                   |  L137  | 1044.45 hPa | 0.01 hPa | 137              |
| 7                   |   L60  | 1046.13 hPa | 0.1 hPa  |  59              |
