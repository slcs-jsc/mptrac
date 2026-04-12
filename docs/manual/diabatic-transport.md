## Set-up for diabatic transport calculations

For diabatic transport, MPTRAC loads the fields $\dot{\zeta}$,
$\zeta$, and $p$ in the original $\eta$ hybrid coordinate system from
the CLaMS dataset. This approach avoids additional interpolation from
$\zeta$ to $p$. However, other fields remain in pressure coordinates
since several MPTRAC modules require a pressure-based formulation.

To ensure compatibility, $\zeta$ and $p$ are converted back and forth
before and after advection. Consequently, vertical model pressure
levels must be explicitly defined. The following configuration is
required for proper reading and transformation of CLaMS data:

```
MET_CLAMS = 1
MET_VERT_COORD = 1
MET_PRESS_LEVEL_DEF = 3
```

To enable vertical advection in $\zeta$ coordinates with diabatic
transport, set the following parameters and include $\zeta$ as an
atmospheric quantity:

```
ADVECT_VERT_COORD = 1
NQ = 1
QNT_NAME[0] = zeta
```

## Predefined pressure level sets

The `MET_PRESS_LEVEL_DEF` parameter allows selection from several
predefined vertical pressure level sets, based on [ECMWF model level
definitions](https://confluence.ecmwf.int/display/UDOC/Model+level+definitions). Additional
near-surface levels have been added (down to ~1045 hPa) to reduce
extrapolation errors.

By default, `MET_PRESS_LEVEL_DEF = -1`, meaning predefined sets are
ignored. Sets 0, 1, and 2 are disabled because their near-surface spacing can
require extrapolation where the local surface pressure exceeds 1013.25 hPa.
Set 6 is the recommended default for ERA5-like model-level data. Set 7 keeps
the coarser ERA-Interim/L60 spacing near the surface and is intended for
ERA-Interim- and JRA-55-like data.

| MET_PRESS_LEVEL_DEF |  Name  | bottom      | top      | number of levels | notes |
| ------------------- | ------ | ----------- | -------- | ---------------- | ----- |
| 0                   |  L137  | disabled    | disabled | disabled         | replaced by set 3 |
| 1                   |   L91  | disabled    | disabled | disabled         | replaced by set 4 |
| 2                   |   L60  | disabled    | disabled | disabled         | replaced by set 5 |
| 3                   |  L137  | 1044.45 hPa | 0.02 hPa | 147              | extended L137 with fine near-surface spacing |
| 4                   |   L91  | 1044.45 hPa | 0.02 hPa | 101              | extended L91 with fine near-surface spacing |
| 5                   |   L60  | 1044.45 hPa | 0.01 hPa |  62              | extended L60 |
| 6                   |  L137  | 1044.45 hPa | 0.01 hPa | 137              | recommended default for ERA5-like data |
| 7                   |   L60  | 1046.13 hPa | 0.1 hPa  |  59              | ERA-Interim/JRA-55-like L60 grid with extrapolated lower levels |

Recommended settings for common meteorological data sources:

| Data source | MET_PRESS_LEVEL_DEF | notes |
| ----------- | ------------------- | ----- |
| ERA5        | 6                   | recommended default |
| ERA-Interim | 7                   | keeps the ERA-Interim/L60-like spacing and 0.1 hPa top |
| MERRA-2     | 6                   | uses a 0.01 hPa top, matching the native upper model domain |
| JRA-3Q      | 6                   | uses a 0.01 hPa top for the extended upper model domain |
| JRA-55      | 7                   | closer to the 60-level, 0.1 hPa-top setup |
