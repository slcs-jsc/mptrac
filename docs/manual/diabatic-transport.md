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
MET_PRESS_LEVEL_DEF = 0
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
ignored. For ERA5 data, it is recommended to use set 6.

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
