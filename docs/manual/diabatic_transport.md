## Set-up for diabatic transport calculations

For diabatic transport calculations, MPTRAC loads the fields $\dot{\zeta}$, $\zeta$ and $p$ in the original $\eta$ hybrid coordinates of the CLaMS data. This is done to avoid an additional interpolation step from $\zeta$ to $p$. However, other fields are kept in pressure coordinates because other modules rely on a formulation in $p$. Before and after advection, pressure and zeta are converted into each other. Therefore, the vertical model pressure levels must also be defined. Reading the CLaMS data and transforming it to the correct levels requires the following settings:

```
MET_CLAMS = 1
MET_VERT_COORD = 1
MET_PRESS_LEVEL_DEF = 0
```

Moreover, to activate the advection in $\zeta$ coordinates and with diabatic transport, the parameter ```ADVECT_VERT_COORD``` must be set and the quantity $\zeta$ must be added to the atmosphere quantities.

```
ADVECT_VERT_COORD = 1
NQ = 1
QNT_NAME[0] = zeta
```
