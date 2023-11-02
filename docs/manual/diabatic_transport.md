## Set-up for diabatic transport calculations

For diabatic transport calculations, MPTRAC loads the fields $\dot{\zeta}$, $\zeta$ and $p$ in the original $\eta$ hybrid coordinates of the CLaMS data. This is done to avoid an additional interpolation step from $\zeta$ to print. However, other fields are kept in pressure coordinates because other modules rely on a formulation in $p$. Before and after advection, pressure and zeta are converted. Therefore, the vertical model pressure levels must also be defined. Reading the CLaMS data and transforming it to the correct levels requires the following settings:

```
CLAMS_MET_DATA = 1
VERT_COORD_MET = 1
PRESS_LEVEL_DEF = 0
```

Moreover, to activate the advection in $\zeta$ coordinates and with diabatic transport, the below parameter ```VERT_COORD_AP``` must be set and the quantity $\zeta$ must be added to the atmosphere quantities.

```
VERT_COORD_AP = 1
NQ = 1
QNT_NAME[0] = zeta
```

Finally, if the advection is coupled with parameterisations such as the convection parameterisation, the following control parameter must be added: 

```
CPL_ZETA_PRESS_MODULES = 1

```


