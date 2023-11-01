# met_subgrid

Calculates the grid-scale standard deviations of the horizontal wind and vertical velocity fields as used in the subgrid-scale wind parameterization. The subgrid-scale wind parametrization is described in Stohl et al. (2005).

* Stohl, A., Forster, C., Frank, A., Seibert, P., and Wotawa, G.: Technical note: The Lagrangian particle dispersion model FLEXPART version 6.2, Atmos. Chem. Phys., 5, 2461–2474, https://doi.org/10.5194/acp-5-2461-2005, 2005. 

```
# calling met_subgrid
$ ./met_subgrid  <ctl> <subgrid.tab> <met0> <met1> [ <met0b> <met1b> ... ]
```

The required control parameters are:
* ctl: the control parameter file
* subgrid.tab: the output file with standard deviations of the wind and velocity (in ASCII table format)
* met*: pairs of meteorological data files from which the standard deviations are calculated
