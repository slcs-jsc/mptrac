# Model physics

## Advection

## Diffusion

## Dry deposition

## Wet deposition

## OH chemistry

MPTRAC provides a module to simulate the chemical reaction of a species with the hydroxyl radical (OH).

### Introduction

Various species are decomposed by chemical reaction with OH, which causes a loss of the mass of the species over time. There are two types of reaction implemented in OH chemistry module. A control parameter oh_chem_reaction should be defined to choose a reaction type. Setting oh_chem_reaction to 2 means bimolecular reaction while oh_chem_reaction to 3 means termolecular reaction. In category of bimolecular reaction, the rate constant is expressed by Arrhenius formula:

```
    k(T) = A × exp (− E / RT)
```

In category of termolecular reaction, it requires a inert component M (e.g. N2 and O2) to stabilize the excited intermediate.


```
    X + OH + M -> XOH + M
```

Therefore, the reaction rate shows a temperature and pressure dependence and changes smoothly between the low-pressure-limiting rate k0 and high-pressure-limiting rate ki. The reaction rate k of this reaction is calculated from:

```
    k = k0 * M / (1. + k0 * M / ki) * pow(0.6, 1. / (1. + c * c))

    c = log10(k0 * M / ki)
```

The molecular density of air is calculated from:

```
    M = 7.243e21 * (p / p0) / T
```

The low and high pressure limits of the reaction rate are given by:

```
    k0 = k0n * pow(T / Trefn, -n)
    ki = kim * pow(T / Trefm, -m)
```

To activate the MPTRAC OH chemistry module, the constants k0n, kin, n, and m need to be specified as control parameters:

```
    OH_CHEM[0] = k0n
    OH_CHEM[1] = n
    OH_CHEM[2] = kim
    OH_CHEM[3] = m
```

The parameters n or m can be set to zero in order to neglect the temperature dependence of the low- or high-pressure limits.

Note that the simulation of loss of mass by exponential decay over time may have to be switched off if the OH chemistry module is activated (set TDEC_TROP and TDEC_STRAT to zero).

### Atmospheric concentrations of OH

The atmospheric concentrations of OH are obtained by linear interpolation of the monthly mean zonal mean climatology of Grooss and Russell [2005], which has been derived from HALOE satellite measurements.

### Reaction rates for sulfur dioxide

The parameters for the reaction of sulfur dioxde (SO2) with OH can be found in the JPL data evaluation [Burkholder et al., 2020]:

``` 
    OH_CHEM[0] = 2.9e-31
    OH_CHEM[1] = 4.1
    OH_CHEM[2] = 1.7e-12
    OH_CHEM[3] = -0.2
```

### References

J. B. Burkholder, S. P. Sander, J. Abbatt, J. R. Barker, C. Cappa, J. D. Crounse, T. S. Dibble, R. E. Huie,  C. E. Kolb, M. J. Kurylo, V. L. Orkin, C. J. Percival, D. M. Wilmouth, and P. H. Wine "Chemical Kinetics and Photochemical Data for Use in Atmospheric Studies, Evaluation No. 19," JPL Publication 19-5, Jet Propulsion Laboratory, Pasadena, http://jpldataeval.jpl.nasa.gov, 2020.

Grooß, J.-U. and Russell III, J. M.: Technical note: A stratospheric climatology for O3, H2O, CH4, NOx, HCl and HF derived from HALOE measurements, Atmos. Chem. Phys., 5, 2797-2807, https://doi.org/10.5194/acp-5-2797-2005, 2005.
