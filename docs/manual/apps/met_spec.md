# met_spec

This app conducts a spectral analysis of the meteorological temperature fields using the fast Fourier transform method. It will provide estimates of zonal wave amplitude, wavelength, and phase as a function of pressure and latitude.

```
Calling met_spec
./met_sepc <ctl> <spec.tab> <met0>
```

The required control parameters are:
* ctl: the control parameter file
* spec.tab: the output file (in form of an ASCII table)
* met0: the metorological input file
