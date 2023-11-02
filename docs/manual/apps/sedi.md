# sedi

This application calculates the sedimentation velocity $v_s$ (in units of m/s) of aerosol or cloud particles for a given particle radius $r_p$ and particle density $\rho_p$. In MPTRAC, $v_s$ is calculated for spherical particles following the method described by Jacobson (1999). Sedimentation velocities are valid for particle Reynolds number $Re_p \le 1$, for which $v_s$ is expected to have accuracies better than 10 % (Hesketh, 1996). Larger particles will require additional corrections, which have not been implemented. See Hoffmann et al. (2022, Sect. 2.3.5) for further details.

* Jacobson, M. Z.: Fundamentals of Atmospheric Modeling, Cambridge University Press, https://doi.org/10.1017/CBO9781139165389, 1999.
* Hesketh, H. E. (Ed.): Air Pollution Control: Traditional Hazardous Pollutants, Revised Edition, CRC Press, ISBN 9781566764131, 1996.
* Hoffmann, L., Baumeister, P. F., Cai, Z., Clemens, J., Griessbach, S., Günther, G., Heng, Y., Liu, M., Haghighi Mood, K., Stein, O., Thomas, N., Vogel, B., Wu, X., and Zou, L.: Massive-Parallel Trajectory Calculations version 2.2 (MPTRAC-2.2): Lagrangian transport simulations on graphics processing units (GPUs), Geosci. Model Dev., 15, 2731–2762, https://doi.org/10.5194/gmd-15-2731-2022, 2022.

```
# calling sedi
$./sedi <p> <T> <r_p> <rho_p>
```

The required control parameters are:
* p: Pressure in hPa
* T: Temperature in K
* r_p: particle radius in micrometers
* rho_p: particle density in kg/m^3
