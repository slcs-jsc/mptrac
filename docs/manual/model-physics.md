# Model physics

## Advection

The advection of an air parcel, i.\,e., the position $\vec{x}(t)$ at time $t$ for a given wind and velocity field $\vec{v}(\vec{x}, t)$, is given by the trajectory equation,

$$
\begin{align}
  \frac{d\vec{x}}{dt}=\vec{v}(\vec{x}, t).
\end{align}
$$

Here, the position $\vec{x}$ is provided in a meteorological coordinate system, $\vec{x}=(\lambda,\phi,p)$, with longitude $\lambda$, latitude $\phi$, and a vertical coordinate (pressure $p$ or zeta $\zeta$). The velocity vector $\vec{v}=(u,v,\omega)$ is composed of the zonal wind component $u$, the meridional wind component $v$, and the vertical velocity ($\omega$ or $\dot{\zeta}$), respectively. 

#### Integration Method
Based on its accuracy and computational efficiency for trajectory calculations (Rößler et al., 2018), per default we apply the explicit midpoint method to solve the trajectory equation,

$$
\begin{align}
  \vec{x}(t+\Delta t) = \vec{x}(t) + \vec{v}(\vec{x}(t) + \frac{\Delta t}{2}\vec{v}\left[\vec{x}(t), t\right], t+\frac{\Delta t}{2})  \Delta t.
\end{align}
$$

As an alternative, the classical 4th-order Runge-Kutta method can be employed for the simulation, which has higher accuracy at higher computational costs. The method is defined by the equations:

$$
\begin{align}
\vec{x}(t+\Delta t) = \vec{x}(t) + \frac{1}{6}\left(k_1 + 2k_2 + 2k_3 + k_4 \right)\Delta t\\
\end{align}
$$

$$
\begin{align}
 \vec{k_1} = \vec{v}(t, \vec{x}(t))\\
 \vec{k_2} = \vec{v}\left(t + \frac{\Delta t}{2}, \vec{x}(t) + \Delta t\frac{k_1}{2}\right)\\
 \vec{k_3} = \vec{v}\left(t + \frac{\Delta t}{2}, \vec{x}(t) + \Delta t\frac{k_2}{2}\right)\\
 \vec{k_4} = \vec{v}\left(t + \Delta t, \vec{x}(t) + \Delta t\,k_3\right)
\end{align}
$$

Finally, MPTRAC can use the Euler method as well. The Euler method has low accuracy but is the fastest option available.  

The integration method is selected with the parameter ADVECT.

| ADVECT | Integration Method | 
|---|---|
|  1 | explicit Euler   | 
|  2 (default) | mid-point  |
|  4 |  classical 4th order Runge-Kutta  |

#### Vertical coordinate and diabatic/kinematic velocities

Besides the pressure coordinate with the pressure tendency as vertical velocity, MPTRAC applies the hybrid vertical coordinate $\zeta$ with associated diabatic vertical velocity $\dot\zeta=\frac{d\zeta}{dt}$ for trajectory calculations. The hybrid coordinate $\zeta$ is defined as shown in the equation below. Near the surface, the hybrid coordinate $\zeta$ follows the orography in the form of a sigma-like coordinate $\sigma=\frac{p}{p_s}$, where $p_s$ is the surface pressure. At higher altitudes, above $\sigma_{r}=0.3$, the zeta coordinate is smoothly transformed into the potential temperature $\theta(p)$ (see also Mahowald et al, 2004). 

$$
\begin{equation}
    \zeta(p) = \begin{cases}
                \theta(p) & \text{if}\,\sigma < \sigma_r \\
                \theta(p) \sin{\left(\frac{1-\sigma(p)}{1-\sigma_r}\right)} & \text{if}\,\sigma \geq \sigma_r
                \end{cases}
\end{equation}
$$

While kinematic velocities are derived from the continuity equation diabatic velocities are obtained by radiative calculations. Diabatic velocities are well-suited for stratospheric transport processes, because of the isentropic nature of transport at this height. However, in the troposphere, where parameterisations of mixing in the boundary layer and for convection are needed, most modules rely on a formulation in pressure coordinates (see e.g. convection, turbulent diffusion below).  MPTRAC allows a coupled mode, where advection is performed in hybrid coordinates ($\zeta$) but other modules can be employed in pressure coordinates.

The vertical coordinate and velocity are selected with the parameter VERT_COORD_AP. The coupling between the different coordinate system is applied with CPL_ZETA_AND_PRESS_MODULES=1.

| VERT_COORD_AP  | CPL_ZETA_AND_PRESS_MODULES | Vertical coordinate | Vertical velocity  | coupling?
|---|---|---|---|---|
|  0 (default)  | 0(default) | $p$ | $\Omega$ kinematic | off|
|  1  | 0 (default) | $\zeta$  | $\dot{\zeta}$ diabatic | off|
|  1   | 1  | $\zeta$ | $\dot{\zeta}$ diabatic| on |

## Turbulent diffusion

Rather complex parameterizations of atmospheric diffusivity are available for the planetary boundary layer. Much less is known on the diffusivities in the free troposphere and in the stratosphere, being in the scope of the MPTRAC model. In our model, the effects of atmospheric diffusion are simulated by adding stochastic perturbations to the positions $\vec{x}$ of the air parcels at each time step $\Delta t$ of the model,

$$
\begin{align}
  \Delta x_i(t+\Delta t) = x_i(t) + \sqrt{2\,D_i\,\Delta t}\xi_i.
\end{align}
$$

This method requires a vector $\vec{\xi}=(\xi_x, \xi_y, \xi_z)$ of random variates to be drawn from the standard normal distribution for each air parcel at each time step. The vector of diffusivities $\vec{d}=(D_x, D_x, D_z)$ is composed of the horizontal diffusivity $D_x$ and the vertical diffusivity $D_z$. The model allows to specify $D_x$ and $D_z$ separately for the troposphere and stratosphere as control parameters. A smooth transition between tropospheric and stratospheric diffusivities is created within a $\pm 1km$ log-pressure altitude range around the tropopause. The tropopause's pressure level is determined by linear interpolation from a monthly mean zonal mean climatology derived from the NCEP/NCAR Reanalysis 1 project (Kalnay et al, 1996). Following choices made for the FLEXPART model (Stohl et al, 2005), default values of $D_x=50m^2s^{-1}$ and $D_z=0$ have been selected for the troposphere and $D_x=0$ and $D_z=0.1m^2s^{-1}$ have been selected for the stratosphere.

Note that for simulations without diffusion, the diffusivity needs to be set to zero everywhere. To specify the diffusivity in the stratosphere and troposphere as desired the following parameters can be used:

|Parameter|Layer|Direction|Default|
|---|---|---|---|
|TURB_DX_TROP|Troposphere|horizontal|$50m^2s^{-1}$|
|TURB_DX_STRAT|Stratosphere|horizontal|$0m^2s^{-1}$|
|TURB_DZ_TROP|Troposphere|vertical|$0m^2s^{-1}$|
|TURB_DZ_STRAT|Stratosphere|vertical|$0.1m^2s^{-1}$|


## Subgrid-scale wind fluctuations

In addition to turbulent diffusion, the effects of unresolved subgrid-scale winds, also referred to as mesoscale wind perturbations, are considered. The starting point for this approach is the separation of the horizontal wind and vertical velocity vector $\vec{v}$ into the grid-scale mean $\bar{\vec{v}}$ and the subgrid-scale perturbations $\vec{v}'$,

$$
\begin{equation}
  \vec{v}=\bar{\vec{v}} + \vec{v}'.
\end{equation}
$$

It is further assumed that the mean wind $\bar{\vec{v}}$ is given by the coarse-grid meteorological data, whereas the wind perturbations $\vec{v}'$ need to be modelled. The sub-grid scale wind perturbations are simulated by
means of a Langevin equation,

$$
\begin{align}
  & v_i'(t+\Delta t) = rv_i'(t) + \sqrt{1-r^2} (f\sigma_i)^2 \xi_i, \\
  & r = 1-2\frac{\Delta t}{\Delta t_{met}}.
\end{align}
$$

Mathematically, this represents a Markov chain or a random walk, which adds temporally correlated stochastic perturbations to the winds over time. The degree of correlation depends on the correlation coefficient $r$, and therefore on the time step $\Delta t$ of the model and the time step $\Delta t_{met}$ of the meteorological data. The variance $(f\sigma_i)^2$ of the random component added at each time step depends on the grid-scale variance $\sigma_i^2$ of the wind and volcitiy data and a scaling factor $f$, which is used for downscaling to the subgrid scale. The scaling factor $f$ needs to be specified as a control parameter of the model. The default value is $f=40$%, following a choice made for the FLEXPART model (Stohl et al., 2005). For each air parcel, the grid-scale variance $\sigma_i^2$ is calculated from the neighbouring eight grid boxes and two-time step of the meteorological data. To make computations more efficient, the grid-scale variance of each parcel is kept in cache. As before, $\vec{\xi}$ is a vector of random variates to be drawn from the standard normal distribution.

The scaling factor f is *not* directly set in MPTRAC. Instead the square $f^2$ can be set for the vertical and horizontal dispersion of air parcels with the parameters TURB_MESOX and TURB_MESOZ. Accordingly their default value is 16%.

## Convection

The convection parameterization implemented in MPTRAC introduces the concept of the Extreme Convection Parametrization (ECP) to efficiently simulate convective effects within Lagrangian transport models. This method is based on the assumption that convective events create well-mixed columns of air, and it leverages convective available potential energy (CAPE), convective inhibition (CIN) and equilibrium level (EL) data from meteorological fields to capture this behavior.

The fundamental idea underlying the ECP approach involves several key steps. First, the meteorological fields are used to interpolate gridded CAPE, CIN and EL values to the horizontal positions of individual air parcels within the Lagrangian model. Then, user-defined thresholds of CAPE$_0$ and CIN$_0$ are applied globally. If the interpolated CAPE value of an air parcel surpasses the CAPE$_0$ threshold or the interpolated CIN value is below the CIN$_0$ threshold and the parcel's position is below the EL, the model triggers a convective event.

During a convective event, air parcels affected by the event are subject to random vertical redistribution within the atmospheric column, spanning from the Earth's surface to the equilibrium level. This vertical redistribution takes into account air density to ensure the creation of a well-mixed column of air. Importantly, this redistribution process adheres to mass conservation principles, maintaining both the number of air parcels and their collective mass constant.

In subsequent time steps of the model, the trajectories of air parcels are continued from their new vertical positions, assigned during the convective mixing event. This dynamic approach enables the model to capture the transport implications of convective events, while also accounting for the vertical mixing effects they induce.

The ECP method within MPTRAC is versatile and accommodates different simulation scenarios based on the chosen threshold CAPE$_0$. By setting CAPE$_0$ to zero, the model implements the "extreme convection" approach, where convection is simulated wherever CAPE exists below the EL. This represents an upper-limit scenario for the effects of unresolved convection in the meteorological fields. On the other hand, by turning off the ECP entirely, only explicitly resolved convective updrafts of the meteorological fields are considered, representing a lower-limit scenario. Intermediate states can be simulated by selecting specific values for threshold CAPE$_0$. The parameter CIN$_0$ can be used to prevent triggered convection in regions where a larger lower-level inversion layer is located, with high CAPE above it.

The frequency of applying the ECP can be customized to the simulation's requirements. It can be implemented every model time step or at user-defined intervals that match typical convective timescales. This adaptability adds to the flexibility of the ECP method, making it suitable for a variety of simulation contexts.

The convection parameterisation is set off per default, i.e. values are set to -999. 

|Parameter|Explanation|Default Value|
|---|---|---|
|CONV_CAPE|CAPE threshold|-999|
|CONV_CIN|CIN threshold|-999|
|CONV_DT|evaluation interval|-999|

## Sedimentation

## Wet deposition

## Dry deposition

## Hydroxyl chemistry

## Exponential decay

## Boundary conditions

When an air parcels reach the upper and lower boundary layer two options are available. First, they can be reflected back or second, their position is set to the lowest or highest height available in the meteorological data and hence will slide along the boundaries until updrafts or downdrafts transport them back into the wider model domain. Those two options can be selected with the parameter REFLECT, which is 1 if air parcels are supposed to be reflected and 0 if they are supposed to be set to the boundary heights. The default value is 0, hence no reflection.



