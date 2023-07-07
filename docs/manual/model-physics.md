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

## Subgrid-scale wind fluctuations

## Convection

## Sedimentation

## Wet deposition

## Dry deposition

## Hydroxyl chemistry

## Exponential decay

## Boundary conditions
