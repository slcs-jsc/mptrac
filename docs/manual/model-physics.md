# Model physics

## Advection

The advection of an air parcel, i.\,e., the position $\vec{x}(t)$ at time $t$ for a given wind and velocity field $\vec{v}(\vec{x}, t)$, is given by the trajectory equation,
$$
\begin{equation}
  \frac{d\vec{x}}{dt}=\vec{v}(\vec{x}, t).
\end{equation}
$$
Here, the position $\vec{x}$ is provided in a meteorological coordinate system, $\vec{x}=(\lambda,\phi,p)$, with longitude $\lambda$, latitude $\phi$, and pressure $p$. The velocity vector $\vec{v}=(u,v,\omega)$ is composed of the zonal wind component $u$, the meridional wind component $v$, and the vertical velocity $\omega$, respectively. Based on its accuracy and computational efficiency for trajectory calculations \citep{roessler18}, we apply the explicit midpoint method to solve the trajectory equation,
$$
\begin{equation}
  \vec{x}(t+\Delta t) = \vec{x}(t) + \Delta t\, \vec{v}\left\{\vec{x}(t) + \frac{\Delta t}{2}\vec{v}\left[\vec{x}(t), t\right], t+\frac{\Delta t}{2}\right\}. \label{eq:midpoint}
\end{equation}
$$

## Turbulent diffusion

## Subgrid-scale wind fluctuations

## Convection

## Sedimentation

## Wet deposition

## Dry deposition

## Hydroxyl chemistry

## Exponential decay

## Boundary conditions
