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
|TURB_DX_TROP|Troposphere|horizontal|$50~\rm{m^2~s^{-1}}$|
|TURB_DX_STRAT|Stratosphere|horizontal|$0~\rm{m^2~s^{-1}}$|
|TURB_DZ_TROP|Troposphere|vertical|$0~\rm{m^2~s^{-1}}$|
|TURB_DZ_STRAT|Stratosphere|vertical|$0.1~\rm{m^2~s^{-1}}$|


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

Mathematically, this represents a Markov chain or a random walk, which adds temporally correlated stochastic perturbations to the winds over time. The degree of correlation depends on the correlation coefficient $r$, and therefore on the time step $\Delta t$ of the model and the time step $\Delta t_{met}$ of the meteorological data. The variance $(f\sigma_i)^2$ of the random component added at each time step depends on the grid-scale variance $\sigma_i^2$ of the wind and velocity data and a scaling factor $f$, which is used for downscaling to the subgrid scale. The scaling factor $f$ needs to be specified as a control parameter of the model. The default value is $f=40$%, following a choice made for the FLEXPART model (Stohl et al., 2005). For each air parcel, the grid-scale variance $\sigma_i^2$ is calculated from the neighbouring eight grid boxes and two-time step of the meteorological data. To make computations more efficient, the grid-scale variance of each parcel is kept in cache. As before, $\vec{\xi}$ is a vector of random variates to be drawn from the standard normal distribution.

The scaling factor f is *not* directly set in MPTRAC. Instead the square $f^2$ can be set for the vertical and horizontal dispersion of air parcels with the parameters TURB_MESOX and TURB_MESOZ. Accordingly their default value is 16%.

## Convection

The convection parameterization implemented in MPTRAC introduces the concept of the Extreme Convection Parametrization (ECP) to efficiently simulate convective effects within Lagrangian transport models. This method is based on the assumption that convective events create well-mixed columns of air, and it leverages convective available potential energy (CAPE), convective inhibition (CIN) and equilibrium level (EL) data from meteorological fields to capture this behavior.

The fundamental idea underlying the ECP approach involves several key steps. First, the meteorological fields are used to interpolate gridded CAPE, CIN and EL values to the horizontal positions of individual air parcels within the Lagrangian model. Then, user-defined thresholds of $\sf{CAPE_0}$ and $\sf{CIN_0}$ are applied globally. If the interpolated CAPE value of an air parcel surpasses the $\sf{CAPE_0}$ threshold or the interpolated CIN value is below the $\sf{CIN_0}$ threshold and the parcel's position is below the EL, the model triggers a convective event.

During a convective event, air parcels affected by the event are subject to random vertical redistribution within the atmospheric column, spanning from the Earth's surface to the equilibrium level. This vertical redistribution takes into account air density to ensure the creation of a well-mixed column of air. Importantly, this redistribution process adheres to mass conservation principles, maintaining both the number of air parcels and their collective mass constant.

In subsequent time steps of the model, the trajectories of air parcels are continued from their new vertical positions, assigned during the convective mixing event. This dynamic approach enables the model to capture the transport implications of convective events, while also accounting for the vertical mixing effects they induce.

The ECP method within MPTRAC is versatile and accommodates different simulation scenarios based on the chosen threshold $\sf{CAPE_0}$. By setting $\sf{CAPE_0}$ to zero, the model implements the "extreme convection" approach, where convection is simulated wherever CAPE exists below the EL. This represents an upper-limit scenario for the effects of unresolved convection in the meteorological fields. On the other hand, by turning off the ECP entirely, only explicitly resolved convective updrafts of the meteorological fields are considered, representing a lower-limit scenario. Intermediate states can be simulated by selecting specific values for threshold $\sf{CAPE_0}$. The parameter $\sf{CIN_0}$ can be used to prevent triggered convection in regions where a larger lower-level inversion layer is located, with high CAPE above it.

The frequency of applying the ECP can be customized to the simulation's requirements. It can be implemented every model time step or at user-defined intervals that match typical convective timescales. This adaptability adds to the flexibility of the ECP method, making it suitable for a variety of simulation contexts.

The convection parameterisation is set off per default, i.e. values are set to -999. 

|Parameter|Explanation|Default Value|
|---|---|---|
|CONV_CAPE|CAPE threshold|-999|
|CONV_CIN|CIN threshold|-999|
|CONV_DT|evaluation interval|-999|

## Sedimentation

In order to take into account the gravitational settling of particles, the sedimentation velocity $v_s$ needs to be calculated. Using $v_s$, the change of the vertical position of the particles over the model time step $\Delta t$ can be calculated. In MPTRAC, $v_s$ is calculated for spherical particles following the method described by [Jacobson (1999)](https://www.cambridge.org/core/books/fundamentals-of-atmospheric-modeling/A6B866737D682B17EE46F8449F76FB2C). In the first step, we calculate the density of dry air,

$$
\begin{equation}
  \rho=\frac{p}{RT},
\end{equation}
$$

using the atmospheric pressure $p$, temperature $T$, and the specific gas constant $R$ of dry air. The dynamic viscosity of air, 

$$
\begin{equation}
  \eta=1.832515 \cdot 10^{-5} \frac{\rm{kg}}{\rm{m}~\rm{s}} \left (\frac{416.16~\sf{K}}{T+120~\rm{K}} \right ) \left (\frac{T}{296.16~\rm{K}} \right )^{1.5},
  \end{equation}
$$

and the thermal velocity of an air molecule,

$$
\begin{equation}
  v=\sqrt{\rho RT}, 
\end{equation}
$$

are used to calculate the mean free path of an air molecule $\lambda$, as well as the Knudsen number $K_n$ for air,

$$
\begin{equation}
  \lambda=\frac{\rho}{2\eta v} \quad \text{and} \quad K_n=\frac{\lambda}{r_p},
\end{equation}
$$

where $r_p$ refers to the particle radius. The Cunningham slipflow correction is calculated from

$$
\begin{equation}
  G=1+K \left [1.249+0.42 \exp \left ( \frac{-0.87}{K} \right ) \right]
\end{equation}
$$
 
Finally, the sedimentation velocity is obtained by means of Stokes law and from the slip-flow correction, 

$$
\begin{equation}
 v_s=\frac{2r_p^2g(\rho_p-\rho)}{9\eta}G
\end{equation}
$$

with particle density $\rho_p$ and conventional standard gravitational acceleration g. Note that $r_p$ and $\rho_p$ can be specified individually for each air parcel. A larger set of parcels can be used to represent a size distribution of aerosol or cloud particles.

## Wet deposition

Wet deposition causes the removal of trace gases and aerosol particles from the atmosphere within or below clouds by mixing with suspended water and following washout through rain, snow, or fog. Wet deposition in MPTRAC is calculated based on the following four steps:

(1) it is determined whether an air parcel is located below a cloud top. The cloud-top pressure $p_{ct}$ is determined from the meteorological data as the highest vertical level where cloud water or ice (i.e., CLWC, CRWC, CIWC, or CSWC) is existent. 
 
(2) the wet deposition parametrization determines an estimate of the subgrid-scale precipitation rate $I_s$, which is needed to calculate the scavenging coefficient $\Lambda$. The precipitation rate $I_s$ (in units of $\rm{mm h^{-1}}$) is calculated from the total column cloud water $c_l$ (in units of $\rm{kg m^{-2}}$) by means of a correlation function reported by Pisso et al. (2019),

$$
\begin{equation}
 I_s= (2c_l)^{1/0.36},
\end{equation}
$$

(3) it is inferred whether the air parcel is located within or below the cloud because scavenging coefficients will be different under these conditions. The position of the air parcel within or below the cloud is determined by interpolating the cloud water content to the position of the air parcel and by testing whether the interpolated values are larger than zero.

(4) the scavenging coefficient $\Lambda$ is calculated based on the precipitation rate $I_s$,

$$
\begin{equation}
 \Lambda=HRT I_s \Delta z_c^{-1},
\end{equation}
$$

where H is Henry’s law constant, R is the universal gas constant, and $\Delta z_c$ is the depth of the cloud layer, which is calculated from $p_{ct}$ and $p_{cb}$. Henry’s law constant is obtained from

$$
\begin{equation}
  H(T)=H^\ominus \exp \left [\frac{\Lambda_{sol}H}{R} \left (\frac{1}{T}-\frac{1}{T^\ominus}  \right ) \right] 
\end{equation}
$$

The constants $H^\ominus$ and $\Delta_{sol}H/R$ with enthalpy of dissolution $\Delta_{sol}H$ at the reference temperature $T^\ominus=298.15~K$ need to be specified as control parameters. Values for a wide range of species are tabulated by Sander (2015). The values of selected species of interest are summarized in the following table are included as default parameters in MPTRAC.


| Species      | $H^\ominus$ (at 298.15 K) | $-\frac{\Delta_{sol}H}{R}$ |
| ------------ | ------------------ | ------|
| $\sf{CF_2Cl_2}$ |  3.0 $\cdot 10^{-5}$ | 3500 |   
| $\sf{CFCl_3}$  | 1.1 $\cdot 10^{-4}$ | 3300 |
| $\sf{CH_4}$    | 1.4 $\cdot 10^{-5}$ | 1600 |
| $\sf{CO}$      | 9.7 $\cdot 10^{-6}$ | 1300 |
| $\sf{CO_2}$    | 3.3 $\cdot 10^{-4}$ | 2400 |
| $\sf{N_2O}$    | 2.4 $\cdot 10^{-4}$ | 2600 |
| $\sf{NH_3}$    | 5.9 $\cdot 10^{-1}$ | 4200 |
| $\sf{HNO_3}$   | 2.1 $\cdot 10^3$    | 8700 |
| $\sf{NO}$      | 1.9 $\cdot 10^{-5}$ | 1600 |
| $\sf{NO_2}$    | 1.2 $\cdot 10^{-4}$ | 2400 |
| $\sf{O_3}$     | 1.0 $\cdot 10^{-4}$ | 2800 |
| $\sf{SF_6}$    | 2.4 $\cdot 10^{-6}$ | 3100 |
| $\sf{SO_2}$    | 1.3 $\cdot 10^{-2}$ | 2900 |

## Dry deposition
Dry deposition leads to a loss of mass of aerosol particles or trace gases by gravitational settling or chemical and physical interactions with the surface of the dry phase. In the parametrization implemented in MPTRAC, dry deposition is calculated for air parcels located in the lowermost $\Delta_p=$ 30 hPa layer above the surface. This corresponds to a layer width of $\Delta z \approx$ 200 m at standard conditions. 

For aerosol particles, the deposition velocity $v_{dep}$ will be calculated as described in Hoffmann et al. (2022) as a function of surface pressure $p$ and temperature $T$ as well as particle radius $r_p$ and particle density $\rho$. For trace gases, the deposition velocity $v_{dep}$ needs to be specified as a control parameter. Currently, this parameter is set to a constant value across the globe for each trace gas. For future applications with a stronger focus on the boundary layer, $v_{dep}$ will need to vary geographically to account for dependence on the surface characteristics and atmospheric conditions. 

For both particles and gases, the loss of mass is calculated based on the deposition velocity $v_{dep}$, the model time step $\Delta t$, and the surface layer width $\Delta z$ from 

$$
\begin{equation}
  m(t+\Delta t)=M(t) \mathrm{exp} \left ( -\frac{\Delta t v_{dep}}{\Delta z} \right )
\end{equation}
$$

## Hydroxyl chemistry

The hydroxyl radical (OH) is an important oxidant in the atmosphere, causing the decomposition of many gas-phase species. The oxidation of different gas-phase species with OH can be classified into two main categories, bimolecular reactions (e.g., reactions of $\sf{CH_4}$ or $\sf{NH_3}$), and termolecular reactions (e.g., CO or $\sf{SO_2}$).

For bimolecular reactions, the rate constant is calculated from Arrhenius law,

$$
\begin{equation}
  k(T)=A \times \mathrm{exp} \left ( -\frac{E}{RT} \right )
\end{equation}
$$

with Avogadro constant $N_A$. For the calculation of the second-order rate constant k see Eq. 25 in Hoffmann et al. (2022).

For the calculation of k the low- and high-pressure limits of the reaction rate constant are given by


$$
\begin{equation}
  k_0=k_0^{298} \left ( \frac{T}{298} \right )^{-n}, \quad
  k_\infty=k_\infty^{298} \left ( \frac{T}{298}\right )^{-m} .
\end{equation}
$$

The constants $k_0^{298}$ and $k_\infty^{298}$ at the reference temperature of 298 K and the exponents n and m need to be specified as control parameters. The exponents can be set to zero in order to neglect the temperature dependence of the low- or high pressure limits of $k_0$ and $k_\infty$. The termolecular reaction rate
coefficients implemented directly into MPTRAC are listed in the table.

| Reaction | A factor | E/R |
| -------- | -------- | --- |
$\sf{CH}_4+\sf{OH} \to \sf{CH}_3+\rm{H_2O}$ | 2.45 $\cdot 10^{-12}$ | 1775 |
$\sf{NH}_3+\sf{OH} \to \sf{H_2O}+\rm{NH}_2$ | 1.7 $\cdot 10^{-12}$ | 710 |
$\sf{O}_3+\sf{OH} \to \sf{HO}_2+\rm{O}_2$ | 1.7 $\cdot 10^{-12}$ | 940 |

Where A is in $\sf{cm^{-3}}\sf{molec^{-1}} \sf{s^{-1}}$ and E/R in K.

Based on the bimolecular reaction rate k=k(T) or the termolecular reaction rate k=k(T, [M]), the loss of mass of the gas-phase species over time is calculated from

$$
\begin{equation}
m(t+\Delta t)=m(t) \mathrm{exp} (-k[OH] \Delta t).
\end{equation}
$$

The hydroxyl radical concentrations [OH] are obtained by bilinear interpolation from the monthly mean zonal mean climatology of Pommrich et al. (2014). This approach is suitable for global simulations covering at least several days, as hydroxyl concentrations may vary significantly between day and nighttime as well as the local atmospheric composition.

## Exponential decay

A rather generic module was implemented in MPTRAC, to simulate the loss of mass of an air parcel over a model time step $\Delta t$ due to any kind of exponential decay process, e.g.,
chemical loss or radioactivity,

$$
\begin{equation}
m(t+\Delta t)=m(t)\mathrm{exp} \left ( -\frac{\Delta t}{t_e} \right ).
\end{equation}
$$

The e-folding lifetime $t_e$ of the species needs to be specified as a control parameter. As typical lifetimes may differ, we implemented an option to specify separate lifetimes for the troposphere and stratosphere. A smooth transition between the tropospheric and stratospheric lifetime is created within a 1 km log-pressure altitude range around the tropopause.

## Boundary conditions

When an air parcels reach the upper and lower boundary layer two options are available. First, they can be reflected back or second, their position is set to the lowest or highest height available in the meteorological data and hence will slide along the boundaries until updrafts or downdrafts transport them back into the wider model domain. Those two options can be selected with the parameter REFLECT, which is 1 if air parcels are supposed to be reflected and 0 if they are supposed to be set to the boundary heights. The default value is 0, hence no reflection.



