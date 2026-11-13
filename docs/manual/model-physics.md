# Model physics

## Advection

The advection of an air parcel, i.e., the position $\vec{x}(t)$ at
time $t$ for a given wind and velocity field $\vec{v}(\vec{x}, t)$, is
given by the trajectory equation,

$$
\begin{align}
  \frac{d\vec{x}}{dt}=\vec{v}(\vec{x}, t).
\end{align}
$$

Here, the position $\vec{x}$ is provided in a meteorological
coordinate system, $\vec{x}=(\lambda,\phi,p)$, with longitude
$\lambda$, latitude $\phi$, and a vertical coordinate (pressure $p$ or
zeta $\zeta$). The velocity vector $\vec{v}=(u,v,\omega)$ is composed
of the zonal wind component $u$, the meridional wind component $v$,
and the vertical velocity ($\omega$ or $\dot{\zeta}$), respectively.

For special-purpose regional UTM-grid simulations, the horizontal coordinates are in UTM easting and northing in metres.

#### Integration Method

Based on its accuracy and computational efficiency for trajectory
calculations (Rößler et al., 2018), per default we apply the explicit
midpoint method to solve the trajectory equation,

$$
\begin{align}
  \vec{x}(t+\Delta t) = \vec{x}(t) + \vec{v}(\vec{x}(t) + \frac{\Delta t}{2}\vec{v}\left[\vec{x}(t), t\right], t+\frac{\Delta t}{2})  \Delta t.
\end{align}
$$

Finally, MPTRAC can use the Euler method as well. The Euler method has
low accuracy but is the fastest option available.

The integration method is selected with the parameter `ADVECT`.

| ADVECT      | Integration Method                | 
|-------------|-----------------------------------|
| 1           | explicit Euler                    | 
| 2 (default) | mid-point                         |

`DT_MOD` must not exceed the meteorological interval `DT_MET`. Normally,
`DT_MOD` should be much smaller than `DT_MET` to resolve temporal changes in
the winds and limit extrapolation near meteorological time boundaries. This is
an accuracy recommendation, not a mathematical stability requirement.

Periodic convection, mixing, and KPP chemistry intervals must be integer
multiples of `DT_MOD` when enabled; incompatible settings are rejected.
Incompatible meteorological sampling and periodic output intervals produce
warnings. Supported nonpositive interval settings retain their documented
special behavior. MPTRAC aligns the simulation start time to the `DT_MOD` grid
automatically. Fractional `DT_MOD` values can still miss events because
periodic triggers compare floating-point remainders with zero exactly.

At the first model step, `trac` estimates a recommended maximum time step from
the horizontal meteorological grid spacing:

- Geographic grids: `abs(lon[1] - lon[0]) * 111132 / 150` seconds, using
  longitude spacing in degrees and a reference speed of 150 m/s. No
  latitude-dependent distance correction is applied.
- Cartesian/UTM grids: `min(dx, dy) / 50` seconds, using the smallest adjacent
  x and y grid spacings in metres and a reference speed of 50 m/s.

The estimate roughly limits horizontal displacement to one grid cell per step
at the assumed speed. Exceeding it produces one warning, not an error. It is
an accuracy guideline rather than a strict CFL stability criterion; actual
accuracy also depends on wind speed, spatial gradients, temporal variability,
and the integration scheme. Strong gradients or extreme storms may require
smaller steps, which can be assessed with timestep sensitivity studies.

The estimate uses nominal meteorological grid spacing, which may differ from
the effective resolution of the underlying atmospheric model. Interpolating
onto a finer grid does not necessarily increase effective resolution or require
a proportionally smaller model time step.

#### Vertical coordinate and diabatic/kinematic velocities

Besides the pressure coordinate with the pressure tendency as vertical
velocity, MPTRAC applies the hybrid vertical coordinate $\zeta$ with
associated diabatic vertical velocity $\dot\zeta=\frac{d\zeta}{dt}$
for trajectory calculations. The hybrid coordinate $\zeta$ is defined
as shown in the equation below. Near the surface, the hybrid
coordinate $\zeta$ follows the orography in the form of a sigma-like
coordinate $\sigma=\frac{p}{p_s}$, where $p_s$ is the surface
pressure. At higher altitudes, above $\sigma_{r}=0.3$, the zeta
coordinate is smoothly transformed into the potential temperature
$\theta(p)$ (see also Mahowald et al, 2004).

$$
\begin{equation}
    \zeta(p) = \begin{cases}
                \theta(p) & \text{if}\,\sigma \leq \sigma_r \\
                \theta(p) \sin{\left(\frac{\pi}{2}\frac{1-\sigma(p)}{1-\sigma_r}\right)} & \text{if}\,\sigma > \sigma_r
                \end{cases}
\end{equation}
$$

While kinematic velocities are derived from the continuity equation
diabatic velocities are obtained by radiative calculations. Diabatic
velocities are well-suited for stratospheric transport processes,
because of the isentropic nature of transport at this height. However,
in the troposphere, where parameterisations of mixing in the boundary
layer and for convection are needed, most modules rely on a
formulation in pressure coordinates (see e.g. convection, turbulent
diffusion below).  MPTRAC allows a coupled mode, where advection is
performed in hybrid coordinates ($\zeta$) but other modules can be
employed in pressure coordinates.

At the beginning of each active particle time step, MPTRAC clamps positions
outside the meteorological vertical domain to the model top or the local
surface pressure. This prevents invalid initial positions from being mirrored
deep into the model domain. After advection and all position-changing physics
modules, boundary overshoots are reflected at the model top and local surface.
Particles with a zero individual time step are not modified.

## Turbulent diffusion

Rather complex parameterizations of atmospheric diffusivity are
available for the planetary boundary layer. Much less is known on the
diffusivities in the free troposphere and in the stratosphere, being
in the scope of the MPTRAC model. In our model, the effects of
atmospheric diffusion are simulated by adding stochastic perturbations
to the positions $\vec{x}$ of the air parcels at each time step
$\Delta t$ of the model. For either horizontal component $x_i$,

$$
\begin{align}
  x_i(t+\Delta t) = x_i(t) + \sqrt{2\,D_x\,|\Delta t|}\,\xi_i.
\end{align}
$$

For vertical displacement in geometric height, the implementation uses

$$
\Delta z = \sqrt{2\,D_z\,|\Delta t|}\,\xi_z
  + \left(\frac{dD_z}{dz}-\frac{D_z}{H}\right)|\Delta t|,
$$

where $H=1000H_0$ is the assumed atmospheric scale height in metres.
The deterministic term accounts for changes in vertical diffusivity and
the decrease of atmospheric density with height. MPTRAC converts
$\Delta z$ to a pressure increment using the linearized hydrostatic
relation $\Delta p=-p\,\Delta z/H$, before applying vertical boundary
conditions.

This method requires a vector $\vec{\xi}=(\xi_x, \xi_y, \xi_z)$ of
random variates to be drawn from the standard normal distribution for
each air parcel at each time step. The vector of diffusivities
$\vec{d}=(D_x, D_x, D_z)$ is composed of the horizontal diffusivity
$D_x$ and the vertical diffusivity $D_z$. The model allows to specify
$D_x$ and $D_z$ separately for the troposphere and stratosphere as
control parameters. A smooth transition between tropospheric and
stratospheric diffusivities is created within a $\pm 1km$ log-pressure
altitude range around the tropopause. The tropopause's pressure level
is determined by linear interpolation from a monthly mean zonal mean
climatology derived from the NCEP/NCAR Reanalysis 1 project (Kalnay et
al, 1996). Following choices made for the FLEXPART model (Stohl et al,
2005), default values of $D_x=50m^2s^{-1}$ and $D_z=0$ have been
selected for the troposphere and $D_x=0$ and $D_z=0.1m^2s^{-1}$ have
been selected for the stratosphere.

Note that for simulations without diffusion, the diffusivity needs to
be set to zero everywhere. To specify the diffusivity in the
stratosphere and troposphere as desired the following parameters can
be used:

|Parameter|Layer|Direction|Default|
|---|---|---|---|
|TURB_DX_TROP|Troposphere|horizontal|$50~\rm{m^2~s^{-1}}$|
|TURB_DX_STRAT|Stratosphere|horizontal|$0~\rm{m^2~s^{-1}}$|
|TURB_DZ_TROP|Troposphere|vertical|$0~\rm{m^2~s^{-1}}$|
|TURB_DZ_STRAT|Stratosphere|vertical|$0.1~\rm{m^2~s^{-1}}$|
|TURB_PBL_TRANS|PBL transition layer|fraction of PBL depth|$0$|
|TURB_PBL_SCHEME|optional PBL-specific turbulence scheme ($0$=off, $1$=closure)|-|$0$|

If `TURB_PBL_SCHEME=1`, MPTRAC applies a simplified Hanna/FLEXPART-style PBL closure only inside the diagnosed PBL, while the fixed-`K` background diffusion remains active outside the PBL. The current closure uses hard reflection at the surface and at the diagnosed PBL top.

The current implementation interprets surface sensible heat flux in ECMWF
convention, that is, positive downward and negative upward. Accordingly,
negative `SHF` drives unstable daytime conditions in the closure.

## Subgrid-scale wind fluctuations

In addition to turbulent diffusion, the effects of unresolved
subgrid-scale winds, also referred to as mesoscale wind perturbations,
are considered. The starting point for this approach is the separation
of the horizontal wind and vertical velocity vector $\vec{v}$ into the
grid-scale mean $\bar{\vec{v}}$ and the subgrid-scale perturbations
$\vec{v}'$,

$$
\begin{equation}
  \vec{v}=\bar{\vec{v}} + \vec{v}'.
\end{equation}
$$

It is further assumed that the mean wind $\bar{\vec{v}}$ is given by
the coarse-grid meteorological data, whereas the wind perturbations
$\vec{v}'$ need to be modelled. The sub-grid scale wind perturbations
are simulated by means of a Langevin equation,

$$
\begin{align}
  & v_i'(t+\Delta t) = rv_i'(t) + \sqrt{1-r^2} \, f_i\sigma_i \, \xi_i, \\
  & r = 1-2\frac{\Delta t}{\Delta t_{met}}.
\end{align}
$$

Mathematically, this represents a Markov chain or a random walk, which
adds temporally correlated stochastic perturbations to the winds over
time. The degree of correlation depends on the correlation coefficient
$r$, and therefore on the time step $\Delta t$ of the model and the
time step $\Delta t_{met}$ of the meteorological data. The stochastic
increment is proportional to the grid-scale standard deviation
$\sigma_i$ of the wind or vertical velocity and to a user-defined
scaling factor $f_i$ for the corresponding component. In MPTRAC these
scaling factors are specified directly through the control parameters
`TURB_MESOX` (horizontal) and `TURB_MESOZ` (vertical). Their default
value is 16%, following a choice made for the FLEXPART model (Stohl et
al., 2005). For each air parcel, the grid-scale standard deviation is
calculated from the neighbouring eight grid boxes and two time steps of
the meteorological data. As before, $\vec{\xi}$ is a vector of random
variates drawn from the standard normal distribution.

## Inter-parcel mixing

Lagrangian air parcels do not exchange their carried quantities unless an
explicit mixing process is applied. The inter-parcel mixing module groups
active particles into configurable geographic longitude, latitude, and
altitude boxes
and relaxes each supported quantity towards its local box mean. The update
for particle $i$ is

$$
q_i^\prime = q_i + \alpha_i (\overline{q}_\alpha-q_i),
$$

where the mixing-parameter-weighted target mean

$$
\overline{q}_\alpha =
\frac{\sum_i \alpha_i q_i}{\sum_i \alpha_i}
$$

conserves the unweighted sum of the mixed quantity within each box and
ensemble. Air-mass-weighted conservation requires additional assumptions,
such as equal represented air masses. Particles with zero individual time
steps do not participate. Boxes containing only one particle are unaffected.
Consequently, the effective mixing also depends on particle density and box
size.

Users specify e-folding times rather than per-call exchange fractions. The
tropospheric weight $w_i$ is used to interpolate the mixing rates,

$$
k_i = \frac{w_i}{\tau_\mathrm{trop}} +
      \frac{1-w_i}{\tau_\mathrm{strat}},
\qquad
\alpha_i = 1-\exp(-k_i\,\Delta t_\mathrm{mix}).
$$

For a fixed box mean, a tracer anomaly is reduced to $1/e$ (about 36.8 %)
of its initial value after one e-folding time; about 63.2 % of the anomaly
has then been removed. Mixing is enabled only when both
`MIXING_TAU_TROP` and `MIXING_TAU_STRAT` are positive. Both values are
given in seconds. Their defaults are `-999`, so inter-parcel mixing is
disabled by default. `MIXING_DT` controls the invocation interval and
defaults to 3600 s. If it is non-positive, mixing is applied every model
time step and `DT_MOD` is used as the interval in the exponential update. A
positive `MIXING_DT` must be an integer multiple of `DT_MOD`.

The default mixing grid covers the globe from -5 to 85 km with
$1^\circ \times 1^\circ \times 1$ km boxes. A useful starting point for
sensitivity studies on this grid is

```text
MIXING_DT = 3600
MIXING_TAU_TROP = 2592000    # 30 days
MIXING_TAU_STRAT = 31557600  # 1 year
```

These are example values rather than universal atmospheric constants. Shorter
e-folding times produce faster relaxation, whereas larger boxes mix particles
over larger spatial scales. The mixing times, box dimensions, and particle
density should therefore be varied together in sensitivity studies.

Mixing currently requires geographic coordinates; Cartesian/UTM mixing is
unsupported. With MPI domain decomposition, box means use only particles on
the local rank. Subdomain boundaries follow the meteorological grid and may
split mixing boxes. Particle migration occurs before mixing in each step.

## Convection

MPTRAC's Extreme Convection Parameterization (ECP) represents unresolved
convection by randomly redistributing particles vertically in a column.
It uses convective available potential energy (CAPE) and equilibrium-level
(EL) pressure from meteorological fields. Convective inhibition (CIN) is
an optional additional filter.

At each scheduled evaluation, ECP acts on an active particle at or below
the EL when its interpolated CAPE is at least `CONV_CAPE`, the EL pressure is
finite and between zero and surface pressure, and, if `CONV_CIN > 0`,
its interpolated CIN is finite and at most `CONV_CIN`. A nonpositive
`CONV_CIN` disables the CIN filter.

During a convective event, particle pressure is redrawn between surface
and EL pressure using the established density-weighted method for a
well-mixed column. The routine does not directly change particle number
or tracer mass. Particles with zero individual timestep do not participate.

In subsequent time steps of the model, the trajectories of air parcels
are continued from their new vertical positions, assigned during the
convective mixing event. This dynamic approach enables the model to
capture the transport implications of convective events, while also
accounting for the vertical mixing effects they induce.

The ECP method within MPTRAC is versatile and accommodates different
simulation scenarios based on the chosen threshold $\sf{CAPE_0}$. By
setting $\sf{CAPE_0}$ to zero, the model implements the "extreme
convection" approach, where convection is simulated wherever CAPE meets
the zero threshold below a valid EL. This represents an upper-limit
scenario for the effects of unresolved convection in the meteorological
fields. Turning off ECP leaves explicitly resolved convective
updrafts and any other enabled vertical processes. Intermediate states
can be simulated by selecting specific values for threshold
$\sf{CAPE_0}$. The parameter $\sf{CIN_0}$ can be used to prevent
triggered convection in regions where a larger lower-level inversion
layer is located, with high CAPE above it.

With `CONV_DT <= 0`, convection is evaluated every model timestep. A
positive `CONV_DT` selects steps whose absolute model time is an exact
multiple of that interval. Choose `CONV_DT` as an integer multiple of
`DT_MOD`; incompatible positive intervals are rejected.

Separately, `CONV_MIX_PBL = 1` enables stochastic pressure redistribution
between the surface and the PBL top. `CONV_PBL_TRANS` extends the upper
boundary by a fraction of the PBL pressure thickness above that top.
This option is separate from the PBL Langevin diffusion scheme and does
not require CAPE or EL. If both convection options are enabled, the
lower-pressure upper boundary is used.

For netCDF and GRIB input, `MET_CAPE = 1` computes CAPE, CIN and EL from
meteorological profiles on geographic grids. Binary meteorological input
contains these fields. With `MET_CAPE = 0`, the netCDF and GRIB readers do
not load EL, so ECP cannot use that combination; PBL-only redistribution
remains available.

Convection is off by default.

|Parameter|Explanation|Default Value|
|---|---|---|
|CONV_CAPE|CAPE threshold [J/kg]; negative disables ECP|-999|
|CONV_CIN|Maximum CIN when positive [J/kg]; nonpositive disables the filter|-999|
|CONV_DT|Evaluation interval [s]; nonpositive evaluates every step|-999|
|CONV_MIX_PBL|Enable PBL pressure redistribution|0|
|CONV_PBL_TRANS|Fractional extension above the PBL top|0|

## Sedimentation

In order to take into account the gravitational settling of particles,
the sedimentation velocity $v_s$ needs to be calculated. Using $v_s$,
the change of the vertical position of the particles over the model
time step $\Delta t$ can be calculated. In MPTRAC, $v_s$ is calculated
for spherical particles following the method described by [Jacobson
(1999)](https://www.cambridge.org/core/books/fundamentals-of-atmospheric-modeling/A6B866737D682B17EE46F8449F76FB2C). In
the first step, we calculate the density of dry air,

$$
\begin{equation}
  \rho=\frac{p}{RT},
\end{equation}
$$

using the atmospheric pressure $p$ in Pa, temperature $T$, and the specific
gas constant $R$ of dry air. Particle pressures are stored in hPa and
converted to Pa for this calculation. The dynamic viscosity of air,

$$
\begin{equation}
  \eta=1.8325 \cdot 10^{-5} \frac{\rm{kg}}{\rm{m}~\rm{s}} \left (\frac{416.16~\sf{K}}{T+120~\rm{K}} \right ) \left (\frac{T}{296.16~\rm{K}} \right )^{1.5},
  \end{equation}
$$

and the thermal velocity of an air molecule,

$$
\begin{equation}
  v=\sqrt{\frac{8k_B T}{\pi m_{\rm air}}}=\sqrt{\frac{8RT}{\pi}},
\end{equation}
$$

where $k_B$ is Boltzmann's constant and $m_{\rm air}$ is the mass of an air
molecule, are used to calculate the mean free path $\lambda$,
as well as the Knudsen number $K_n$ for air,

$$
\begin{equation}
  \lambda=\frac{2\eta}{\rho v} \quad \text{and} \quad K_n=\frac{\lambda}{r_p},
\end{equation}
$$

where $r_p$ refers to the particle radius. The Cunningham slipflow
correction is calculated from

$$
\begin{equation}
  G=1+K_n \left [1.249+0.42 \exp \left ( \frac{-0.87}{K_n} \right ) \right]
\end{equation}
$$
 
Finally, the sedimentation velocity is obtained by means of Stokes law
and from the slip-flow correction,

$$
\begin{equation}
 v_s=\frac{2r_p^2g(\rho_p-\rho)}{9\eta}G
\end{equation}
$$

with particle density $\rho_p$ and conventional standard gravitational
acceleration $g$. The radius $r_p$ is specified in micrometres and
$\rho_p$ in kg/m³. Both can be specified individually for each air parcel;
a larger set of parcels can represent a size distribution. Invalid,
nonpositive, or nonfinite radius or density gives zero settling velocity.
The spherical-particle Stokes–Cunningham approximation is intended for
small particle Reynolds numbers, approximately $Re_p \leq 1$.

The signed particle timestep $\Delta t$ changes pressure by
$\Delta p=p v_s\Delta t/(1000 H_0)$, where $p$ is in hPa, $v_s$ in m/s,
and the fixed scale height $H_0=7$ km. This linearized approximation
increases pressure for positive settling velocity in forward time.
Sedimentation does not directly change tracer mass; dry deposition is
handled separately.

## Wet deposition

Wet deposition removes trace gases and aerosol particles within or below
clouds. MPTRAC calculates it in four steps:

(1) A parcel must be below the diagnosed cloud top. Cloud-top and
cloud-bottom pressures $p_{ct}$ and $p_{cb}$ are diagnosed from cloud
cover and cloud liquid, rain, ice, or snow water content.
 
(2) Column-integrated condensate $c_l$ (kg m$^{-2}$) provides an
empirical precipitation proxy $I_s$ (mm h$^{-1}$), obtained by inverting
the relationship reported by Pisso et al. (2019). The default values are
$a_{\rm pre}=0.5$ kg m$^{-2}$(mm h$^{-1}$)$^{-b_{\rm pre}}$
and dimensionless $b_{\rm pre}=0.36$:

$$
\begin{equation}
 I_s=\left(\frac{c_l}{a_{\rm pre}}\right)^{1/b_{\rm pre}}.
\end{equation}
$$

(3) Positive local cloud liquid, rain, ice, or snow water content selects
in-cloud scavenging; otherwise below-cloud scavenging is used. Wet
deposition requires $I_s\geq0.01$ mm h$^{-1}$.

(4) Each branch independently selects empirical aerosol scavenging if
$A>0$, otherwise Henry-law gas scavenging if $H^\ominus>0$. An
unconfigured branch has zero scavenging. The rates are

$$
\begin{equation}
 \Lambda_{\rm aerosol}=A I_s^B\eta,\qquad
 \Lambda_{\rm gas}=\frac{H(T)RT\eta}{\Delta z_c}
                    \frac{I_s}{3.6\times10^6},
\end{equation}
$$

where $A$ has units s$^{-1}$(mm h$^{-1}$)$^{-B}$, $B$ and the retention
factor $\eta$ are dimensionless, $H$ has units mol m$^{-3}$ Pa$^{-1}$,
$R$ is the universal gas constant, and $\Delta z_c$ is the cloud-layer
depth in metres, calculated from $p_{ct}$ and $p_{cb}$. The factor
$3.6\times10^6$ converts mm h$^{-1}$ to m s$^{-1}$. Henry's constant is

$$
\begin{equation}
  H(T)=H^\ominus\exp\left[C\left(\frac{1}{T}-\frac{1}{T^\ominus}\right)\right],
\end{equation}
$$

with $C=-\Delta_{\rm sol}H/R$ in kelvin and
$T^\ominus=298.15$ K. Selected species defaults are listed below;
Sander (2023) provides a current Henry-law compilation. Setting
`WET_DEPO_SO2_PH` to a positive pH applies an SO$_2$ dissociation
correction to in-cloud Henry-law scavenging only. For negative individual
timesteps, MPTRAC mathematically reverses the mass update
$m_{\rm new}=m_{\rm old}\exp(-\Lambda\Delta t)$.


| Species      | $H^\ominus$ (mol m$^{-3}$ Pa$^{-1}$ at 298.15 K) | $C$ (K) |
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

Dry deposition reduces particle-associated tracer mass or volume mixing
ratio; sedimentation separately changes particle position. Set
`DRY_DEPO_MODE=0` (the default) to disable dry deposition, `1` to use the
prescribed velocity `DRY_DEPO_VDEP` (m/s), or `2` to calculate an aerosol
velocity from particle radius `rp` and density `rhop` using `sedi()`.
`DRY_DEPO_VDEP` does not activate deposition and is ignored in mode 2.
Both active modes use the lowermost `DRY_DEPO_DZ` metres above the local
surface (default 30 m). This replaces the former pressure-thickness
parameter `DRY_DEPO_DP`; the two layers and their results differ.

In mode 2, the aerosol velocity $v_{dep}$ is calculated from particle
pressure, interpolated atmospheric temperature, radius, and density as
described in Hoffmann et al. (2022). Nonpositive calculated velocities
produce no dry-deposition loss. In mode 1, $v_{dep}$ is the prescribed
constant velocity, including when aerosol quantities are present.

For both modes, the loss uses the deposition velocity $v_{dep}$, each
particle's signed time step $\Delta t$, and the surface layer depth
$\Delta z=\texttt{DRY\_DEPO\_DZ}$:

$$
\begin{equation}
  m(t+\Delta t)=m(t) \mathrm{exp} \left ( -\frac{\Delta t v_{dep}}{\Delta z} \right )
\end{equation}
$$

For negative time steps, the exponential mathematically reverses the
forward loss and may increase airborne mass.
The layer boundary is evaluated in pressure coordinates with the fixed
scale height $H_0$: $p_{\mathrm{limit}}=p_s\exp[-\Delta z/(1000H_0)]$,
where $H_0$ is in km, $\Delta z$ is in m, and $p_s$ is local surface
pressure. This is layer-averaged loss, not explicit surface contact.

## Hydroxyl chemistry

The hydroxyl radical (OH) is an important oxidant in the atmosphere,
causing the decomposition of many gas-phase species. The oxidation of
different gas-phase species with OH can be classified into two main
categories, bimolecular reactions (e.g., reactions of $\sf{CH_4}$ or
$\sf{NH_3}$), and termolecular reactions (e.g., CO or $\sf{SO_2}$).

For bimolecular reactions, the rate constant is calculated from
Arrhenius law,

$$
\begin{equation}
  k(T)=A \times \mathrm{exp} \left ( -\frac{E}{RT} \right )
\end{equation}
$$

with Avogadro constant $N_A$. For the calculation of the second-order
rate constant k see Eq. 25 in Hoffmann et al. (2022).

For the calculation of k the low- and high-pressure limits of the
reaction rate constant are given by

$$
\begin{equation}
  k_0=k_0^{298} \left ( \frac{T}{298} \right )^{-n}, \quad
  k_\infty=k_\infty^{298} \left ( \frac{T}{298}\right )^{-m} .
\end{equation}
$$

The constants $k_0^{298}$ and $k_\infty^{298}$ at the reference
temperature of 298 K and the exponents n and m need to be specified as
control parameters. The exponents can be set to zero in order to
neglect the temperature dependence of the low- or high pressure limits
of $k_0$ and $k_\infty$. The termolecular reaction rate coefficients
implemented directly into MPTRAC are listed in the table.

| Reaction | A factor | E/R |
| -------- | -------- | --- |
$\sf{CH}_4+\sf{OH} \to \sf{CH}_3+\rm{H_2O}$ | 2.45 $\cdot 10^{-12}$ | 1775 |
$\sf{NH}_3+\sf{OH} \to \sf{H_2O}+\rm{NH}_2$ | 1.7 $\cdot 10^{-12}$ | 710 |
$\sf{O}_3+\sf{OH} \to \sf{HO}_2+\rm{O}_2$ | 1.7 $\cdot 10^{-12}$ | 940 |

Where A is in $\sf{cm^{-3}}\sf{molec^{-1}} \sf{s^{-1}}$ and E/R in K.

Based on the bimolecular reaction rate k=k(T) or the termolecular
reaction rate k=k(T, [M]), the loss of mass of the gas-phase species
over time is calculated from

$$
\begin{equation}
m(t+\Delta t)=m(t) \mathrm{exp} (-k[OH] \Delta t).
\end{equation}
$$

The hydroxyl radical concentrations [OH] are obtained by bilinear
interpolation from the monthly mean zonal mean climatology of Pommrich
et al. (2014). This approach is suitable for global simulations
covering at least several days, as hydroxyl concentrations may vary
significantly between day and nighttime as well as the local
atmospheric composition.

## Exponential decay

A rather generic module was implemented in MPTRAC, to simulate the
loss of mass of an air parcel over a model time step $\Delta t$ due to
any kind of exponential decay process, e.g., chemical loss or
radioactivity,

$$
\begin{equation}
m(t+\Delta t)=m(t)\mathrm{exp} \left ( -\frac{\Delta t}{t_e} \right ).
\end{equation}
$$

The e-folding lifetime $t_e$ of the species needs to be specified as a
control parameter. As typical lifetimes may differ, we implemented an
option to specify separate lifetimes for the troposphere and
stratosphere. A smooth transition between the tropospheric and
stratospheric lifetime is created within a 1 km log-pressure altitude
range around the tropopause.

### Radionuclide decay and deposition

For the explicitly supported radionuclides, `RADIO_DECAY` controls decay in
both airborne activity and deposited ground inventories. The half-lives used
for both phases are 22.3 years for Pb-210, 53.22 days for Be-7, 30.05 years
for Cs-137, and 8.02 days for I-131.

The activity $A$ of a parcel is reduced over a model time step according to

$$
A(t+\Delta t)=A(t)\exp[-(k_\mathrm{dry}+k_\mathrm{wet})\Delta t].
$$

This equation describes deposition loss only. Radioactive decay is applied
separately when `RADIO_DECAY=1`, using the half-lives defined in `mptrac.h`.
Rn-222 and Xe-133 undergo airborne decay but are not included in the
deposition inventories.

Radioactive dry deposition uses the same `DRY_DEPO_DZ` metre layer above
the local surface, independently of `DRY_DEPO_MODE`. For its thickness
$\Delta z$, the rate is
$k_\mathrm{dry}=v_\mathrm{d}/\Delta z$. The fixed deposition velocities are
0.002 m s$^{-1}$ for Pb-210 and Cs-137, 0.001 m s$^{-1}$ for Be-7, and
0.0005 m s$^{-1}$ for aerosol-bound I-131.

Wet deposition is applied below the diagnosed cloud top when the precipitation
proxy

$$
I_s=\left(\frac{c_l}{p_0}\right)^{1/p_1}
$$

is at least 0.01. Here, $c_l$ is total column cloud water and $p_0$ and $p_1$
are `WET_DEPO_PRE[0]` and `WET_DEPO_PRE[1]`. The wet rates are
$10^{-4} I_s \eta$ s$^{-1}$ for Pb-210 and Cs-137 and
$5\times10^{-5} I_s \eta$ s$^{-1}$ for Be-7 and I-131. The retention factor
$\eta$ uses `WET_DEPO_IC_RET_RATIO` in cold in-cloud conditions and
`WET_DEPO_BC_RET_RATIO` below cold clouds, with temperature interpolation for
mixed-phase clouds.

The half-lives, dry deposition velocities, and wet deposition coefficients are
defined as `RADIO_HALF_LIFE_*`, `RADIO_DRY_VDEP_*`, and
`RADIO_WET_COEFF_*` macros in `mptrac.h`. Each definition is guarded by
`#ifndef` and can be replaced at compile time, for example with
`make DEFINES="-DRADIO_DRY_VDEP_CS137=0.0015"`.

Deposited activities are accumulated on the grid selected by the `GRID_*`
parameters and can be written with the `DEPO_*` output parameters. Deposited
ground inventories are not contained in atmospheric restart files. A run
started from such a restart therefore begins with zero ground inventory; only
deposition occurring after the restart is reported.

The radionuclide deposition module currently requires a latitude/longitude
meteorological grid (`MET_COORD_TYPE=0`). Domain-decomposed simulations are
not yet supported; enabling `RADIO_DEPO` in a DD build terminates with an
error rather than producing incomplete ground inventories.

## Boundary conditions

When an air parcels reach the upper and lower boundary layer two
options are available. First, they can be reflected back or second,
their position is set to the lowest or highest height level available in the
meteorological data and hence will slide along the boundaries until
updrafts or downdrafts transport them back into the wider model
domain. Those two options can be selected with the parameter REFLECT,
which is 1 if air parcels are supposed to be reflected and 0 if they
are supposed to be set to the boundary heights. The default value is
0, hence no reflection.
