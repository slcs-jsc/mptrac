/*
  This file is part of MPTRAC.
  
  MPTRAC is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  MPTRAC is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with MPTRAC. If not, see <http://www.gnu.org/licenses/>.
  
  Copyright (C) 2013-2021 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  MPTRAC library declarations.
  
  \mainpage
  
  Massive-Parallel Trajectory Calculations (MPTRAC) is a Lagrangian
  particle dispersion model for the troposphere and stratosphere.
  
  This reference manual provides information on the algorithms
  and data structures used in the code. Further information can be found at:

  https://github.com/slcs-jsc/mptrac
*/

/* ------------------------------------------------------------
   Includes...
   ------------------------------------------------------------ */

#include <ctype.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_statistics.h>
#include <math.h>
#include <netcdf.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#ifdef MPI
#include "mpi.h"
#endif

#ifdef _OPENACC
#include "openacc.h"
#include "curand.h"
#endif

/* ------------------------------------------------------------
   Constants...
   ------------------------------------------------------------ */

/*! Standard gravity [m/s^2]. */
#define G0 9.80665

/*! Scale height [km]. */
#define H0 7.0

/*! Boltzmann constant [kg m^2/(K s^2)]. */
#define KB 1.3806504e-23

/*! Molar mass of dry air [g/mol]. */
#define MA 28.9644

/*! Molar mass of water vapor [g/mol]. */
#define MH2O 18.01528

/*! Molar mass of ozone [g/mol]. */
#define MO3 48.00

/*! Standard pressure [hPa]. */
#define P0 1013.25

/*! Standard temperature [K]. */
#define T0 273.15

/*! Specific gas constant of dry air [J/(kg K)]. */
#define RA 287.058

/*! Ideal gas constant [J/(mol K)]. */
#define RI 8.3144598

/*! Mean radius of Earth [km]. */
#define RE 6367.421

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Maximum length of ASCII data lines. */
#define LEN 5000

/*! Maximum number of atmospheric data points. */
#define NP 10000000

/*! Maximum number of quantities per data point. */
#define NQ 12

/*! Maximum number of pressure levels for meteorological data. */
#define EP 112

/*! Maximum number of longitudes for meteorological data. */
#define EX 1201

/*! Maximum number of latitudes for meteorological data. */
#define EY 601

/*! Maximum number of longitudes for gridded data. */
#define GX 720

/*! Maximum number of latitudes for gridded data. */
#define GY 360

/*! Maximum number of altitudes for gridded data. */
#define GZ 100

/*! Maximum number of data points for ensemble analysis. */
#define NENS 2000

/*! Maximum number of OpenMP threads. */
#define NTHREADS 512

/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

/*! Allocate and clear memory. */
#define ALLOC(ptr, type, n)				 \
  if((ptr=calloc((size_t)(n), sizeof(type)))==NULL)      \
    ERRMSG("Out of memory!");

/*! Convert degrees to zonal distance. */
#define DEG2DX(dlon, lat)					\
  ((dlon) * M_PI * RE / 180. * cos((lat) / 180. * M_PI))

/*! Convert degrees to meridional distance. */
#define DEG2DY(dlat)				\
  ((dlat) * M_PI * RE / 180.)

/*! Convert pressure change to vertical distance. */
#define DP2DZ(dp, p)				\
  (- (dp) * H0 / (p))

/*! Convert zonal distance to degrees. */
#define DX2DEG(dx, lat)						\
  (((lat) < -89.999 || (lat) > 89.999) ? 0			\
   : (dx) * 180. / (M_PI * RE * cos((lat) / 180. * M_PI)))

/*! Convert meridional distance to degrees. */
#define DY2DEG(dy)				\
  ((dy) * 180. / (M_PI * RE))

/*! Convert vertical distance to pressure change. */
#define DZ2DP(dz, p)				\
  (-(dz) * (p) / H0)

/*! Compute Cartesian distance between two vectors. */
#define DIST(a, b) sqrt(DIST2(a, b))

/*! Compute squared distance between two vectors. */
#define DIST2(a, b)                                                     \
  ((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))

/*! Compute dot product of two vectors. */
#define DOTP(a, b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

/*! Print error message and quit program. */
#define ERRMSG(msg) {							\
    printf("\nError (%s, %s, l%d): %s\n\n",				\
	   __FILE__, __func__, __LINE__, msg);				\
    exit(EXIT_FAILURE);							\
  }

/*! Compute floating point modulo. */
#define FMOD(x, y)				\
  ((x) - (int) ((x) / (y)) * (y))

/*! Read binary data. */
#define FREAD(ptr, type, size, out) {					\
    if(fread(ptr, sizeof(type), size, out)!=size)			\
      ERRMSG("Error while reading!");					\
  }

/*! Write binary data. */
#define FWRITE(ptr, type, size, out) {					\
    if(fwrite(ptr, sizeof(type), size, out)!=size)			\
      ERRMSG("Error while writing!");					\
  }

/*! Compute linear interpolation. */
#define LIN(x0, y0, x1, y1, x)			\
  ((y0)+((y1)-(y0))/((x1)-(x0))*((x)-(x0)))

/*! Execute netCDF library command and check result. */
#define NC(cmd) {				     \
    if((cmd)!=NC_NOERR)				     \
      ERRMSG(nc_strerror(cmd));			     \
  }

/*! Compute norm of a vector. */
#define NORM(a) sqrt(DOTP(a, a))

/*! Print macro for debugging. */
#define PRINT(format, var)						\
  printf("Print (%s, %s, l%d): %s= "format"\n",				\
	 __FILE__, __func__, __LINE__, #var, var);

/*! Convert altitude to pressure. */
#define P(z) (P0 * exp(-(z) / H0))

/*! Compute relative humidty. */
#define RH(p, t, h2o) (0.263 * 100. * (p) * MH2O / MA * (h2o)	\
		       / exp(17.67 * ((t) - T0) / ((t) - 29.65)))

/*! Compute square. */
#define SQR(x) ((x)*(x))

/*! Compute potential temperature. */
#define THETA(p, t) ((t) * pow(1000. / (p), 0.286))

/*! Get string tokens. */
#define TOK(line, tok, format, var) {					\
    if(((tok)=strtok((line), " \t"))) {					\
      if(sscanf(tok, format, &(var))!=1) continue;			\
    } else ERRMSG("Error while reading!");				\
  }

/*! Compute virtual temperature. */
#define TVIRT(t, h2o) ((t) * (1.0 + 0.609133 * (h2o) * MH2O / MA))

/*! Print warning message. */
#define WARN(msg) {							\
    printf("\nWarning (%s, %s, l%d): %s\n\n",				\
	   __FILE__, __func__, __LINE__, msg);				\
  }

/*! Convert pressure to altitude. */
#define Z(p) (H0 * log(P0 / (p)))

/* ------------------------------------------------------------
   Timers...
   ------------------------------------------------------------ */

/*! Starts a timer. */
#define START_TIMER(id)	timer(#id, id, 1)

/*! Stops a timer. */
#define STOP_TIMER(id) timer(#id, id, 2)

/*! Prints a timer name and its time. */
#define PRINT_TIMER(id)	timer(#id, id, 3)

/*! Maximum number of timers. */
#define NTIMER 20

/*! Timer for initalization. */
#define TIMER_INIT 1

/*! Timer for file input. */
#define TIMER_INPUT 2

/*! Timer for file output. */
#define TIMER_OUTPUT 3

/*! Timer for advection module. */
#define TIMER_ADVECT 4

/*! Timer for decay module. */
#define TIMER_DECAY 5

/*! Timer for mesoscale diffusion module. */
#define TIMER_DIFFMESO 6

/*! Timer for turbulent diffusion module. */
#define TIMER_DIFFTURB 7

/*! Timer for isosurface module module. */
#define TIMER_ISOSURF 8

/*! Timer for interpolation meteorological data. */
#define TIMER_METEO 9

/*! Timer for position module. */
#define TIMER_POSITION 10

/*! Timer for sedimentation module. */
#define TIMER_SEDI 11

/*! Timer for OH chemistry module. */
#define TIMER_OHCHEM 12

/*! Timer for dry deposition module. */
#define TIMER_DRYDEPO 13

/*! Timer for wet deposition module. */
#define TIMER_WETDEPO 14

/*! Timer for total runtime. */
#define TIMER_TOTAL 15

/* ------------------------------------------------------------
   NVIDIA Tools Extension (NVTX)...
   ------------------------------------------------------------ */

#ifdef USE_NVTX
#include "nvToolsExt.h"

/*! Light blue color code (computation on CPUs). */
#define NVTX_CPU 0xFFADD8E6

/*! Dark blue color code (computation on GPUs). */
#define NVTX_GPU 0xFF00008B

/*! Yellow color code (data transfer from CPUs to GPUs). */
#define NVTX_H2D 0xFFFFFF00

/*! Orange color code (data transfer from GPUs to CPUs). */
#define NVTX_D2H 0xFFFF8800

/*! Light red color code (reading data). */
#define NVTX_READ 0xFFFFCCCB

/*! Dark red color code (writing data). */
#define NVTX_WRITE 0xFF8B0000

/*! Gray color code (other). */
#define NVTX_MISC 0xFF808080

/*! Macro for calling nvtxRangePushEx. */
#define RANGE_PUSH(range_title, range_color) {		\
    nvtxEventAttributes_t eventAttrib = {0};		\
    eventAttrib.version = NVTX_VERSION;			\
    eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;	\
    eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;	\
    eventAttrib.colorType = NVTX_COLOR_ARGB;		\
    eventAttrib.color = range_color;			\
    eventAttrib.message.ascii = range_title;		\
    nvtxRangePushEx(&eventAttrib);			\
  }

/*! Macro for calling nvtxRangePop. */
#define RANGE_POP {				\
    nvtxRangePop();				\
  }
#else

/* Empty definitions of RANGE_PUSH and RANGE_POP... */
#define RANGE_PUSH(range_title, range_color) {}
#define RANGE_POP {}
#endif

/* ------------------------------------------------------------
   Structs...
   ------------------------------------------------------------ */

/*! Control parameters. */
typedef struct {

  /*! Number of quantities. */
  int nq;

  /*! Quantity names. */
  char qnt_name[NQ][LEN];

  /*! Quantity units. */
  char qnt_unit[NQ][LEN];

  /*! Quantity output format. */
  char qnt_format[NQ][LEN];

  /*! Quantity array index for ensemble IDs. */
  int qnt_ens;

  /*! Quantity array index for mass. */
  int qnt_m;

  /*! Quantity array index for particle density. */
  int qnt_rho;

  /*! Quantity array index for particle radius. */
  int qnt_r;

  /*! Quantity array index for surface pressure. */
  int qnt_ps;

  /*! Quantity array index for tropopause pressure. */
  int qnt_pt;

  /*! Quantity array index for geopotential height. */
  int qnt_z;

  /*! Quantity array index for pressure. */
  int qnt_p;

  /*! Quantity array index for temperature. */
  int qnt_t;

  /*! Quantity array index for zonal wind. */
  int qnt_u;

  /*! Quantity array index for meridional wind. */
  int qnt_v;

  /*! Quantity array index for vertical velocity. */
  int qnt_w;

  /*! Quantity array index for water vapor vmr. */
  int qnt_h2o;

  /*! Quantity array index for ozone vmr. */
  int qnt_o3;

  /*! Quantity array index for cloud liquid water content. */
  int qnt_lwc;

  /*! Quantity array index for cloud ice water content. */
  int qnt_iwc;

  /*! Quantity array index for cloud top pressure. */
  int qnt_pc;

  /*! Quantity array index for nitric acid vmr. */
  int qnt_hno3;

  /*! Quantity array index for hydroxyl number concentrations. */
  int qnt_oh;

  /*! Quantity array index for relative humidty. */
  int qnt_rh;

  /*! Quantity array index for potential temperature. */
  int qnt_theta;

  /*! Quantity array index for horizontal wind. */
  int qnt_vh;

  /*! Quantity array index for vertical velocity. */
  int qnt_vz;

  /*! Quantity array index for potential vorticity. */
  int qnt_pv;

  /*! Quantity array index for T_ice. */
  int qnt_tice;

  /*! Quantity array index for T_STS. */
  int qnt_tsts;

  /*! Quantity array index for T_NAT. */
  int qnt_tnat;

  /*! Quantity array index for station flag. */
  int qnt_stat;

  /*! Direction flag (1=forward calculation, -1=backward calculation). */
  int direction;

  /*! Start time of simulation [s]. */
  double t_start;

  /*! Stop time of simulation [s]. */
  double t_stop;

  /*! Time step of simulation [s]. */
  double dt_mod;

  /*! Time step of meteorological data [s]. */
  double dt_met;

  /*! Stride for longitudes. */
  int met_dx;

  /*! Stride for latitudes. */
  int met_dy;

  /*! Stride for pressure levels. */
  int met_dp;

  /*! Smoothing for longitudes. */
  int met_sx;

  /*! Smoothing for latitudes. */
  int met_sy;

  /*! Smoothing for pressure levels. */
  int met_sp;

  /*! Number of target pressure levels. */
  int met_np;

  /*! Target pressure levels [hPa]. */
  double met_p[EP];

  /*! Tropopause definition
     (0=none, 1=clim, 2=cold point, 3=WMO_1st, 4=WMO_2nd). */
  int met_tropo;

  /*! Surface geopotential data file. */
  char met_geopot[LEN];

  /*! Time step for sampling of meteo data along trajectories [s]. */
  double met_dt_out;

  /*! Command to stage meteo data. */
  char met_stage[LEN];

  /*! Isosurface parameter
     (0=none, 1=pressure, 2=density, 3=theta, 4=balloon). */
  int isosurf;

  /*! Balloon position filename. */
  char balloon[LEN];

  /*! Horizontal turbulent diffusion coefficient (troposphere) [m^2/s]. */
  double turb_dx_trop;

  /*! Horizontal turbulent diffusion coefficient (stratosphere) [m^2/s]. */
  double turb_dx_strat;

  /*! Vertical turbulent diffusion coefficient (troposphere) [m^2/s]. */
  double turb_dz_trop;

  /*! Vertical turbulent diffusion coefficient (stratosphere) [m^2/s]. */
  double turb_dz_strat;

  /*! Horizontal scaling factor for mesoscale wind fluctuations. */
  double turb_mesox;

  /*! Vertical scaling factor for mesoscale wind fluctuations. */
  double turb_mesoz;

  /*! Species. */
  char species[LEN];

  /*! Molar mass [g/mol]. */
  double molmass;

  /*! Life time of particles (troposphere) [s]. */
  double tdec_trop;

  /*! Life time of particles (stratosphere)  [s]. */
  double tdec_strat;

  /*! Coefficients for OH chemistry (k0, n, kinf, m). */
  double oh_chem[4];

  /*! Coefficients for dry deposition (v). */
  double dry_depo[1];

  /*! Coefficients for wet deposition (A, B, H). */
  double wet_depo[4];

  /*! H2O volume mixing ratio for PSC analysis. */
  double psc_h2o;

  /*! HNO3 volume mixing ratio for PSC analysis. */
  double psc_hno3;

  /*! Basename of atmospheric data files. */
  char atm_basename[LEN];

  /*! Gnuplot file for atmospheric data. */
  char atm_gpfile[LEN];

  /*! Time step for atmospheric data output [s]. */
  double atm_dt_out;

  /*! Time filter for atmospheric data output (0=no, 1=yes). */
  int atm_filter;

  /*! Particle index stride for atmospheric data files. */
  int atm_stride;

  /*! Type of atmospheric data files (0=ASCII, 1=binary, 2=netCDF). */
  int atm_type;

  /*! Basename of CSI data files. */
  char csi_basename[LEN];

  /*! Time step for CSI data output [s]. */
  double csi_dt_out;

  /*! Observation data file for CSI analysis. */
  char csi_obsfile[LEN];

  /*! Minimum observation index to trigger detection. */
  double csi_obsmin;

  /*! Minimum column density to trigger detection [kg/m^2]. */
  double csi_modmin;

  /*! Number of altitudes of gridded CSI data. */
  int csi_nz;

  /*! Lower altitude of gridded CSI data [km]. */
  double csi_z0;

  /*! Upper altitude of gridded CSI data [km]. */
  double csi_z1;

  /*! Number of longitudes of gridded CSI data. */
  int csi_nx;

  /*! Lower longitude of gridded CSI data [deg]. */
  double csi_lon0;

  /*! Upper longitude of gridded CSI data [deg]. */
  double csi_lon1;

  /*! Number of latitudes of gridded CSI data. */
  int csi_ny;

  /*! Lower latitude of gridded CSI data [deg]. */
  double csi_lat0;

  /*! Upper latitude of gridded CSI data [deg]. */
  double csi_lat1;

  /*! Basename of grid data files. */
  char grid_basename[LEN];

  /*! Gnuplot file for gridded data. */
  char grid_gpfile[LEN];

  /*! Time step for gridded data output [s]. */
  double grid_dt_out;

  /*! Sparse output in grid data files (0=no, 1=yes). */
  int grid_sparse;

  /*! Number of altitudes of gridded data. */
  int grid_nz;

  /*! Lower altitude of gridded data [km]. */
  double grid_z0;

  /*! Upper altitude of gridded data [km]. */
  double grid_z1;

  /*! Number of longitudes of gridded data. */
  int grid_nx;

  /*! Lower longitude of gridded data [deg]. */
  double grid_lon0;

  /*! Upper longitude of gridded data [deg]. */
  double grid_lon1;

  /*! Number of latitudes of gridded data. */
  int grid_ny;

  /*! Lower latitude of gridded data [deg]. */
  double grid_lat0;

  /*! Upper latitude of gridded data [deg]. */
  double grid_lat1;

  /*! Basename for profile output file. */
  char prof_basename[LEN];

  /*! Observation data file for profile output. */
  char prof_obsfile[LEN];

  /*! Number of altitudes of gridded profile data. */
  int prof_nz;

  /*! Lower altitude of gridded profile data [km]. */
  double prof_z0;

  /*! Upper altitude of gridded profile data [km]. */
  double prof_z1;

  /*! Number of longitudes of gridded profile data. */
  int prof_nx;

  /*! Lower longitude of gridded profile data [deg]. */
  double prof_lon0;

  /*! Upper longitude of gridded profile data [deg]. */
  double prof_lon1;

  /*! Number of latitudes of gridded profile data. */
  int prof_ny;

  /*! Lower latitude of gridded profile data [deg]. */
  double prof_lat0;

  /*! Upper latitude of gridded profile data [deg]. */
  double prof_lat1;

  /*! Basename of ensemble data file. */
  char ens_basename[LEN];

  /*! Basename of station data file. */
  char stat_basename[LEN];

  /*! Longitude of station [deg]. */
  double stat_lon;

  /*! Latitude of station [deg]. */
  double stat_lat;

  /*! Search radius around station [km]. */
  double stat_r;

} ctl_t;

/*! Atmospheric data. */
typedef struct {

  /*! Number of air pacels. */
  int np;

  /*! Time [s]. */
  double time[NP];

  /*! Pressure [hPa]. */
  double p[NP];

  /*! Longitude [deg]. */
  double lon[NP];

  /*! Latitude [deg]. */
  double lat[NP];

  /*! Quantity data (for various, user-defined attributes). */
  double q[NQ][NP];

} atm_t;

/*! Cache data. */
typedef struct {

  /*! Zonal wind perturbation [m/s]. */
  float up[NP];

  /*! Meridional wind perturbation [m/s]. */
  float vp[NP];

  /*! Vertical velocity perturbation [hPa/s]. */
  float wp[NP];

  /*! Isosurface variables. */
  double iso_var[NP];

  /*! Isosurface balloon pressure [hPa]. */
  double iso_ps[NP];

  /*! Isosurface balloon time [s]. */
  double iso_ts[NP];

  /*! Isosurface balloon number of data points. */
  int iso_n;

  /*! Cache for reference time of wind standard deviations. */
  double tsig[EX][EY][EP];

  /*! Cache for zonal wind standard deviations. */
  float usig[EX][EY][EP];

  /*! Cache for meridional wind standard deviations. */
  float vsig[EX][EY][EP];

  /*! Cache for vertical velocity standard deviations. */
  float wsig[EX][EY][EP];

} cache_t;

/*! Meteorological data. */
typedef struct {

  /*! Time [s]. */
  double time;

  /*! Number of longitudes. */
  int nx;

  /*! Number of latitudes. */
  int ny;

  /*! Number of pressure levels. */
  int np;

  /*! Longitude [deg]. */
  double lon[EX];

  /*! Latitude [deg]. */
  double lat[EY];

  /*! Pressure [hPa]. */
  double p[EP];

  /*! Surface pressure [hPa]. */
  float ps[EX][EY];

  /*! Geopotential height at the surface [km]. */
  float zs[EX][EY];

  /*! Tropopause pressure [hPa]. */
  float pt[EX][EY];

  /*! Cloud top pressure [hPa]. */
  float pc[EX][EY];

  /*! Total column cloud water [kg/m^2]. */
  float cl[EX][EY];

  /*! Geopotential height at model levels [km]. */
  float z[EX][EY][EP];

  /*! Temperature [K]. */
  float t[EX][EY][EP];

  /*! Zonal wind [m/s]. */
  float u[EX][EY][EP];

  /*! Meridional wind [m/s]. */
  float v[EX][EY][EP];

  /*! Vertical wind [hPa/s]. */
  float w[EX][EY][EP];

  /*! Potential vorticity [PVU]. */
  float pv[EX][EY][EP];

  /*! Water vapor volume mixing ratio [1]. */
  float h2o[EX][EY][EP];

  /*! Ozone volume mixing ratio [1]. */
  float o3[EX][EY][EP];

  /*! Cloud liquid water content [kg/kg]. */
  float lwc[EX][EY][EP];

  /*! Cloud ice water content [kg/kg]. */
  float iwc[EX][EY][EP];

  /*! Pressure on model levels [hPa]. */
  float pl[EX][EY][EP];

} met_t;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Convert Cartesian coordinates to geolocation. */
void cart2geo(
  double *x,
  double *z,
  double *lon,
  double *lat);

/*! Check if x is finite. */
#ifdef _OPENACC
#pragma acc routine (check_finite)
#endif
int check_finite(
  const double x);

/*! Climatology of HNO3 volume mixing ratios. */
#ifdef _OPENACC
#pragma acc routine (clim_hno3)
#endif
double clim_hno3(
  double t,
  double lat,
  double p);

/*! Climatology of OH number concentrations. */
#ifdef _OPENACC
#pragma acc routine (clim_oh)
#endif
double clim_oh(
  double t,
  double lat,
  double p);

/*! Climatology of tropopause pressure. */
#ifdef _OPENACC
#pragma acc routine (clim_tropo)
#endif
double clim_tropo(
  double t,
  double lat);

/*! Get day of year from date. */
void day2doy(
  int year,
  int mon,
  int day,
  int *doy);

/*! Get date from day of year. */
void doy2day(
  int year,
  int doy,
  int *mon,
  int *day);

/*! Convert geolocation to Cartesian coordinates. */
void geo2cart(
  double z,
  double lon,
  double lat,
  double *x);

/*! Get meteorological data for given timestep. */
void get_met(
  ctl_t * ctl,
  char *metbase,
  double t,
  met_t ** met0,
  met_t ** met1);

/*! Get meteorological data for timestep. */
void get_met_help(
  double t,
  int direct,
  char *metbase,
  double dt_met,
  char *filename);

/*! Replace template strings in filename. */
void get_met_replace(
  char *orig,
  char *search,
  char *repl);

/*! Spatial interpolation of meteorological data. */
#ifdef _OPENACC
#pragma acc routine (intpol_met_space_3d)
#endif
void intpol_met_space_3d(
  met_t * met,
  float array[EX][EY][EP],
  double p,
  double lon,
  double lat,
  double *var,
  int *ci,
  double *cw,
  int init);

/*! Spatial interpolation of meteorological data. */
#ifdef _OPENACC
#pragma acc routine (intpol_met_space_2d)
#endif
void intpol_met_space_2d(
  met_t * met,
  float array[EX][EY],
  double lon,
  double lat,
  double *var,
  int *ci,
  double *cw,
  int init);

/*! Temporal interpolation of meteorological data. */
#ifdef _OPENACC
#pragma acc routine (intpol_met_time_3d)
#endif
void intpol_met_time_3d(
  met_t * met0,
  float array0[EX][EY][EP],
  met_t * met1,
  float array1[EX][EY][EP],
  double ts,
  double p,
  double lon,
  double lat,
  double *var,
  int *ci,
  double *cw,
  int init);

/*! Temporal interpolation of meteorological data. */
#ifdef _OPENACC
#pragma acc routine (intpol_met_time_2d)
#endif
void intpol_met_time_2d(
  met_t * met0,
  float array0[EX][EY],
  met_t * met1,
  float array1[EX][EY],
  double ts,
  double lon,
  double lat,
  double *var,
  int *ci,
  double *cw,
  int init);

/*! Convert seconds to date. */
void jsec2time(
  double jsec,
  int *year,
  int *mon,
  int *day,
  int *hour,
  int *min,
  int *sec,
  double *remain);

/*! Find array index for irregular grid. */
#ifdef _OPENACC
#pragma acc routine (locate_irr)
#endif
int locate_irr(
  double *xx,
  int n,
  double x);

/*! Find array index for regular grid. */
#ifdef _OPENACC
#pragma acc routine (locate_reg)
#endif
int locate_reg(
  double *xx,
  int n,
  double x);

/*! Read atmospheric data. */
int read_atm(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

/*! Read control parameters. */
void read_ctl(
  const char *filename,
  int argc,
  char *argv[],
  ctl_t * ctl);

/*! Read meteorological data file. */
int read_met(
  ctl_t * ctl,
  char *filename,
  met_t * met);

/*! Calculate cloud properties. */
void read_met_cloud(
  met_t * met);

/*! Extrapolate meteorological data at lower boundary. */
void read_met_extrapolate(
  met_t * met);

/*! Calculate geopotential heights. */
void read_met_geopot(
  met_t * met);

/*! Read and convert 3D variable from meteorological data file. */
int read_met_help_3d(
  int ncid,
  char *varname,
  char *varname2,
  met_t * met,
  float dest[EX][EY][EP],
  float scl);

/*! Read and convert 2D variable from meteorological data file. */
int read_met_help_2d(
  int ncid,
  char *varname,
  char *varname2,
  met_t * met,
  float dest[EX][EY],
  float scl);

/*! Convert meteorological data from model levels to pressure levels. */
void read_met_ml2pl(
  ctl_t * ctl,
  met_t * met,
  float var[EX][EY][EP]);

/*! Create meteorological data with periodic boundary conditions. */
void read_met_periodic(
  met_t * met);

/*! Calculate potential vorticity. */
void read_met_pv(
  met_t * met);

/*! Downsampling of meteorological data. */
void read_met_sample(
  ctl_t * ctl,
  met_t * met);

/*! Read surface data. */
void read_met_surface(
  int ncid,
  met_t * met);

/*! Calculate tropopause pressure. */
void read_met_tropo(
  ctl_t * ctl,
  met_t * met);

/*! Read a control parameter from file or command line. */
double scan_ctl(
  const char *filename,
  int argc,
  char *argv[],
  const char *varname,
  int arridx,
  const char *defvalue,
  char *value);

/*! Calculate sedimentation velocity. */
#ifdef _OPENACC
#pragma acc routine (sedi)
#endif
double sedi(
  double p,
  double T,
  double r_p,
  double rho_p);

/*! Spline interpolation. */
void spline(
  double *x,
  double *y,
  int n,
  double *x2,
  double *y2,
  int n2);

/*! Calculate standard deviation. */
#ifdef _OPENACC
#pragma acc routine (stddev)
#endif
double stddev(
  double *data,
  int n);

/*! Convert date to seconds. */
void time2jsec(
  int year,
  int mon,
  int day,
  int hour,
  int min,
  int sec,
  double remain,
  double *jsec);

/*! Measure wall-clock time. */
void timer(
  const char *name,
  int id,
  int mode);

/*! Write atmospheric data. */
void write_atm(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm,
  double t);

/*! Write CSI data. */
void write_csi(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm,
  double t);

/*! Write ensemble data. */
void write_ens(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm,
  double t);

/*! Write gridded data. */
void write_grid(
  const char *filename,
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double t);

/*! Write profile data. */
void write_prof(
  const char *filename,
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double t);

/*! Write station data. */
void write_station(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm,
  double t);
