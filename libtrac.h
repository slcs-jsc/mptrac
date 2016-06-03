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
  
  Copright (C) 2013-2015 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  MPTRAC library declarations.
  
  \mainpage
  
  Massive-Parallel Trajectory Calculations (MPTRAC) is a Lagrangian
  particle dispersion model for the troposphere and stratosphere.
  
  This reference manual provides information on the algorithms
  and data structures used in the code. Further information can be found at:
  http://www.fz-juelich.de/ias/jsc/mptrac
*/

#include <ctype.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include <math.h>
#include <netcdf.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

/*! Allocate and clear memory. */
#define ALLOC(ptr, type, n)                             \
  if((ptr=calloc((size_t)(n), sizeof(type)))==NULL)      \
    ERRMSG("Out of memory!");

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
	   __FILE__, __FUNCTION__, __LINE__, msg);			\
    exit(EXIT_FAILURE);							\
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
	 __FILE__, __FUNCTION__, __LINE__, #var, var);

/*! Convert altitude to pressure. */
#define P(z) (P0*exp(-(z)/H0))

/*! Get string tokens. */
#define TOK(line, tok, format, var) {					\
    if(((tok)=strtok((line), " \t"))) {					\
      if(sscanf(tok, format, &(var))!=1) continue;			\
    } else ERRMSG("Error while reading!");				\
  }

/*! Convert pressure to altitude. */
#define Z(p) (H0*log(P0/(p)))

/*! Starts a timer. */
#define START_TIMER(id)	timer(#id, id, 1)

/*! Stops a timer. */
#define STOP_TIMER(id) timer(#id, id, 2)

/*! Prints a timer name and its time. */
#define PRINT_TIMER(id)	timer(#id, id, 3)

/* ------------------------------------------------------------
   Constants...
   ------------------------------------------------------------ */

/*! Standard gravity [m/s^2]. */
#define G0 9.80665

/*! Scale height [km]. */
#define H0 7.0

/*! Reference pressure [hPa]. */
#define P0 1013.25

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
#define NQ 5

/*! Maximum number of pressure levels for meteorological data. */
#define EP 66

/*! Maximum number of longitudes for meteorological data. */
#define EX 361

/*! Maximum number of latitudes for meteorological data. */
#define EY 181

/*! Maximum number of longitudes for gridded data. */
#define GX 720

/*! Maximum number of latitudes for gridded data. */
#define GY 360

/*! Maximum number of altitudes for gridded data. */
#define GZ 1

/*! Maximum number of OpenMP threads. */
#define NTHREADS 128

/*! Maximum number of timers. */
#define NTIMER 20

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

  /*! Quantity array index for mass. */
  int qnt_m;

  /*! Quantity array index for particle density. */
  int qnt_rho;

  /*! Quantity array index for particle radius. */
  int qnt_r;

  /*! Quantity array index for surface pressure. */
  int qnt_ps;

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

  /*! Quantity array index for potential temperature. */
  int qnt_theta;

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

  /*! Isosurface parameter (0=none, 1=pressure, 2=density, 3=theta). */
  int isosurf;

  /*! Horizontal turbulent diffusion coefficient (troposphere) [m^2/s]. */
  double turb_dx_trop;

  /*! Horizontal turbulent diffusion coefficient (stratosphere) [m^2/s]. */
  double turb_dx_strat;

  /*! Vertical turbulent diffusion coefficient (troposphere) [m^2/s]. */
  double turb_dz_trop;

  /*! Vertical turbulent diffusion coefficient (stratosphere) [m^2/s]. */
  double turb_dz_strat;

  /*! Scaling factor for mesoscale wind fluctuations. */
  double turb_meso;

  /*! Life time of particles (troposphere) [s]. */
  double tdec_trop;

  /*! Life time of particles (stratosphere)  [s]. */
  double tdec_strat;

  /*! Basename of atmospheric data files. */
  char atm_basename[LEN];

  /*! Gnuplot file for atmospheric data. */
  char atm_gpfile[LEN];

  /*! Time step for atmospheric data output [s]. */
  double atm_dt_out;

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

  /*! Basename for sample output file. */
  char sample_basename[LEN];

  /*! Observation data file for sample output. */
  char sample_obsfile[LEN];

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

  /*! Quantitiy data (for various, user-defined attributes). */
  double q[NQ][NP];

  /*! Zonal wind perturbation [m/s]. */
  double up[NP];

  /*! Meridional wind perturbation [m/s]. */
  double vp[NP];

  /*! Vertical velocity perturbation [hPa/s]. */
  double wp[NP];

} atm_t;

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
  double ps[EX][EY];

  /*! Temperature [K]. */
  float t[EX][EY][EP];

  /*! Zonal wind [m/s]. */
  float u[EX][EY][EP];

  /*! Meridional wind [m/s]. */
  float v[EX][EY][EP];

  /*! Vertical wind [hPa/s]. */
  float w[EX][EY][EP];

  /*! Water vapor volume mixing ratio [1]. */
  float h2o[EX][EY][EP];

  /*! Ozone volume mixing ratio [1]. */
  float o3[EX][EY][EP];

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

/*! Convert degrees to horizontal distance. */
double deg2dx(
  double dlon,
  double lat);

/*! Convert degrees to horizontal distance. */
double deg2dy(
  double dlat);

/*! Convert pressure to vertical distance. */
double dp2dz(
  double dp,
  double p);

/*! Convert horizontal distance to degrees. */
double dx2deg(
  double dx,
  double lat);

/*! Convert horizontal distance to degrees. */
double dy2deg(
  double dy);

/*! Convert vertical distance to pressure. */
double dz2dp(
  double dz,
  double p);

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
  met_t * met0,
  met_t * met1);

/*! Get meteorological data for timestep. */
void get_met_help(
  double t,
  int direct,
  char *metbase,
  double dt_met,
  char *filename);

/*! Auxilary function for interpolation of meteorological data. */
void intpol_met_2d(
  double array[EX][EY],
  int ix,
  int iy,
  double wx,
  double wy,
  double *var);

/*! Auxilary function for interpolation of meteorological data. */
void intpol_met_3d(
  float array[EX][EY][EP],
  int ip,
  int ix,
  int iy,
  double wp,
  double wx,
  double wy,
  double *var);

/*! Spatial interpolation of meteorological data. */
void intpol_met_space(
  met_t * met,
  double p,
  double lon,
  double lat,
  double *ps,
  double *t,
  double *u,
  double *v,
  double *w,
  double *h2o,
  double *o3);

/*! Temporal interpolation of meteorological data. */
void intpol_met_time(
  met_t * met0,
  met_t * met1,
  double ts,
  double p,
  double lon,
  double lat,
  double *ps,
  double *t,
  double *u,
  double *v,
  double *w,
  double *h2o,
  double *o3);

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

/*! Find array index. */
int locate(
  double *xx,
  int n,
  double x);

/*! Read atmospheric data. */
void read_atm(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl);

/*! Read control parameters. */
void read_ctl(
  const char *filename,
  int argc,
  char *argv[],
  ctl_t * ctl);

/*! Read meteorological data file. */
void read_met(
  char *filename,
  met_t * met);

/*! Extrapolate meteorological data at lower boundary. */
void read_met_extrapolate(
  met_t * met);

/*! Read and convert variable from meteorological data file. */
void read_met_help(
  int ncid,
  char *varname,
  met_t * met,
  int np,
  float dest[EX][EY][EP],
  float scl);

/*! Create meteorological data with periodic boundary conditions. */
void read_met_periodic(
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

/* Get tropopause pressure... */
double tropopause(
  double t,
  double lat);

/*! Write atmospheric data. */
void write_atm(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl,
  double t);

/*! Write CSI data. */
void write_csi(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl,
  double t);

/*! Write gridded data. */
void write_grid(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl,
  double t);

/*! Write samples of model data. */
void write_sample(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl,
  double t);

/*! Write station data. */
void write_station(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl,
  double t);

/*! Write trajectory data. */
void write_traj(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl,
  double t);
