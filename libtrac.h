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

/*! Creates a timer and return its ID. */
#define CREATE_TIMER(name) timer(name, 0, 0)

/*! Starts a timer. */
#define START_TIMER(id)	timer(NULL, id, 1)

/*! Stops a timer. */
#define STOP_TIMER(id) timer(NULL, id, 2)

/*! Prints a timer name and its time. */
#define PRINT_TIMER(id)	timer(NULL, id, 3)

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
#define NP 20000000

/*! Maximum number of quantities per data point. */
#define NQ 7

/*! Maximum number of pressure levels for meteorological data. */
#define EP 91

/*! Maximum number of longitudes for meteorological data. */
#define EX 1440

/*! Maximum number of latitudes for meteorological data. */
#define EY 721

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

/*! Maximum number of data points for Spearman correlation coefficient. */
#define SPEARMAN_NP 1500000

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

  /*! Quantity units. */
  char qnt_format[NQ][LEN];

  /*! Quantity array index for mass. */
  int qnt_mass;

  /*! Quantity array index for particle density. */
  int qnt_rho_p;

  /*! Quantity array index for particle radius. */
  int qnt_r_p;

  /*! Quantity array index for station flag. */
  int qnt_station;

  /*! Quantity array index for temperature. */
  int qnt_temp;

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

  /*! Factor to reduce horizontal resolution of meteorological data. */
  int red_met;

  /*! Horizontal turbulent diffusion coefficient [m^2/s]. */
  double turb_dx;

  /*! Vertical turbulent diffusion coefficient [m^2/s]. */
  double turb_dz;

  /*! Scaling factor for mesoscale wind fluctuations. */
  double turb_meso;

  /*! Half life time of particles [days]. */
  double t12;

  /*! Basename of atmospheric data files. */
  char atm_basename[LEN];

  /*! Time step for atmospheric data output [s]. */
  double atm_dt_out;

  /*! Input data format (0=ASCII, 1=netCDF). */
  int atm_iformat;

  /*! Output data format (0=ASCII, 1=netCDF). */
  int atm_oformat;

  /*! Basename of CSI data files. */
  char csi_basename[LEN];

  /*! Time step for CSI data output [s]. */
  double csi_dt_out;

  /*! Time step for CSI data update [s]. */
  double csi_dt_update;

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

  /*! minimum model column density requried for spearman recognition */
  double spearman_modmin;

  /*! minimum mean observed value for a spearman recognition */
  double spearman_obsmin;

  /*! Basename of station data files. */
  char stat_basename[LEN];

  /*! Time step for station data output [s]. */
  double stat_dt_out;

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

  /*! Vertical velocity perturbation [m/s]. */
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

  /*! Temperature [K]. */
  float t[EX][EY][EP];

  /*! Zonal wind [m/s]. */
  float u[EX][EY][EP];

  /*! Meridional wind [m/s]. */
  float v[EX][EY][EP];

  /*! Vertical wind [hPa/s]. */
  float w[EX][EY][EP];

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

/*! Extrapolate meteorological data at lower boundary. */
void extrapolate_met(
  met_t * met);

/*! Convert geolocation to Cartesian coordinates. */
void geo2cart(
  double z,
  double lon,
  double lat,
  double *x);

/*! Get meteorological data for given timestep. */
void get_met(
  double t,
  int direct,
  char *metbase,
  double dt_met,
  int reduce,
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
void intpol_met_help(
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
  double *t,
  double *u,
  double *v,
  double *w);

/*! Temporal interpolation of meteorological data. */
void intpol_met_time(
  met_t * met0,
  met_t * met1,
  double ts,
  double p,
  double lon,
  double lat,
  double *t,
  double *u,
  double *v,
  double *w);

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
  const char *dirname,
  const char *filename,
  atm_t * atm,
  ctl_t * ctl);

/*! Read atmospheric data from ASCII file. */
void read_atm_from_ascii(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl);

/*! Read atmospheric data from netCDF file. */
void read_atm_from_netcdf(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl);

/*! Read control parameters. */
void read_ctl(
  const char *dirname,
  const char *filename,
  int argc,
  char *argv[],
  ctl_t * ctl);

/*! Read meteorological data file. */
void read_met(
  char *filename,
  met_t * met);

/*! Read and convert variable from meteorological data file. */
void read_met_help(
  int ncid,
  char *varname,
  met_t * met,
  int np,
  float dest[EX][EY][EP],
  float scl);

/*! Reduce spatial resolution of meteorological data. */
void reduce_met(
  met_t * met,
  int dx,
  int dy,
  int dp);

/*! Read a control parameter from file or command line. */
double scan_ctl(
  const char *dirname,
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
int timer(
  const char *name,
  int id,
  int mode);

/*! Write atmospheric data. */
void write_atm(
  const char *dirname,
  const char *filename,
  atm_t * atm,
  ctl_t * ctl);

/*! Write atmospheric data to NetCDF file. */
void write_atm_to_netcdf(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl);

/*! Write atmospheric data to ASCII file. */
void write_atm_to_ascii(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl);

/*! Write CSI data. */
void write_csi(
  const char *dirname,
  const char *filename,
  atm_t * atm,
  ctl_t * ctl,
  double t,
  double dt,
  int write);

/*! Write gridded data. */
void write_grid(
  const char *dirname,
  const char *filename,
  atm_t * atm,
  ctl_t * ctl,
  double t,
  double dt);

/*! Write station data. */
void write_station(
  const char *dirname,
  const char *filename,
  atm_t * atm,
  ctl_t * ctl,
  double t,
  double dt);
