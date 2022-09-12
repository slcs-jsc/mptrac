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
  
  Copyright (C) 2013-2022 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  MPTRAC library declarations.
*/

/*! 
  \mainpage
  
  Massive-Parallel Trajectory Calculations (MPTRAC) is a Lagrangian
  particle dispersion model for the free troposphere and stratosphere.
  
  This reference manual provides information on the algorithms
  and data structures used in the code.

  Further information can be found at: https://github.com/slcs-jsc/mptrac
*/

#ifndef LIBTRAC_H
#define LIBTRAC_H

/* ------------------------------------------------------------
   Includes...
   ------------------------------------------------------------ */

#include <ctype.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
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

#ifdef ZFP
#include "zfp.h"
#endif

#ifdef ZSTD
#include "zstd.h"
#endif

/* ------------------------------------------------------------
   Constants...
   ------------------------------------------------------------ */

/*! Specific heat of dry air at constant pressure [J/(kg K)]. */
#define CPD 1003.5

/*! Ratio of the specific gas constant of dry air and water vapor [1]. */
#define EPS (MH2O / MA)

/*! Standard gravity [m/s^2]. */
#define G0 9.80665

/*! Scale height [km]. */
#define H0 7.0

/*! Latent heat of vaporization of water [J/kg]. */
#define LV 2501000.

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

/*! Specific gas constant of dry air [J/(kg K)]. */
#define RA (1e3 * RI / MA)

/*! Mean radius of Earth [km]. */
#define RE 6367.421

/*! Ideal gas constant [J/(mol K)]. */
#define RI 8.3144598

/*! Standard temperature [K]. */
#define T0 273.15

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Maximum length of ASCII data lines. */
#ifndef LEN
#define LEN 5000
#endif

/*! Maximum number of atmospheric data points. */
#ifndef NP
#define NP 10000000
#endif

/*! Maximum number of quantities per data point. */
#ifndef NQ
#define NQ 15
#endif

/*! Maximum number of pressure levels for meteo data. */
#ifndef EP
#define EP 140
#endif

/*! Maximum number of longitudes for meteo data. */
#ifndef EX
#define EX 1201
#endif

/*! Maximum number of latitudes for meteo data. */
#ifndef EY
#define EY 601
#endif

/*! Maximum number of longitudes for gridded data. */
#ifndef GX
#define GX 720
#endif

/*! Maximum number of latitudes for gridded data. */
#ifndef GY
#define GY 360
#endif

/*! Maximum number of altitudes for gridded data. */
#ifndef GZ
#define GZ 100
#endif

/*! Maximum number of data points for ensemble analysis. */
#ifndef NENS
#define NENS 2000
#endif

/*! Maximum number of OpenMP threads. */
#ifndef NTHREADS
#define NTHREADS 512
#endif

/*! Maximum number of latitudes for climatological data. */
#ifndef CY
#define CY 250
#endif

/*! Maximum number of pressure levels for climatological data. */
#ifndef CP
#define CP 60
#endif

/*! Maximum number of time steps for climatological data. */
#ifndef CT
#define CT 12
#endif

/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

/*! Allocate and clear memory. */
#ifdef _OPENACC
#define ALLOC(ptr, type, n)				\
  if(acc_get_num_devices(acc_device_nvidia) <= 0)	\
    ERRMSG("Not running on a GPU device!");		\
  if((ptr=calloc((size_t)(n), sizeof(type)))==NULL)	\
    ERRMSG("Out of memory!");
#else
#define ALLOC(ptr, type, n)				 \
  if((ptr=calloc((size_t)(n), sizeof(type)))==NULL)      \
    ERRMSG("Out of memory!");
#endif

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
#define DIST(a, b) \
  sqrt(DIST2(a, b))

/*! Compute squared distance between two vectors. */
#define DIST2(a, b)                                                     \
  ((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))

/*! Compute dot product of two vectors. */
#define DOTP(a, b) \
  (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

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

/*! Initialize cache variables for interpolation. */
#define INTPOL_INIT						\
  double cw[3] = {0.0, 0.0, 0.0}; int ci[3] = {0, 0, 0};

/*! 2-D interpolation of a meteo variable. */
#define INTPOL_2D(var, init)						\
  intpol_met_time_2d(met0, met0->var, met1, met1->var,			\
		     atm->time[ip], atm->lon[ip], atm->lat[ip],		\
		     &var, ci, cw, init);

/*! 3-D interpolation of a meteo variable. */
#define INTPOL_3D(var, init)						\
  intpol_met_time_3d(met0, met0->var, met1, met1->var,			\
		     atm->time[ip], atm->p[ip],				\
		     atm->lon[ip], atm->lat[ip],			\
		     &var, ci, cw, init);

/*! Spatial interpolation of all meteo data. */
#define INTPOL_SPACE_ALL(p, lon, lat) {					\
  intpol_met_space_3d(met, met->z, p, lon, lat, &z, ci, cw, 1);		\
  intpol_met_space_3d(met, met->t, p, lon, lat, &t, ci, cw, 0);		\
  intpol_met_space_3d(met, met->u, p, lon, lat, &u, ci, cw, 0);		\
  intpol_met_space_3d(met, met->v, p, lon, lat, &v, ci, cw, 0);		\
  intpol_met_space_3d(met, met->w, p, lon, lat, &w, ci, cw, 0);		\
  intpol_met_space_3d(met, met->pv, p, lon, lat, &pv, ci, cw, 0);	\
  intpol_met_space_3d(met, met->h2o, p, lon, lat, &h2o, ci, cw, 0);	\
  intpol_met_space_3d(met, met->o3, p, lon, lat, &o3, ci, cw, 0);	\
  intpol_met_space_3d(met, met->lwc, p, lon, lat, &lwc, ci, cw, 0);	\
  intpol_met_space_3d(met, met->iwc, p, lon, lat, &iwc, ci, cw, 0);	\
  intpol_met_space_2d(met, met->ps, lon, lat, &ps, ci, cw, 0);		\
  intpol_met_space_2d(met, met->ts, lon, lat, &ts, ci, cw, 0);		\
  intpol_met_space_2d(met, met->zs, lon, lat, &zs, ci, cw, 0);		\
  intpol_met_space_2d(met, met->us, lon, lat, &us, ci, cw, 0);		\
  intpol_met_space_2d(met, met->vs, lon, lat, &vs, ci, cw, 0);		\
  intpol_met_space_2d(met, met->pbl, lon, lat, &pbl, ci, cw, 0);	\
  intpol_met_space_2d(met, met->pt, lon, lat, &pt, ci, cw, 0);		\
  intpol_met_space_2d(met, met->tt, lon, lat, &tt, ci, cw, 0);		\
  intpol_met_space_2d(met, met->zt, lon, lat, &zt, ci, cw, 0);		\
  intpol_met_space_2d(met, met->h2ot, lon, lat, &h2ot, ci, cw, 0);	\
  intpol_met_space_2d(met, met->pct, lon, lat, &pct, ci, cw, 0);	\
  intpol_met_space_2d(met, met->pcb, lon, lat, &pcb, ci, cw, 0);	\
  intpol_met_space_2d(met, met->cl, lon, lat, &cl, ci, cw, 0);		\
  intpol_met_space_2d(met, met->plcl, lon, lat, &plcl, ci, cw, 0);	\
  intpol_met_space_2d(met, met->plfc, lon, lat, &plfc, ci, cw, 0);	\
  intpol_met_space_2d(met, met->pel, lon, lat, &pel, ci, cw, 0);	\
  intpol_met_space_2d(met, met->cape, lon, lat, &cape, ci, cw, 0);	\
  intpol_met_space_2d(met, met->cin, lon, lat, &cin, ci, cw, 0);	\
  }

/*! Temporal interpolation of all meteo data. */
#define INTPOL_TIME_ALL(time, p, lon, lat) {				\
  intpol_met_time_3d(met0, met0->z, met1, met1->z, time, p, lon, lat, &z, ci, cw, 1); \
  intpol_met_time_3d(met0, met0->t, met1, met1->t, time, p, lon, lat, &t, ci, cw, 0); \
  intpol_met_time_3d(met0, met0->u, met1, met1->u, time, p, lon, lat, &u, ci, cw, 0); \
  intpol_met_time_3d(met0, met0->v, met1, met1->v, time, p, lon, lat, &v, ci, cw, 0); \
  intpol_met_time_3d(met0, met0->w, met1, met1->w, time, p, lon, lat, &w, ci, cw, 0); \
  intpol_met_time_3d(met0, met0->pv, met1, met1->pv, time, p, lon, lat, &pv, ci, cw, 0); \
  intpol_met_time_3d(met0, met0->h2o, met1, met1->h2o, time, p, lon, lat, &h2o, ci, cw, 0); \
  intpol_met_time_3d(met0, met0->o3, met1, met1->o3, time, p, lon, lat, &o3, ci, cw, 0); \
  intpol_met_time_3d(met0, met0->lwc, met1, met1->lwc, time, p, lon, lat, &lwc, ci, cw, 0); \
  intpol_met_time_3d(met0, met0->iwc, met1, met1->iwc, time, p, lon, lat, &iwc, ci, cw, 0); \
  intpol_met_time_2d(met0, met0->ps, met1, met1->ps, time, lon, lat, &ps, ci, cw, 0); \
  intpol_met_time_2d(met0, met0->ts, met1, met1->ts, time, lon, lat, &ts, ci, cw, 0); \
  intpol_met_time_2d(met0, met0->zs, met1, met1->zs, time, lon, lat, &zs, ci, cw, 0); \
  intpol_met_time_2d(met0, met0->us, met1, met1->us, time, lon, lat, &us, ci, cw, 0); \
  intpol_met_time_2d(met0, met0->vs, met1, met1->vs, time, lon, lat, &vs, ci, cw, 0); \
  intpol_met_time_2d(met0, met0->pbl, met1, met1->pbl, time, lon, lat, &pbl, ci, cw, 0); \
  intpol_met_time_2d(met0, met0->pt, met1, met1->pt, time, lon, lat, &pt, ci, cw, 0); \
  intpol_met_time_2d(met0, met0->tt, met1, met1->tt, time, lon, lat, &tt, ci, cw, 0); \
  intpol_met_time_2d(met0, met0->zt, met1, met1->zt, time, lon, lat, &zt, ci, cw, 0); \
  intpol_met_time_2d(met0, met0->h2ot, met1, met1->h2ot, time, lon, lat, &h2ot, ci, cw, 0); \
  intpol_met_time_2d(met0, met0->pct, met1, met1->pct, time, lon, lat, &pct, ci, cw, 0); \
  intpol_met_time_2d(met0, met0->pcb, met1, met1->pcb, time, lon, lat, &pcb, ci, cw, 0); \
  intpol_met_time_2d(met0, met0->cl, met1, met1->cl, time, lon, lat, &cl, ci, cw, 0); \
  intpol_met_time_2d(met0, met0->plcl, met1, met1->plcl, time, lon, lat, &plcl, ci, cw, 0); \
  intpol_met_time_2d(met0, met0->plfc, met1, met1->plfc, time, lon, lat, &plfc, ci, cw, 0); \
  intpol_met_time_2d(met0, met0->pel, met1, met1->pel, time, lon, lat, &pel, ci, cw, 0); \
  intpol_met_time_2d(met0, met0->cape, met1, met1->cape, time, lon, lat, &cape, ci, cw, 0); \
  intpol_met_time_2d(met0, met0->cin, met1, met1->cin, time, lon, lat, &cin, ci, cw, 0); \
  }

/*! Calculate lapse rate between pressure levels. */
#define LAPSE(p1, t1, p2, t2)						\
  (1e3 * G0 / RA * ((t2) - (t1)) / ((t2) + (t1))			\
   * ((p2) + (p1)) / ((p2) - (p1)))

/*! Compute linear interpolation. */
#define LIN(x0, y0, x1, y1, x)			\
  ((y0)+((y1)-(y0))/((x1)-(x0))*((x)-(x0)))

/*! Execute netCDF library command and check result. */
#define NC(cmd) {				     \
  int nc_result=(cmd);				     \
  if(nc_result!=NC_NOERR)			     \
    ERRMSG("%s", nc_strerror(nc_result));	     \
}

/*! Compute nearest neighbor interpolation. */
#define NN(x0, y0, x1, y1, x)				\
  (fabs((x) - (x0)) <= fabs((x) - (x1)) ? (y0) : (y1))

/*! Compute norm of a vector. */
#define NORM(a) \
  sqrt(DOTP(a, a))

/*! Convert altitude to pressure. */
#define P(z)					\
  (P0 * exp(-(z) / H0))

/*! Compute saturation pressure over water (WMO, 2018). */
#define PSAT(t)							\
  (6.112 * exp(17.62 * ((t) - T0) / (243.12 + (t) - T0)))

/*! Compute saturation pressure over ice (Marti and Mauersberger, 1993). */
#define PSICE(t)				\
  (0.01 * pow(10., -2663.5 / (t) + 12.537))

/*! Calculate partial water vapor pressure. */
#define PW(p, h2o)					\
  ((p) * GSL_MAX((h2o), 0.1e-6)				\
   / (1. + (1. - EPS) * GSL_MAX((h2o), 0.1e-6)))

/*! Compute relative humidity over water. */
#define RH(p, t, h2o)				\
  (PW(p, h2o) / PSAT(t) * 100.)

/*! Compute relative humidity over ice. */
#define RHICE(p, t, h2o)			\
  (PW(p, h2o) / PSICE(t) * 100.)

/*! Compute density of air. */
#define RHO(p, t)				\
  (100. * (p) / (RA * (t)))

/*! Set atmospheric quantity value. */
#define SET_ATM(qnt, val)			\
  if (ctl->qnt >= 0)				\
    atm->q[ctl->qnt][ip] = val;

/*! Set atmospheric quantity index. */
#define SET_QNT(qnt, name, unit)			\
  if (strcasecmp(ctl->qnt_name[iq], name) == 0) {	\
    ctl->qnt = iq;					\
    sprintf(ctl->qnt_unit[iq], unit);			\
  } else

/*! Compute specific humidity from water vapor volume mixing ratio. */
#define SH(h2o)					\
  (EPS * GSL_MAX((h2o), 0.1e-6))

/*! Compute square. */
#define SQR(x)					\
  ((x)*(x))

/*! Swap macro. */
#define SWAP(x, y, type)				\
  do {type tmp = x; x = y; y = tmp;} while(0);

/*! Calculate dew point temperature (WMO, 2018). */
#define TDEW(p, h2o)				\
  (T0 + 243.12 * log(PW((p), (h2o)) / 6.112)	\
   / (17.62 - log(PW((p), (h2o)) / 6.112)))

/*! Calculate frost point temperature (Marti and Mauersberger, 1993). */
#define TICE(p, h2o)					\
  (-2663.5 / (log10(100. * PW((p), (h2o))) - 12.537))

/*! Compute potential temperature. */
#define THETA(p, t)				\
  ((t) * pow(1000. / (p), 0.286))

/*! Compute virtual potential temperature. */
#define THETAVIRT(p, t, h2o)				\
  (TVIRT(THETA((p), (t)), GSL_MAX((h2o), 0.1e-6)))

/*! Get string tokens. */
#define TOK(line, tok, format, var) {					\
    if(((tok)=strtok((line), " \t"))) {					\
      if(sscanf(tok, format, &(var))!=1) continue;			\
    } else ERRMSG("Error while reading!");				\
  }

/*! Compute virtual temperature. */
#define TVIRT(t, h2o)					\
  ((t) * (1. + (1. - EPS) * GSL_MAX((h2o), 0.1e-6)))

/*! Convert pressure to altitude. */
#define Z(p)					\
  (H0 * log(P0 / (p)))

/*! Calculate geopotential height difference. */
#define ZDIFF(lnp0, t0, h2o0, lnp1, t1, h2o1)				\
  (RI / MA / G0 * 0.5 * (TVIRT((t0), (h2o0)) + TVIRT((t1), (h2o1)))	\
   * ((lnp0) - (lnp1)))

/*! Calculate zeta vertical coordinate. */
#define ZETA(ps, p, t)							\
  (((p) / (ps) <= 0.3 ? 1. :						\
    sin(M_PI / 2. * (1. - (p) / (ps)) / (1. - 0.3)))			\
   * THETA((p), (t)))

/* ------------------------------------------------------------
   Log messages...
   ------------------------------------------------------------ */

/*! Level of log messages (0=none, 1=basic, 2=detailed, 3=debug). */
#ifndef LOGLEV
#define LOGLEV 2
#endif

/*! Print log message. */
#define LOG(level, ...) {						\
    if(level >= 2)							\
      printf("  ");							\
    if(level <= LOGLEV) {						\
      printf(__VA_ARGS__);						\
      printf("\n");							\
    }									\
  }

/*! Print warning message. */
#define WARN(...) {							\
    printf("\nWarning (%s, %s, l%d): ", __FILE__, __func__, __LINE__);	\
    LOG(0, __VA_ARGS__);						\
  }

/*! Print error message and quit program. */
#define ERRMSG(...) {							\
    printf("\nError (%s, %s, l%d): ", __FILE__, __func__, __LINE__);	\
    LOG(0, __VA_ARGS__);						\
    exit(EXIT_FAILURE);							\
  }

/*! Print macro for debugging. */
#define PRINT(format, var)						\
  printf("Print (%s, %s, l%d): %s= "format"\n",				\
	 __FILE__, __func__, __LINE__, #var, var);

/* ------------------------------------------------------------
   Timers...
   ------------------------------------------------------------ */

/*! Maximum number of timers. */
#define NTIMER 100

/*! Print timers. */
#define PRINT_TIMERS				\
  timer("END", "END", 1);

/*! Select timer. */
#define SELECT_TIMER(id, group, color) {				\
    NVTX_POP;								\
    NVTX_PUSH(id, color);						\
    timer(id, group, 0);						\
  }

/*! Start timers. */
#define START_TIMERS				\
  NVTX_PUSH("START", NVTX_CPU);

/*! Stop timers. */
#define STOP_TIMERS				\
  NVTX_POP;

/* ------------------------------------------------------------
   NVIDIA Tools Extension (NVTX)...
   ------------------------------------------------------------ */

#ifdef NVTX
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

/*! Macro for calling nvtxRangePushEx. */
#define NVTX_PUSH(range_title, range_color) {		\
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
#define NVTX_POP {				\
    nvtxRangePop();				\
  }
#else

/* Empty definitions of NVTX_PUSH and NVTX_POP... */
#define NVTX_PUSH(range_title, range_color) {}
#define NVTX_POP {}
#endif

/* ------------------------------------------------------------
   Thrust...
   ------------------------------------------------------------ */

/*! Wrapper to Thrust sorting function. */
void thrustSortWrapper(
  double *__restrict__ c,
  int n,
  int *__restrict__ index);

/* ------------------------------------------------------------
   Structs...
   ------------------------------------------------------------ */

/*! Control parameters. */
typedef struct {

  ///////////////////////////////////////////////////////////////

  /*! Vertical coordinate of air parcels (0=pressure, 1=zeta) */
  int vert_coord_ap;

  /*! Vertical coordinate of input meteo data (0=automatic, 1=eta). */
  int vert_coord_met;

  /*! Vertical velocity (0=kinematic, 1=diabatic) */
  int vert_vel;

  /*! Read MPTRAC or CLaMS meteo data. (0=MPTRAC, 1=CLaMS) */
  int clams_met_data;

  ///////////////////////////////////////////////////////////////

  /*! Chunk size hint for nc__open. */
  size_t chunkszhint;

  /*! Read mode for nc__open. */
  int read_mode;

  /*! Number of quantities. */
  int nq;

  /*! Quantity names. */
  char qnt_name[NQ][LEN];

  /*! Quantity units. */
  char qnt_unit[NQ][LEN];

  /*! Quantity output format. */
  char qnt_format[NQ][LEN];

  /*! Quantity array index for air parcel IDs. */
  int qnt_idx;

  /*! Quantity array index for ensemble IDs. */
  int qnt_ens;

  /*! Quantity array index for station flag. */
  int qnt_stat;

  /*! Quantity array index for mass. */
  int qnt_m;

  /*! Quantity array index for volume mixing ratio. */
  int qnt_vmr;

  /*! Quantity array index for particle radius. */
  int qnt_rp;

  /*! Quantity array index for particle density. */
  int qnt_rhop;

  /*! Quantity array index for surface pressure. */
  int qnt_ps;

  /*! Quantity array index for surface temperature. */
  int qnt_ts;

  /*! Quantity array index for surface geopotential height. */
  int qnt_zs;

  /*! Quantity array index for surface zonal wind. */
  int qnt_us;

  /*! Quantity array index for surface meridional wind. */
  int qnt_vs;

  /*! Quantity array index for boundary layer pressure. */
  int qnt_pbl;

  /*! Quantity array index for tropopause pressure. */
  int qnt_pt;

  /*! Quantity array index for tropopause temperature. */
  int qnt_tt;

  /*! Quantity array index for tropopause geopotential height. */
  int qnt_zt;

  /*! Quantity array index for tropopause water vapor vmr. */
  int qnt_h2ot;

  /*! Quantity array index for geopotential height. */
  int qnt_z;

  /*! Quantity array index for pressure. */
  int qnt_p;

  /*! Quantity array index for temperature. */
  int qnt_t;

  /*! Quantity array index for density of air. */
  int qnt_rho;

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
  int qnt_pct;

  /*! Quantity array index for cloud bottom pressure. */
  int qnt_pcb;

  /*! Quantity array index for total column cloud water. */
  int qnt_cl;

  /*! Quantity array index for pressure at lifted condensation level (LCL). */
  int qnt_plcl;

  /*! Quantity array index for pressure at level of free convection (LCF). */
  int qnt_plfc;

  /*! Quantity array index for pressure at equilibrium level (EL). */
  int qnt_pel;

  /*! Quantity array index for convective available potential energy (CAPE). */
  int qnt_cape;

  /*! Quantity array index for convective inhibition (CIN). */
  int qnt_cin;

  /*! Quantity array index for nitric acid vmr. */
  int qnt_hno3;

  /*! Quantity array index for hydroxyl number concentrations. */
  int qnt_oh;

  /*! Quantity array index for saturation pressure over water. */
  int qnt_psat;

  /*! Quantity array index for saturation pressure over ice. */
  int qnt_psice;

  /*! Quantity array index for partial water vapor pressure. */
  int qnt_pw;

  /*! Quantity array index for specific humidity. */
  int qnt_sh;

  /*! Quantity array index for relative humidity over water. */
  int qnt_rh;

  /*! Quantity array index for relative humidity over ice. */
  int qnt_rhice;

  /*! Quantity array index for potential temperature. */
  int qnt_theta;

  /*! Quantity array index for zeta vertical coordinate. */
  int qnt_zeta;

  /*! Quantity array index for virtual temperature. */
  int qnt_tvirt;

  /*! Quantity array index for lapse rate. */
  int qnt_lapse;

  /*! Quantity array index for horizontal wind. */
  int qnt_vh;

  /*! Quantity array index for vertical velocity. */
  int qnt_vz;

  /*! Quantity array index for potential vorticity. */
  int qnt_pv;

  /*! Quantity array index for dew point temperature. */
  int qnt_tdew;

  /*! Quantity array index for T_ice. */
  int qnt_tice;

  /*! Quantity array index for T_STS. */
  int qnt_tsts;

  /*! Quantity array index for T_NAT. */
  int qnt_tnat;

  /*! Direction flag (1=forward calculation, -1=backward calculation). */
  int direction;

  /*! Start time of simulation [s]. */
  double t_start;

  /*! Stop time of simulation [s]. */
  double t_stop;

  /*! Time step of simulation [s]. */
  double dt_mod;

  /*! Basename for meteo data. */
  char metbase[LEN];

  /*! Time step of meteo data [s]. */
  double dt_met;

  /*! Type of meteo data files (0=netCDF, 1=binary, 2=pack, 3=zfp, 4=zstd). */
  int met_type;

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

  /*! FWHM of horizontal Gaussian used for detrending [km]. */
  double met_detrend;

  /*! Number of target pressure levels. */
  int met_np;

  /*! Target pressure levels [hPa]. */
  double met_p[EP];

  /*! Longitudinal smoothing of geopotential heights. */
  int met_geopot_sx;

  /*! Latitudinal smoothing of geopotential heights. */
  int met_geopot_sy;

  /*! Tropopause definition
     (0=none, 1=clim, 2=cold point, 3=WMO_1st, 4=WMO_2nd, 5=dynamical). */
  int met_tropo;

  /*! WMO tropopause lapse rate [K/km]. */
  double met_tropo_lapse;

  /*! WMO tropopause layer width (number of levels). */
  int met_tropo_nlev;

  /*! WMO tropopause separation layer lapse rate [K/km]. */
  double met_tropo_lapse_sep;

  /*! WMO tropopause separation layer width (number of levels). */
  int met_tropo_nlev_sep;

  /*! Dyanmical tropopause potential vorticity threshold [PVU]. */
  double met_tropo_pv;

  /*! Dynamical tropopause potential temperature threshold [K]. */
  double met_tropo_theta;

  /*! Tropopause interpolation method (0=linear, 1=spline). */
  int met_tropo_spline;

  /*! Cloud data (0=none, 1=LWC+IWC, 2=RWC+SWC, 3=all). */
  double met_cloud;

  /*! Time step for sampling of meteo data along trajectories [s]. */
  double met_dt_out;

  /*! Preload meteo data into disk cache (0=no, 1=yes). */
  int met_cache;

  /*! Time step for sorting of particle data [s]. */
  double sort_dt;

  /*! Isosurface parameter
     (0=none, 1=pressure, 2=density, 3=theta, 4=balloon). */
  int isosurf;

  /*! Balloon position filename. */
  char balloon[LEN];

  /*! Advection scheme (0=midpoint, 1=Runge-Kutta). */
  int advect;

  /*! Reflection of particles at top and bottom boundary (0=no, 1=yes). */
  int reflect;

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

  /*! CAPE threshold for convection module [J/kg]. */
  double conv_cape;

  /*! CIN threshold for convection module [J/kg]. */
  double conv_cin;

  /*! Maximum vertical velocity for convection module [m/s]. */
  double conv_wmax;

  /*! Limit vertical velocity based on CAPE (0=no, 1=yes). */
  double conv_wcape;

  /*! Time interval for convection module [s]. */
  double conv_dt;

  /*! Lower level for mixing (0=particle pressure, 1=surface). */
  int conv_mix_bot;

  /*! Upper level for mixing (0=particle pressure, 1=EL). */
  int conv_mix_top;

  /*! Boundary conditions mass per particle [kg]. */
  double bound_mass;

  /*! Boundary conditions mass per particle trend [kg/s]. */
  double bound_mass_trend;

  /*! Boundary conditions volume mixing ratio [ppv]. */
  double bound_vmr;

  /*! Boundary conditions volume mixing ratio trend [ppv/s]. */
  double bound_vmr_trend;

  /*! Boundary conditions minimum longitude [deg]. */
  double bound_lat0;

  /*! Boundary conditions maximum longitude [deg]. */
  double bound_lat1;

  /*! Boundary conditions bottom pressure [hPa]. */
  double bound_p0;

  /*! Boundary conditions top pressure [hPa]. */
  double bound_p1;

  /*! Boundary conditions delta to surface pressure [hPa]. */
  double bound_dps;

  /*! Species. */
  char species[LEN];

  /*! Molar mass [g/mol]. */
  double molmass;

  /*! Life time of particles (troposphere) [s]. */
  double tdec_trop;

  /*! Life time of particles (stratosphere)  [s]. */
  double tdec_strat;

  /*! Filename of OH climatology. */
  char clim_oh_filename[LEN];

  /*! Reaction type for OH chemistry (0=none, 2=bimolecular, 3=termolecular). */
  int oh_chem_reaction;

  /*! Coefficients for OH reaction rate (A, E/R or k0, n, kinf, m). */
  double oh_chem[4];

  /*! Beta parameter for diurnal variablity of OH. */
  double oh_chem_beta;

  /*! Coefficients for dry deposition (v). */
  double dry_depo[1];

  /*! Coefficient A for wet deposition below cloud (exponential form). */
  double wet_depo_bc_a;

  /*! Coefficient B for wet deposition below cloud (exponential form). */
  double wet_depo_bc_b;

  /*! Coefficient A for wet deposition in cloud (exponential form). */
  double wet_depo_ic_a;

  /*! Coefficient B for wet deposition in cloud (exponential form). */
  double wet_depo_ic_b;

  /*! Coefficients for wet deposition in cloud (Henry's law: Hb, Cb, pH). */
  double wet_depo_ic_h[3];

  /*! Coefficients for wet deposition below cloud (Henry's law: Hb, Cb). */
  double wet_depo_bc_h[2];

  /*! Coefficients for wet deposition in cloud: retention ratio. */
  double wet_depo_ic_ret_ratio;

  /*! Coefficients for wet deposition below cloud: retention ratio. */
  double wet_depo_bc_ret_ratio;

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

  /*! Time filter for atmospheric data output (0=none, 1=missval, 2=remove). */
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

  /*! Basename of sample data file. */
  char sample_basename[LEN];

  /*! Observation data file for sample output. */
  char sample_obsfile[LEN];

  /*! Horizontal radius for sample output [km]. */
  double sample_dx;

  /*! Layer width for sample output [km]. */
  double sample_dz;

  /*! Basename of station data file. */
  char stat_basename[LEN];

  /*! Longitude of station [deg]. */
  double stat_lon;

  /*! Latitude of station [deg]. */
  double stat_lat;

  /*! Search radius around station [km]. */
  double stat_r;

  /*! Start time for station output [s]. */
  double stat_t0;

  /*! Stop time for station output [s]. */
  double stat_t1;

} ctl_t;

/*! Atmospheric data. */
typedef struct {

  /*! Number of air parcels. */
  int np;

  /*! Time [s]. */
  double time[NP];

  /*! Pressure [hPa]. */
  double p[NP];

  /*! Zeta [K]. */
  double zeta[NP];

  /*! Longitude [deg]. */
  double lon[NP];

  /*! Latitude [deg]. */
  double lat[NP];

  /*! Quantity data (for various, user-defined attributes). */
  double q[NQ][NP];

} atm_t;

/*! Cache data. */
typedef struct {

  /*! Isosurface variables. */
  double iso_var[NP];

  /*! Isosurface balloon pressure [hPa]. */
  double iso_ps[NP];

  /*! Isosurface balloon time [s]. */
  double iso_ts[NP];

  /*! Isosurface balloon number of data points. */
  int iso_n;

  /*! Wind perturbations [m/s]. */
  float uvwp[NP][3];

} cache_t;

/*! Climatological data. */
typedef struct {

  /*! Number of OH data timesteps. */
  int oh_nt;

  /*! Number of OH data latitudes. */
  int oh_ny;

  /*! Number of OH data pressure levels. */
  int oh_np;

  /*! OH data time steps [s]. */
  double oh_time[CT];

  /*! OH data latitudes [deg]. */
  double oh_lat[CY];

  /*! OH data pressure levels [hPa]. */
  double oh_p[CP];

  /*! OH data concentration [molec/cm^3]. */
  double oh[CT][CP][CY];

  /*! Number of H2O2 data timesteps. */
  int h2o2_nt;

  /*! Number of H2O2 data latitudes. */
  int h2o2_ny;

  /*! Number of H2O2 data pressure levels. */
  int h2o2_np;

  /*! H2O2 data time steps [s]. */
  double h2o2_time[CT];

  /*! H2O2 data latitudes [deg]. */
  double h2o2_lat[CY];

  /*! H2O2 data pressure levels [hPa]. */
  double h2o2_p[CP];

  /*! H2O2 data concentration [molec/cm^3]. */
  double h2o2[CT][CP][CY];

} clim_t;

/*! Meteo data. */
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

  /*! Surface temperature [K]. */
  float ts[EX][EY];

  /*! Surface geopotential height [km]. */
  float zs[EX][EY];

  /*! Surface zonal wind [m/s]. */
  float us[EX][EY];

  /*! Surface meridional wind [m/s]. */
  float vs[EX][EY];

  /*! Boundary layer pressure [hPa]. */
  float pbl[EX][EY];

  /*! Tropopause pressure [hPa]. */
  float pt[EX][EY];

  /*! Tropopause temperature [K]. */
  float tt[EX][EY];

  /*! Tropopause geopotential height [km]. */
  float zt[EX][EY];

  /*! Tropopause water vapor vmr [ppv]. */
  float h2ot[EX][EY];

  /*! Cloud top pressure [hPa]. */
  float pct[EX][EY];

  /*! Cloud bottom pressure [hPa]. */
  float pcb[EX][EY];

  /*! Total column cloud water [kg/m^2]. */
  float cl[EX][EY];

  /*! Pressure at lifted condensation level (LCL) [hPa]. */
  float plcl[EX][EY];

  /*! Pressure at level of free convection (LFC) [hPa]. */
  float plfc[EX][EY];

  /*! Pressure at equilibrium level [hPa]. */
  float pel[EX][EY];

  /*! Convective available potential energy [J/kg]. */
  float cape[EX][EY];

  /*! Convective inhibition [J/kg]. */
  float cin[EX][EY];

  /*! Geopotential height at model levels [km]. */
  float z[EX][EY][EP];

  /*! Temperature [K]. */
  float t[EX][EY][EP];

  /*! Zonal wind [m/s]. */
  float u[EX][EY][EP];

  /*! Meridional wind [m/s]. */
  float v[EX][EY][EP];

  /*! Vertical velocity [hPa/s]. */
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

  /*! Cache for wind data. */
  float uvw[EX][EY][EP][3];

} met_t;

/*! HNO3 Volume Climatological. */
typedef struct {
    double secs[12];
    double clim_hno3_lats[18];
    double clim_hno3_ps[10];
    double clim_hno3_var[12][18][10];
}clim_hno3_t;

/*! HNO3 Volume Climatological initial data. */
__attribute__((unused))
static clim_hno3_t clim_hno3_init_data = {
        {
                1209600.00, 3888000.00, 6393600.00,
                9072000.00, 11664000.00, 14342400.00,
                16934400.00, 19612800.00, 22291200.00,
                24883200.00, 27561600.00, 30153600.00
        },
        {
                -85, -75, -65, -55, -45, -35, -25, -15, -5,
                5, 15, 25, 35, 45, 55, 65, 75, 85
        },
        {
                4.64159, 6.81292, 10, 14.678, 21.5443,
                31.6228, 46.4159, 68.1292, 100, 146.78
        },
        {
                {{0.782, 1.65, 2.9, 4.59, 6.71, 8.25, 7.16, 5.75, 2.9, 1.74},
                        {0.529, 1.64, 2.76, 4.55, 6.58, 8, 6.99, 5.55, 2.68, 1.57},
                        {0.723, 1.55, 2.73, 4.48, 6.32, 7.58, 7.05, 5.16, 2.49, 1.54},
                        {0.801, 1.56, 2.74, 4.52, 6.23, 7.35, 6.68, 4.4, 1.97, 1.23},
                        {0.818, 1.62, 2.77, 4.38, 5.98, 6.84, 5.83, 3.05, 1.15, 0.709},
                        {0.901, 1.73, 2.78, 4.21, 5.63, 6.16, 4.68, 1.87, 0.617, 0.37},
                        {0.997, 1.8, 2.79, 4.09, 4.88, 4.96, 3.12, 1.22, 0.311, 0.244},
                        {1, 1.71, 2.51, 3.4, 3.74, 3.39, 2.25, 0.845, 0.204, 0.222},
                        {0.997, 1.7, 2.36, 2.88, 3.01, 2.25, 1.77, 0.608, 0.163, 0.181},
                        {0.991, 1.79, 2.57, 3.06, 3.08, 2.15, 1.81, 0.59, 0.168, 0.104},
                        {0.974, 1.86, 2.84, 3.8, 3.93, 3.79, 2.91, 1.02, 0.152, 0.0985},
                        {0.85, 1.86, 3.3, 5.24, 6.55, 6.86, 5.12, 1.93, 0.378, 0.185},
                        {0.783, 1.89, 3.85, 6.6, 8.56, 8.66, 6.95, 3.95, 1.47, 0.745},
                        {0.883, 2.05, 4.34, 7.54, 9.68, 9.77, 8.19, 5.72, 3.15, 1.77},
                        {1.4, 2.44, 4.72, 8.07, 10.5, 10.9, 9.28, 6.95, 4.47, 2.49},
                        {1.7, 2.43, 4.24, 7.43, 10.4, 11.2, 9.72, 8.15, 5.7, 2.97},
                        {2.06, 2.27, 3.68, 6.77, 10.3, 10.3, 9.05, 9.1, 6.73, 3.14},
                        {2.33, 2.39, 3.51, 6.45, 10.3, 9.88, 8.57, 9.42, 7.22, 3.19}},
                {{0.947, 2.21, 3.81, 5.69, 7.55, 8.63, 7.53, 5.98, 3.03, 1.64},
                        {0.642, 2, 3.4, 5.49, 7.5, 8.52, 7.53, 5.83, 2.74, 1.42},
                        {0.756, 1.83, 3.18, 5.11, 7.24, 8.63, 7.66, 5.5, 2.45, 1.33},
                        {0.837, 1.75, 3.06, 5, 6.79, 8.08, 7.05, 4.42, 1.81, 1.05},
                        {0.86, 1.73, 2.96, 4.68, 6.38, 7.38, 6.09, 2.92, 1.06, 0.661},
                        {0.926, 1.78, 2.89, 4.37, 5.74, 6.14, 4.59, 1.78, 0.561, 0.332},
                        {0.988, 1.78, 2.75, 3.95, 4.64, 4.49, 2.85, 1.13, 0.271, 0.184},
                        {0.999, 1.7, 2.44, 3.27, 3.57, 3.03, 2.06, 0.736, 0.181, 0.189},
                        {0.971, 1.67, 2.23, 2.63, 2.83, 2.15, 1.74, 0.554, 0.157, 0.167},
                        {0.985, 1.72, 2.34, 2.69, 2.81, 2.11, 1.78, 0.592, 0.152, 0.101},
                        {0.95, 1.72, 2.57, 3.44, 3.84, 3.89, 2.91, 0.976, 0.135, 0.114},
                        {0.819, 1.64, 2.93, 4.75, 6.02, 6.93, 5.2, 1.83, 0.347, 0.191},
                        {0.731, 1.58, 3.3, 5.95, 7.81, 8.32, 6.93, 3.83, 1.47, 0.875},
                        {0.77, 1.75, 3.74, 6.67, 8.76, 9.41, 8.19, 5.78, 3.32, 2.11},
                        {1.08, 2.17, 4.24, 7.13, 9.2, 10.3, 9.03, 6.87, 4.65, 3.01},
                        {1.43, 2.49, 4.31, 7, 9.14, 10.6, 9.34, 7.6, 5.86, 3.64},
                        {1.5, 2.68, 4.32, 6.75, 8.78, 10.6, 9.05, 7.65, 6.27, 4.07},
                        {1.73, 2.91, 4.33, 6.67, 8.73, 10.6, 8.5, 7.54, 6.63, 4.17}},
                {{1.43, 3.07, 5.22, 7.54, 9.78, 10.4, 10.1, 7.26, 3.61, 1.69},
                        {0.989, 2.69, 4.76, 7.19, 9.44, 9.94, 9.5, 6.74, 3.24, 1.52},
                        {0.908, 2.23, 4.11, 6.48, 8.74, 9.41, 8.58, 5.8, 2.66, 1.3},
                        {0.923, 1.99, 3.61, 5.83, 7.84, 8.6, 7.55, 4.57, 1.87, 0.98},
                        {0.933, 1.9, 3.31, 5.28, 7.1, 7.84, 6.44, 3.18, 1.1, 0.642},
                        {0.982, 1.88, 3.1, 4.76, 6.16, 6.57, 5.16, 2.04, 0.598, 0.33},
                        {1.02, 1.82, 2.88, 4.12, 4.71, 4.54, 3.03, 1.22, 0.268, 0.174},
                        {0.992, 1.7, 2.51, 3.33, 3.62, 2.87, 2.05, 0.705, 0.161, 0.169},
                        {0.969, 1.69, 2.2, 2.62, 2.84, 2.13, 1.78, 0.529, 0.146, 0.186},
                        {0.945, 1.69, 2.27, 2.64, 2.83, 2.2, 1.83, 0.561, 0.139, 0.121},
                        {0.922, 1.65, 2.48, 3.33, 3.83, 4.09, 2.92, 0.973, 0.117, 0.135},
                        {0.886, 1.59, 2.66, 4.26, 5.51, 6.57, 5.09, 1.79, 0.342, 0.194},
                        {0.786, 1.5, 2.78, 5.01, 6.8, 7.83, 6.65, 3.62, 1.45, 1},
                        {0.745, 1.55, 3.05, 5.49, 7.44, 8.6, 7.8, 5.28, 2.95, 2.12},
                        {0.938, 1.76, 3.4, 5.82, 7.8, 9.04, 8.43, 6.15, 3.85, 2.82},
                        {0.999, 2, 3.66, 5.95, 7.94, 9.27, 8.8, 6.93, 4.87, 3.54},
                        {1.13, 2.23, 3.86, 5.82, 7.65, 9, 8.82, 7.17, 5.72, 4.08},
                        {1.23, 2.33, 3.94, 5.74, 7.48, 8.9, 8.84, 7.35, 6.3, 4.42}},
                {{1.55, 3.2, 6.25, 10, 12.9, 12.9, 11.9, 7.96, 3.96, 1.75},
                        {1.32, 3.27, 6.32, 9.99, 12.7, 12.4, 11.3, 7.51, 3.66, 1.58},
                        {1.25, 3.08, 5.77, 8.71, 11.2, 11.2, 9.84, 6.52, 3.23, 1.5},
                        {1.18, 2.59, 4.76, 7.46, 9.61, 9.66, 8.42, 5.06, 2.25, 1.09},
                        {1.09, 2.24, 3.99, 6.4, 8.33, 8.54, 7.08, 3.69, 1.36, 0.727},
                        {1.06, 2.07, 3.52, 5.52, 7.06, 7.26, 5.83, 2.46, 0.732, 0.409},
                        {1.07, 1.91, 3.09, 4.63, 5.21, 4.9, 3.68, 1.43, 0.326, 0.198},
                        {1.03, 1.74, 2.63, 3.54, 3.78, 2.89, 2.09, 0.743, 0.175, 0.12},
                        {0.959, 1.71, 2.32, 2.77, 2.99, 2.24, 1.76, 0.519, 0.149, 0.172},
                        {0.931, 1.68, 2.32, 2.74, 2.99, 2.46, 1.88, 0.578, 0.156, 0.157},
                        {0.933, 1.66, 2.49, 3.42, 3.99, 4.12, 2.93, 1.02, 0.181, 0.138},
                        {0.952, 1.64, 2.6, 4, 5.15, 6.07, 4.84, 1.78, 0.407, 0.286},
                        {0.84, 1.54, 2.68, 4.47, 5.97, 7.13, 6.23, 3.25, 1.38, 1.02},
                        {0.714, 1.44, 2.73, 4.68, 6.28, 7.68, 7.21, 4.82, 2.55, 1.96},
                        {0.838, 1.57, 2.96, 4.93, 6.55, 8.08, 7.74, 5.77, 3.32, 2.52},
                        {0.823, 1.65, 3.11, 5.09, 6.89, 8.36, 8.31, 6.59, 4.1, 3.04},
                        {0.886, 1.83, 3.42, 5.33, 6.92, 8.36, 8.63, 7.21, 4.82, 3.46},
                        {1.07, 2.12, 3.74, 5.54, 6.98, 8.41, 8.75, 7.41, 5.16, 3.62}},
                {{1.13, 2.59, 7.49, 13.5, 15.4, 12.9, 11.3, 8.62, 4.18, 1.63},
                        {0.973, 2.79, 7.23, 12.8, 15.2, 13.3, 11.6, 8.42, 4.06, 1.57},
                        {1.46, 3.44, 6.78, 10.4, 12.7, 12.1, 10.5, 7.04, 3.59, 1.63},
                        {1.52, 3.38, 6.04, 9.08, 11, 10.3, 8.9, 5.7, 2.77, 1.37},
                        {1.32, 2.65, 4.75, 7.49, 9.32, 8.89, 7.42, 4.27, 1.7, 0.88},
                        {1.19, 2.2, 3.88, 6.36, 8.03, 7.81, 6.19, 2.94, 0.948, 0.527},
                        {1.14, 1.96, 3.28, 5.26, 6.12, 5.8, 4.47, 1.66, 0.388, 0.229},
                        {1.07, 1.82, 2.82, 3.92, 4.03, 3.15, 2.31, 0.871, 0.183, 0.0972},
                        {0.978, 1.77, 2.53, 3.04, 3.1, 2.36, 1.76, 0.575, 0.16, 0.126},
                        {0.962, 1.72, 2.49, 3.01, 3.22, 2.72, 2, 0.716, 0.162, 0.183},
                        {0.968, 1.7, 2.6, 3.57, 4.28, 4.35, 3.09, 1.2, 0.262, 0.18},
                        {0.977, 1.68, 2.71, 4.03, 5.17, 6.01, 4.81, 1.81, 0.473, 0.343},
                        {0.819, 1.58, 2.75, 4.37, 5.8, 6.9, 5.96, 2.95, 1.19, 0.964},
                        {0.672, 1.44, 2.69, 4.42, 5.92, 7.26, 6.79, 4.32, 2.22, 1.83},
                        {0.783, 1.42, 2.65, 4.45, 6.04, 7.57, 7.39, 5.4, 2.94, 2.25},
                        {0.757, 1.43, 2.7, 4.54, 6.14, 7.65, 7.51, 5.95, 3.42, 2.39},
                        {0.758, 1.57, 3.04, 4.88, 6.24, 7.85, 7.58, 6.35, 3.81, 2.52},
                        {0.835, 1.72, 3.35, 5.24, 6.5, 8.1, 7.67, 6.51, 4, 2.6}},
                {{1.5, 2.12, 7.64, 10.5, 5.59, 2.14, 2.2, 3.5, 4.71, 3.26},
                        {1.32, 2.14, 7.23, 12, 9.3, 5.3, 5.11, 5.37, 5.12, 3.05},
                        {1.53, 2.92, 6.9, 11.9, 13.5, 11.3, 9.91, 7.18, 4.75, 2.65},
                        {1.66, 3.48, 6.25, 9.53, 11.3, 10.3, 9.01, 5.76, 2.99, 1.67},
                        {1.54, 3.03, 5.21, 8.03, 9.66, 8.98, 7.5, 4.64, 2.11, 1.13},
                        {1.32, 2.39, 4.03, 6.74, 8.52, 8.05, 6.4, 3.48, 1.2, 0.639},
                        {1.17, 2.08, 3.35, 5.52, 6.86, 6.54, 5.08, 1.97, 0.462, 0.217},
                        {1.07, 1.92, 3.01, 4.24, 4.47, 3.77, 2.77, 1.07, 0.213, 0.0694},
                        {0.992, 1.88, 2.76, 3.39, 3.32, 2.52, 1.8, 0.713, 0.192, 0.136},
                        {0.992, 1.8, 2.63, 3.34, 3.46, 2.95, 2.09, 0.9, 0.242, 0.194},
                        {0.987, 1.77, 2.67, 3.64, 4.37, 4.36, 3, 1.27, 0.354, 0.229},
                        {0.979, 1.74, 2.77, 3.99, 5.12, 5.75, 4.53, 1.75, 0.555, 0.302},
                        {0.832, 1.6, 2.78, 4.32, 5.53, 6.67, 5.69, 2.59, 0.982, 0.66},
                        {0.696, 1.41, 2.64, 4.31, 5.65, 7.14, 6.56, 3.8, 1.75, 1.41},
                        {0.788, 1.36, 2.59, 4.3, 5.73, 7.35, 7.04, 4.82, 2.41, 1.8},
                        {0.761, 1.43, 2.61, 4.28, 5.64, 7.37, 7.11, 5.37, 2.68, 1.9},
                        {0.701, 1.44, 2.82, 4.64, 5.76, 7.63, 7.07, 5.74, 2.98, 1.88},
                        {0.763, 1.5, 2.95, 4.97, 6.08, 7.88, 7.12, 5.98, 3.21, 1.91}},
                {{3.58, 2.59, 6.49, 5.84, 1.63, 0.282, 0.647, 0.371, 1.36, 2.33},
                        {3.09, 2.38, 6.37, 7.66, 4.06, 1.23, 1.8, 1.65, 2.32, 2.78},
                        {2.31, 2.84, 5.58, 9.63, 11, 9.02, 8.2, 6.23, 4.17, 3.08},
                        {1.61, 3.16, 5.72, 9.13, 11.4, 10.4, 9.15, 6.18, 3.52, 2.3},
                        {1.32, 2.8, 4.79, 7.44, 9.43, 8.83, 7.41, 4.9, 2.38, 1.38},
                        {1.14, 2.36, 3.94, 6.41, 8.38, 8.17, 6.53, 3.76, 1.31, 0.656},
                        {1.05, 2.1, 3.36, 5.45, 7.07, 6.98, 5.44, 2.22, 0.52, 0.176},
                        {1.02, 2, 3.05, 4.33, 4.74, 4.21, 3.2, 1.26, 0.277, 0.0705},
                        {1.01, 1.96, 2.9, 3.53, 3.46, 2.69, 1.89, 0.859, 0.254, 0.12},
                        {1.01, 1.86, 2.7, 3.46, 3.59, 3.03, 2.14, 1, 0.34, 0.199},
                        {1.02, 1.81, 2.67, 3.68, 4.39, 4.3, 2.93, 1.35, 0.477, 0.25},
                        {0.991, 1.79, 2.82, 4.05, 5.08, 5.5, 4.21, 1.74, 0.605, 0.259},
                        {0.844, 1.73, 2.87, 4.38, 5.49, 6.47, 5.5, 2.44, 0.85, 0.422},
                        {0.729, 1.57, 2.76, 4.43, 5.73, 7.13, 6.43, 3.52, 1.38, 0.913},
                        {0.819, 1.46, 2.69, 4.45, 5.92, 7.47, 7.05, 4.52, 2, 1.4},
                        {0.783, 1.47, 2.71, 4.48, 5.92, 7.46, 7.16, 5.08, 2.35, 1.56},
                        {0.735, 1.51, 2.96, 4.84, 5.92, 7.77, 7.2, 5.54, 2.56, 1.61},
                        {0.8, 1.61, 3.14, 5.2, 6.26, 8.08, 7.27, 5.72, 2.75, 1.62}},
                {{5, 4.43, 5.53, 5.35, 2.33, 0.384, 0.663, 0.164, 0.692, 1.4},
                        {3.62, 3.79, 4.77, 5.94, 4.12, 1.36, 1.3, 0.973, 1.37, 1.73},
                        {2.11, 2.7, 4.12, 7.14, 9.03, 7.74, 7.12, 5.44, 3.73, 2.6},
                        {1.13, 2.32, 4.12, 6.97, 9.86, 9.69, 8.85, 6.22, 3.59, 2.14},
                        {0.957, 2.28, 4.11, 6.47, 8.66, 8.78, 7.33, 4.94, 2.44, 1.38},
                        {0.881, 2.1, 3.65, 5.94, 7.98, 8.29, 6.69, 3.95, 1.36, 0.672},
                        {0.867, 1.96, 3.26, 5.23, 6.94, 7.2, 5.63, 2.41, 0.578, 0.19},
                        {0.953, 1.94, 2.98, 4.23, 4.83, 4.52, 3.38, 1.34, 0.293, 0.181},
                        {1.01, 1.91, 2.77, 3.35, 3.3, 2.62, 1.99, 0.905, 0.245, 0.107},
                        {1.03, 1.81, 2.57, 3.29, 3.43, 2.87, 2.13, 0.988, 0.306, 0.185},
                        {1.02, 1.78, 2.58, 3.59, 4.19, 4, 2.72, 1.29, 0.389, 0.224},
                        {1.01, 1.84, 2.84, 4.06, 4.9, 5.08, 3.71, 1.64, 0.529, 0.232},
                        {0.902, 1.84, 2.98, 4.43, 5.5, 6.28, 5.18, 2.35, 0.734, 0.341},
                        {0.785, 1.68, 2.93, 4.67, 5.95, 7.3, 6.52, 3.48, 1.24, 0.754},
                        {0.847, 1.62, 2.94, 4.86, 6.38, 7.99, 7.5, 4.64, 1.93, 1.23},
                        {0.8, 1.6, 2.94, 4.95, 6.62, 8.16, 7.91, 5.43, 2.43, 1.45},
                        {0.82, 1.76, 3.37, 5.47, 6.82, 8.24, 7.73, 5.79, 2.69, 1.5},
                        {0.988, 2.05, 3.87, 6.01, 7.18, 8.41, 7.7, 5.93, 2.89, 1.55}},
                {{1.52, 2.7, 3.79, 4.95, 3.8, 1.51, 1.11, 0.784, 1.1, 1.56},
                        {1.19, 2.16, 3.34, 4.76, 4.61, 2.93, 2.07, 1.65, 1.63, 1.74},
                        {0.804, 1.65, 2.79, 4.63, 6.64, 6.95, 6.68, 5.11, 3.3, 2.09},
                        {0.86, 1.8, 3.25, 5.3, 7.91, 8.76, 8.28, 6.01, 3.39, 1.83},
                        {0.859, 1.95, 3.54, 5.64, 7.88, 8.55, 7.3, 4.88, 2.3, 1.22},
                        {0.809, 1.88, 3.38, 5.45, 7.47, 8.02, 6.69, 3.98, 1.35, 0.646},
                        {0.822, 1.81, 3.11, 4.9, 6.62, 6.96, 5.63, 2.47, 0.614, 0.169},
                        {0.92, 1.83, 2.8, 3.93, 4.56, 4.4, 3.25, 1.31, 0.295, 0.0587},
                        {0.986, 1.83, 2.6, 3.13, 3.08, 2.53, 1.94, 0.886, 0.244, 0.0815},
                        {0.997, 1.74, 2.5, 3.16, 3.24, 2.67, 2.05, 0.939, 0.281, 0.147},
                        {1.01, 1.75, 2.57, 3.55, 4.1, 3.81, 2.53, 1.21, 0.354, 0.197},
                        {1.04, 1.88, 2.9, 4.16, 4.95, 4.96, 3.48, 1.63, 0.502, 0.163},
                        {0.967, 1.95, 3.17, 4.72, 5.85, 6.5, 5.34, 2.53, 0.748, 0.303},
                        {0.846, 1.83, 3.23, 5.15, 6.62, 7.82, 6.85, 3.79, 1.36, 0.714},
                        {0.91, 1.81, 3.35, 5.55, 7.32, 8.55, 7.88, 5.03, 2.13, 1.1},
                        {0.87, 1.94, 3.6, 5.97, 7.98, 9.14, 8.71, 6.04, 2.73, 1.41},
                        {1.04, 2.36, 4.22, 6.57, 8.5, 9.53, 9.22, 6.71, 3.2, 1.56},
                        {1.36, 2.84, 4.72, 6.94, 8.81, 9.87, 9.59, 7.1, 3.43, 1.65}},
                {{0.704, 1.4, 2.03, 3.08, 4.64, 4.24, 2.55, 1.57, 1.99, 1.91},
                        {0.484, 1.38, 2.08, 3.54, 5.11, 4.98, 3.73, 2.57, 2.29, 1.84},
                        {0.749, 1.57, 2.63, 4.17, 6.15, 6.97, 6.64, 5.11, 3.35, 1.97},
                        {0.864, 1.69, 3.16, 4.87, 7.13, 8.33, 7.87, 5.9, 3.17, 1.56},
                        {0.861, 1.79, 3.28, 5.2, 7.29, 8.32, 7.38, 4.9, 2.23, 1.11},
                        {0.835, 1.79, 3.19, 4.99, 6.72, 7.58, 6.45, 3.68, 1.25, 0.616},
                        {0.847, 1.8, 3.07, 4.66, 6.12, 6.6, 5.21, 2.18, 0.554, 0.21},
                        {0.941, 1.78, 2.68, 3.68, 4.28, 4.18, 2.97, 1.15, 0.238, 0.0968},
                        {0.98, 1.78, 2.48, 2.99, 2.96, 2.35, 1.88, 0.747, 0.207, 0.105},
                        {0.978, 1.74, 2.51, 3.07, 3.12, 2.36, 1.95, 0.777, 0.216, 0.146},
                        {1.01, 1.79, 2.63, 3.53, 3.95, 3.47, 2.38, 1.08, 0.265, 0.178},
                        {1.06, 1.94, 3.02, 4.43, 5.19, 5.01, 3.68, 1.71, 0.429, 0.14},
                        {0.99, 2.02, 3.38, 5.22, 6.56, 6.91, 5.56, 2.75, 0.816, 0.353},
                        {0.923, 2.05, 3.66, 5.98, 7.78, 8.5, 7.23, 4.26, 1.67, 0.802},
                        {1.08, 2.27, 4.17, 6.8, 8.89, 9.55, 8.59, 5.64, 2.58, 1.2},
                        {1.12, 2.5, 4.52, 7.22, 9.76, 10.3, 9.72, 6.79, 3.32, 1.52},
                        {1.2, 2.64, 4.81, 7.64, 10.5, 11.4, 10.6, 7.65, 3.87, 1.73},
                        {1.4, 2.91, 5.01, 7.75, 10.7, 11.6, 11.1, 8.02, 4.04, 1.8}},
                {{0.75, 1.49, 2.39, 3.39, 4.93, 5.94, 5.03, 2.75, 2.27, 1.78},
                        {0.508, 1.52, 2.38, 3.82, 5.34, 6.13, 5.6, 3.31, 2.42, 1.73},
                        {0.715, 1.56, 2.7, 4.39, 6.18, 6.96, 7.1, 5.04, 3.01, 1.75},
                        {0.813, 1.62, 2.94, 4.65, 6.53, 7.65, 7.52, 5.49, 2.75, 1.41},
                        {0.802, 1.68, 2.97, 4.64, 6.37, 7.53, 7.01, 4.56, 1.9, 0.955},
                        {0.816, 1.75, 3.01, 4.59, 6.15, 7.06, 6.15, 3.38, 1.11, 0.61},
                        {0.867, 1.78, 2.92, 4.35, 5.69, 6.05, 4.73, 1.91, 0.519, 0.269},
                        {0.932, 1.7, 2.55, 3.44, 4.03, 3.98, 2.74, 1.08, 0.247, 0.132},
                        {0.937, 1.74, 2.51, 3.09, 3.11, 2.34, 1.84, 0.67, 0.189, 0.121},
                        {0.942, 1.75, 2.63, 3.3, 3.27, 2.21, 1.87, 0.663, 0.171, 0.147},
                        {0.959, 1.8, 2.82, 3.78, 4.03, 3.37, 2.53, 1.04, 0.199, 0.146},
                        {1.01, 1.9, 3.13, 4.76, 5.63, 5.6, 4.31, 1.83, 0.367, 0.172},
                        {0.989, 2.04, 3.64, 6, 7.62, 7.6, 6, 3.35, 1.05, 0.448},
                        {1.02, 2.28, 4.32, 7.19, 9.21, 9.16, 7.64, 4.97, 2.2, 0.948},
                        {1.26, 2.77, 5.2, 8.31, 10.5, 10.4, 9.01, 6.37, 3.46, 1.56},
                        {1.31, 2.76, 5.23, 8.49, 11.2, 11.3, 10.1, 7.27, 3.98, 1.76},
                        {1.26, 2.5, 5.14, 8.85, 12.3, 12.3, 11.2, 8.13, 4.45, 1.97},
                        {1.35, 2.49, 5.26, 9.16, 13, 12.8, 11.8, 8.57, 4.72, 2.05}},
                {{0.759, 1.54, 2.54, 4.22, 6.26, 7.44, 7.14, 4.99, 2.84, 1.89},
                        {0.508, 1.55, 2.5, 4.29, 6.29, 7.29, 7.07, 5.03, 2.77, 1.74},
                        {0.699, 1.56, 2.62, 4.17, 6.08, 7.38, 7.04, 5.17, 2.81, 1.65},
                        {0.778, 1.5, 2.65, 4.35, 6.07, 7.28, 6.84, 4.8, 2.28, 1.28},
                        {0.772, 1.55, 2.71, 4.3, 5.76, 6.91, 6.2, 3.69, 1.45, 0.837},
                        {0.836, 1.67, 2.78, 4.21, 5.56, 6.41, 5.33, 2.47, 0.807, 0.488},
                        {0.937, 1.79, 2.78, 4.12, 5.17, 5.38, 3.89, 1.47, 0.392, 0.256},
                        {0.97, 1.75, 2.52, 3.39, 3.83, 3.63, 2.48, 0.968, 0.212, 0.198},
                        {0.968, 1.74, 2.5, 3.11, 3.2, 2.34, 1.79, 0.629, 0.169, 0.173},
                        {0.98, 1.8, 2.69, 3.42, 3.4, 2.18, 1.81, 0.606, 0.164, 0.138},
                        {0.975, 1.84, 2.96, 4.08, 4.12, 3.5, 2.79, 1.02, 0.145, 0.133},
                        {0.96, 1.94, 3.27, 5.17, 6.26, 6.35, 4.88, 1.91, 0.329, 0.189},
                        {0.954, 2.06, 3.8, 6.53, 8.46, 8.32, 6.53, 3.83, 1.32, 0.6},
                        {1, 2.34, 4.58, 7.71, 9.68, 9.75, 7.96, 5.45, 2.84, 1.39},
                        {1.24, 2.65, 5.14, 8.51, 10.7, 10.6, 8.96, 6.51, 3.83, 1.85},
                        {1.34, 2.44, 4.99, 8.63, 11.6, 11.4, 10.1, 7.84, 4.77, 2.24},
                        {1.33, 2.1, 4.76, 8.78, 12.2, 11.7, 10.8, 8.68, 5.15, 2.35},
                        {1.42, 2.04, 4.68, 8.92, 12.7, 12, 11.2, 8.99, 5.32, 2.33}}
        }
};

/*! Tropopause Pressure Climatological struct. */
typedef struct {
    double secs[12];
    double lats[73];
    double tps[12][73];
} clim_tropo_t;

/*! Tropopause Pressure Climatological initial data. */
__attribute__((unused))
static clim_tropo_t clim_tropo_init_data = {
        {
                1209600.00, 3888000.00, 6393600.00,
                9072000.00, 11664000.00, 14342400.00,
                16934400.00, 19612800.00, 22291200.00,
                24883200.00, 27561600.00, 30153600.00
        },
        { -90, -87.5, -85, -82.5, -80, -77.5, -75, -72.5, -70, -67.5,
          -65, -62.5, -60, -57.5, -55, -52.5, -50, -47.5, -45, -42.5,
          -40, -37.5, -35, -32.5, -30, -27.5, -25, -22.5, -20, -17.5,
          -15, -12.5, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 12.5,
          15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5,
          45, 47.5, 50, 52.5, 55, 57.5, 60, 62.5, 65, 67.5, 70, 72.5,
          75, 77.5, 80, 82.5, 85, 87.5, 90
        },
        { {324.1, 325.6, 325, 324.3, 322.5, 319.7, 314, 307.2, 301.8, 299.6,
           297.1, 292.2, 285.6, 276.1, 264, 248.9, 231.9, 213.5, 194.4,
           175.3, 157, 140.4, 126.7, 116.3, 109.5, 105.4, 103, 101.4, 100.4,
           99.69, 99.19, 98.84, 98.56, 98.39, 98.39, 98.42, 98.44, 98.54,
           98.68, 98.81, 98.89, 98.96, 99.12, 99.65, 101.4, 105.4, 113.5, 128,
           152.1, 184.7, 214, 234.1, 247.3, 255.8, 262.6, 267.7, 271.7, 275,
           277.2, 279, 280.1, 280.4, 280.6, 280.1, 279.3, 278.3, 276.8, 275.8,
           275.3, 275.6, 275.4, 274.1, 273.5},
          {337.3, 338.7, 337.8, 336.4, 333, 328.8, 321.1, 312.6, 306.6, 303.7,
           300.2, 293.8, 285.4, 273.8, 259.6, 242.7, 224.4, 205.2, 186, 167.5,
           150.3, 135, 122.8, 113.9, 108.2, 104.7, 102.5, 101.1, 100.2, 99.42,
           98.88, 98.52, 98.25, 98.09, 98.07, 98.1, 98.12, 98.2, 98.25, 98.27,
           98.26, 98.27, 98.36, 98.79, 100.2, 104.2, 113.7, 131.2, 159.5, 193,
           220.4, 238.1, 250.2, 258.1, 264.7, 269.7, 273.7, 277.3, 280.2, 282.8,
           284.9, 286.5, 288.1, 288.8, 289, 288.5, 287.2, 286.3, 286.1, 287.2,
           287.5, 286.2, 285.8},
          {335, 336, 335.7, 335.1, 332.3, 328.1, 320.6, 311.8, 305.1, 301.9,
           297.6, 290, 280.4, 268.3, 254.6, 239.6, 223.9, 207.9, 192.2, 176.9,
           161.7, 146.4, 132.2, 120.6, 112.3, 107.2, 104.3, 102.4, 101.3,
           100.4, 99.86, 99.47, 99.16, 98.97, 98.94, 98.97, 99, 99.09, 99.2,
           99.31, 99.35, 99.41, 99.51, 99.86, 101.1, 104.9, 114.3, 131, 156.8,
           186.3, 209.3, 224.6, 236.8, 246.3, 254.9, 262.3, 268.8, 274.8,
           279.9, 284.6, 288.6, 291.6, 294.9, 297.5, 299.8, 301.8, 303.1,
           304.3, 304.9, 306, 306.6, 306.2, 306},
          {306.2, 306.7, 305.7, 307.1, 307.3, 306.4, 301.8, 296.2, 292.4,
           290.3, 287.1, 280.9, 273.4, 264.3, 254.1, 242.8, 231, 219, 207.2,
           195.5, 183.3, 169.7, 154.7, 138.7, 124.1, 113.6, 107.8, 104.7,
           102.8, 101.7, 100.9, 100.4, 100, 99.79, 99.7, 99.66, 99.68, 99.79,
           99.94, 100.2, 100.5, 100.9, 101.4, 102.1, 103.4, 107, 115.2, 129.1,
           148.7, 171, 190.8, 205.6, 218.4, 229.4, 239.6, 248.6, 256.5,
           263.7, 270.3, 276.6, 282.6, 288.1, 294.5, 300.4, 306.3, 311.4,
           315.1, 318.3, 320.3, 322.2, 322.8, 321.5, 321.1},
          {266.5, 264.9, 260.8, 261, 262, 263, 261.3, 259.7, 259.2, 259.8,
           260.1, 258.6, 256.7, 253.6, 249.5, 243.9, 237.4, 230, 222.1, 213.9,
           205, 194.4, 180.4, 161.8, 140.7, 122.9, 112.1, 106.7, 104.1, 102.7,
           101.8, 101.4, 101.1, 101, 101, 101, 101.1, 101.2, 101.5, 101.9,
           102.4, 103, 103.8, 104.9, 106.8, 110.1, 115.6, 124, 135.2, 148.9,
           165.2, 181.3, 198, 211.8, 223.5, 233.8, 242.9, 251.5, 259, 266.2,
           273.1, 279.2, 286.2, 292.8, 299.6, 306, 311.1, 315.5, 318.8, 322.6,
           325.3, 325.8, 325.8},
          {220.1, 218.1, 210.8, 207.2, 207.6, 210.5, 211.4, 213.5, 217.3,
           222.4, 227.9, 232.8, 237.4, 240.8, 242.8, 243, 241.5, 238.6, 234.2,
           228.5, 221, 210.7, 195.1, 172.9, 147.8, 127.6, 115.6, 109.9, 107.1,
           105.7, 105, 104.8, 104.8, 104.9, 105, 105.1, 105.3, 105.5, 105.8,
           106.4, 107, 107.6, 108.1, 108.8, 110, 111.8, 114.2, 117.4, 121.6,
           127.9, 137.3, 151.2, 169.5, 189, 205.8, 218.9, 229.1, 237.8, 245,
           251.5, 257.1, 262.3, 268.2, 274, 280.4, 286.7, 292.4, 297.9, 302.9,
           308.5, 312.2, 313.1, 313.3},
          {187.4, 184.5, 173.3, 166.1, 165.4, 167.8, 169.6, 173.6, 179.6,
           187.9, 198.9, 210, 220.5, 229.2, 235.7, 239.9, 241.8, 241.6, 239.6,
           235.8, 229.4, 218.6, 200.9, 175.9, 149.4, 129.4, 118.3, 113.1,
           110.8, 109.7, 109.3, 109.4, 109.7, 110, 110.2, 110.4, 110.5, 110.7,
           111, 111.4, 111.8, 112.1, 112.3, 112.7, 113.2, 113.9, 115, 116.4,
           117.9, 120.4, 124.1, 130.9, 142.2, 159.6, 179.6, 198.5, 212.9,
           224.2, 232.7, 239.1, 243.8, 247.7, 252.4, 257.3, 263.2, 269.5,
           275.4, 281.1, 286.3, 292, 296.3, 298.2, 298.8},
          {166, 166.4, 155.7, 148.3, 147.1, 149, 152.1, 157, 163.6, 172.4,
           185.3, 199.2, 212.6, 224, 233.2, 239.6, 243.3, 244.6, 243.6, 240.3,
           233.9, 222.6, 203.7, 177, 149.5, 129.7, 119, 114, 111.7, 110.7,
           110.3, 110.3, 110.6, 110.9, 111.1, 111.3, 111.5, 111.6, 111.9,
           112.2, 112.5, 112.6, 112.8, 113, 113.4, 114, 115.1, 116.5, 118.3,
           120.9, 124.4, 130.2, 139.4, 154.6, 173.8, 193.1, 208.1, 220.4,
           230.1, 238.2, 244.7, 249.5, 254.5, 259.3, 264.5, 269.4, 273.7,
           278.2, 282.6, 287.4, 290.9, 292.5, 293},
          {171.9, 172.8, 166.2, 162.3, 161.4, 162.5, 165.2, 169.6, 175.3,
           183.1, 193.8, 205.9, 218.3, 229.6, 238.5, 244.3, 246.9, 246.7,
           243.8, 238.4, 230.2, 217.9, 199.6, 174.9, 148.9, 129.8, 119.5,
           114.8, 112.3, 110.9, 110.3, 110.1, 110.2, 110.3, 110.4, 110.5,
           110.6, 110.8, 111, 111.4, 111.8, 112, 112.2, 112.4, 112.9, 113.6,
           114.7, 116.3, 118.4, 121.9, 127.1, 136.1, 149.8, 168.4, 186.9,
           203.3, 217, 229.1, 238.7, 247, 254, 259.3, 264.3, 268.3, 272.5,
           276.6, 280.4, 284.4, 288.4, 293.3, 297.2, 298.7, 299.1},
          {191.6, 192.2, 189, 188.1, 190.2, 193.7, 197.8, 202.9, 208.5,
           215.6, 224.2, 233.1, 241.2, 247.3, 250.8, 251.3, 248.9, 244.2,
           237.3, 228.4, 217.2, 202.9, 184.5, 162.5, 140.7, 124.8, 116.2,
           111.8, 109.4, 107.9, 107, 106.7, 106.6, 106.6, 106.7, 106.7,
           106.8, 107, 107.4, 108, 108.7, 109.3, 109.8, 110.4, 111.2,
           112.4, 114.2, 116.9, 121.1, 127.9, 139.3, 155.2, 173.6, 190.7,
           206.1, 220.1, 232.3, 243, 251.8, 259.2, 265.7, 270.6, 275.3,
           279.3, 283.3, 286.9, 289.7, 292.8, 296.1, 300.5, 303.9, 304.8,
           305.1},
          {241.5, 239.6, 236.8, 237.4, 239.4, 242.3, 244.2, 246.4, 249.2,
           253.6, 258.6, 262.7, 264.8, 264.2, 260.6, 254.1, 245.5, 235.3,
           223.9, 211.7, 198.3, 183.1, 165.6, 147.1, 130.5, 118.7, 111.9,
           108.1, 105.8, 104.3, 103.4, 102.8, 102.5, 102.4, 102.5, 102.5,
           102.5, 102.7, 103.1, 103.8, 104.6, 105.4, 106.1, 107, 108.2,
           109.9, 112.8, 117.5, 126, 140.4, 161, 181.9, 201.2, 216.8, 230.4,
           241.8, 251.4, 259.9, 266.9, 272.8, 277.4, 280.4, 282.9, 284.6,
           286.1, 287.4, 288.3, 289.5, 290.9, 294.2, 296.9, 297.5, 297.6},
          {301.2, 300.3, 296.6, 295.4, 295, 294.3, 291.2, 287.4, 284.9, 284.7,
           284.1, 281.5, 277.1, 270.4, 261.7, 250.6, 237.6, 223.1, 207.9, 192,
           175.8, 158.8, 142.1, 127.6, 116.8, 109.9, 106, 103.6, 102.1, 101.1,
           100.4, 99.96, 99.6, 99.37, 99.32, 99.32, 99.31, 99.46, 99.77, 100.2,
           100.7, 101.3, 101.8, 102.7, 104.1, 106.8, 111.9, 121, 136.7, 160,
           186.9, 209.9, 228.1, 241.2, 251.5, 259.5, 265.7, 270.9, 274.8, 278,
           280.3, 281.8, 283, 283.3, 283.7, 283.8, 283, 282.2, 281.2, 281.4,
           281.7, 281.1, 281.2}
        }
};

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
  clim_t * clim,
  double t,
  double lat,
  double p);

/*! Climatology of OH number concentrations with diurnal variation. */
#ifdef _OPENACC
#pragma acc routine (clim_oh_diurnal)
#endif
double clim_oh_diurnal(
  ctl_t * ctl,
  clim_t * clim,
  double t,
  double p,
  double lon,
  double lat);

/*! Initialization function for OH climatology. */
void clim_oh_init(
  ctl_t * ctl,
  clim_t * clim);

/*! Apply diurnal correction to OH climatology. */
double clim_oh_init_help(
  double beta,
  double time,
  double lat);

/*! Climatology of tropopause pressure. */
#ifdef _OPENACC
#pragma acc routine (clim_tropo)
#endif
double clim_tropo(
  double t,
  double lat,
  clim_tropo_t *clim_tropo_obj);

/*! Pack or unpack array. */
void compress_pack(
  char *varname,
  float *array,
  size_t nxy,
  size_t nz,
  int decompress,
  FILE * inout);

/*! Compress or decompress array with zfp. */
#ifdef ZFP
void compress_zfp(
  char *varname,
  float *array,
  int nx,
  int ny,
  int nz,
  int precision,
  double tolerance,
  int decompress,
  FILE * inout);
#endif

/*! Compress or decompress array with zstd. */
#ifdef ZSTD
void compress_zstd(
  char *varname,
  float *array,
  size_t n,
  int decompress,
  FILE * inout);
#endif

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

/*! Get meteo data for given time step. */
void get_met(
  ctl_t * ctl,
  double t,
  met_t ** met0,
  met_t ** met1);

/*! Get meteo data for time step. */
void get_met_help(
  ctl_t * ctl,
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

/*! Spatial interpolation of meteo data. */
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

/*! Spatial interpolation of meteo data. */
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

/*! Spatial interpolation of meteo data. */
#ifdef _OPENACC
#pragma acc routine (intpol_met_space_uvw)
#endif
void intpol_met_space_uvw(
  met_t * met,
  double p,
  double lon,
  double lat,
  double *u,
  double *v,
  double *w,
  int *ci,
  double *cw,
  int init);

/*! Temporal interpolation of meteo data. */
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

/*! Temporal interpolation of meteo data. */
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

/*! Temporal interpolation of meteo data. */
#ifdef _OPENACC
#pragma acc routine (intpol_met_time_uvw)
#endif
void intpol_met_time_uvw(
  met_t * met0,
  met_t * met1,
  double ts,
  double p,
  double lon,
  double lat,
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

/*! Calculate moist adiabatic lapse rate. */
#ifdef _OPENACC
#pragma acc routine (lapse_rate)
#endif
double lapse_rate(
  double t,
  double h2o);

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

/*! Calculate NAT existence temperature. */
#ifdef _OPENACC
#pragma acc routine (nat_temperature)
#endif
double nat_temperature(
  double p,
  double h2o,
  double hno3);

/*! Parallel quicksort. */
void quicksort(
  double arr[],
  int brr[],
  int low,
  int high);

/*! Partition function for quicksort. */
int quicksort_partition(
  double arr[],
  int brr[],
  int low,
  int high);

/*! Read atmospheric data. */
int read_atm(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

/*! Read climatological data. */
void read_clim(
  ctl_t * ctl,
  clim_t * clim);

/*! Read control parameters. */
void read_ctl(
  const char *filename,
  int argc,
  char *argv[],
  ctl_t * ctl);

/*! Read meteo data file. */
int read_met(
  char *filename,
  ctl_t * ctl,
  met_t * met);

/*! Read 2-D meteo variable. */
void read_met_bin_2d(
  FILE * out,
  met_t * met,
  float var[EX][EY],
  char *varname);

/*! Read 3-D meteo variable. */
void read_met_bin_3d(
  FILE * in,
  ctl_t * ctl,
  met_t * met,
  float var[EX][EY][EP],
  char *varname,
  int precision,
  double tolerance);

/*! Calculate convective available potential energy. */
void read_met_cape(met_t *met, clim_tropo_t *clim_tropo_obj);

/*! Calculate cloud properties. */
void read_met_cloud(
  met_t * met);

/*! Apply detrending method to temperature and winds. */
void read_met_detrend(
  ctl_t * ctl,
  met_t * met);

/*! Extrapolate meteo data at lower boundary. */
void read_met_extrapolate(
  met_t * met);

/*! Calculate geopotential heights. */
void read_met_geopot(
  ctl_t * ctl,
  met_t * met);

/*! Read coordinates of meteo data. */
void read_met_grid(
  char *filename,
  int ncid,
  ctl_t * ctl,
  met_t * met);

/*! Read meteo data on vertical levels. */
void read_met_levels(
  int ncid,
  ctl_t * ctl,
  met_t * met);

/*! Convert meteo data from model levels to pressure levels. */
void read_met_ml2pl(
  ctl_t * ctl,
  met_t * met,
  float var[EX][EY][EP]);

/*! Read and convert 2D variable from meteo data file. */
int read_met_nc_2d(
  int ncid,
  char *varname,
  char *varname2,
  met_t * met,
  float dest[EX][EY],
  float scl,
  int init);

/*! Read and convert 3D variable from meteo data file. */
int read_met_nc_3d(
  int ncid,
  char *varname,
  char *varname2,
  met_t * met,
  float dest[EX][EY][EP],
  float scl,
  int init);

/*! Calculate pressure of the boundary layer. */
void read_met_pbl(
  met_t * met);

/*! Create meteo data with periodic boundary conditions. */
void read_met_periodic(
  met_t * met);

/*! Calculate potential vorticity. */
void read_met_pv(
  met_t * met);

/*! Downsampling of meteo data. */
void read_met_sample(
  ctl_t * ctl,
  met_t * met);

/*! Read surface data. */
void read_met_surface(
  int ncid,
  met_t * met);

/*! Calculate tropopause data. */
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
  double rp,
  double rhop);

/*! Spline interpolation. */
void spline(
  double *x,
  double *y,
  int n,
  double *x2,
  double *y2,
  int n2,
  int method);

/*! Calculate standard deviation. */
#ifdef _OPENACC
#pragma acc routine (stddev)
#endif
float stddev(
  float *data,
  int n);

/*! Calculate solar zenith angle. */
#ifdef _OPENACC
#pragma acc routine (sza)
#endif
double sza(
  double sec,
  double lon,
  double lat);

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
  const char *group,
  int output);

/*! Get weighting factor based on tropopause distance. */
#ifdef _OPENACC
#pragma acc routine (tropo_weight)
#endif
double tropo_weight(double t, double lat, double p, clim_tropo_t *clim_tropo);

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

/*! Read meteo data file. */
int write_met(
  char *filename,
  ctl_t * ctl,
  met_t * met);

/*! Write 2-D meteo variable. */
void write_met_bin_2d(
  FILE * out,
  met_t * met,
  float var[EX][EY],
  char *varname);

/*! Write 3-D meteo variable. */
void write_met_bin_3d(
  FILE * out,
  ctl_t * ctl,
  met_t * met,
  float var[EX][EY][EP],
  char *varname,
  int precision,
  double tolerance);

/*! Write profile data. */
void write_prof(
  const char *filename,
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double t);

/*! Write sample data. */
void write_sample(
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

#endif /* LIBTRAC_H */
