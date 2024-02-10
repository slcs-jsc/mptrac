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
  
  Copyright (C) 2013-2024 Forschungszentrum Juelich GmbH
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
#include <stdint.h>
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

/*! Avogadro constant [1/mol]. */
#ifndef AVO
#define AVO 6.02214076e23
#endif

/*! Specific heat of dry air at constant pressure [J/(kg K)]. */
#ifndef CPD
#define CPD 1003.5
#endif

/*! Ratio of the specific gas constant of dry air and water vapor [1]. */
#ifndef EPS
#define EPS (MH2O / MA)
#endif

/*! Standard gravity [m/s^2]. */
#ifndef G0
#define G0 9.80665
#endif

/*! Scale height [km]. */
#ifndef H0
#define H0 7.0
#endif

/*! Latent heat of vaporization of water [J/kg]. */
#ifndef LV
#define LV 2501000.
#endif

/*! Boltzmann constant [kg m^2/(K s^2)]. */
#ifndef KB
#define KB 1.3806504e-23
#endif

/*! Molar mass of dry air [g/mol]. */
#ifndef MA
#define MA 28.9644
#endif

/*! Molar mass of water vapor [g/mol]. */
#ifndef MH2O
#define MH2O 18.01528
#endif

/*! Molar mass of ozone [g/mol]. */
#ifndef MO3
#define MO3 48.00
#endif

/*! Standard pressure [hPa]. */
#ifndef P0
#define P0 1013.25
#endif

/*! Specific gas constant of dry air [J/(kg K)]. */
#ifndef RA
#define RA (1e3 * RI / MA)
#endif

/*! Mean radius of Earth [km]. */
#ifndef RE
#define RE 6367.421
#endif

/*! Ideal gas constant [J/(mol K)]. */
#ifndef RI
#define RI 8.3144598
#endif

/*! Standard temperature [K]. */
#ifndef T0
#define T0 273.15
#endif

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

/*! Maximum number of data points for CSI calculation. */
#ifndef NCSI
#define NCSI 1000000
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

/*! Maximum number of data points for ensemble analysis. */
#ifndef NENS
#define NENS 2000
#endif

/*! Maximum number of observation data points. */
#ifndef NOBS
#define NOBS 10000000
#endif

/*! Maximum number of OpenMP threads. */
#ifndef NTHREADS
#define NTHREADS 512
#endif

/*! Maximum number of latitudes for climatological data. */
#ifndef CY
#define CY 250
#endif

/*! Maximum number of total column ozone data for climatological data. */
#ifndef CO3
#define CO3 30
#endif

/*! Maximum number of pressure levels for climatological data. */
#ifndef CP
#define CP 70
#endif

/*! Maximum number of solar zenith angles for climatological data. */
#ifndef CSZA
#define CSZA 50
#endif

/*! Maximum number of time steps for climatological data. */
#ifndef CT
#define CT 12
#endif

/*! Maximum number of data points of climatological time series. */
#ifndef CTS
#define CTS 1000
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

/*! Get 2-D array index. */
#define ARRAY_2D(ix, iy, ny)			\
  ((ix) * (ny) + (iy))

/*! Get 3-D array index. */
#define ARRAY_3D(ix, iy, ny, iz, nz)		\
  (((ix)*(ny) + (iy)) * (nz) + (iz))

/*! Arrhenius equation for temperature dependence of reaction rates. */
#define ARRHENIUS(a, b, t)			\
  ((a) * exp( -(b) / (t)))

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
#define DIST(a, b)				\
  sqrt(DIST2(a, b))

/*! Compute squared distance between two vectors. */
#define DIST2(a, b)                                                     \
  ((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))

/*! Compute dot product of two vectors. */
#define DOTP(a, b)				\
  (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

/*! Compute floating point modulo. */
#define FMOD(x, y)				\
  ((x) - (int) ((x) / (y)) * (y))

/*! Read binary data. */
#define FREAD(ptr, type, size, in) {					\
    if(fread(ptr, sizeof(type), size, in)!=size)			\
      ERRMSG("Error while reading!");					\
  }

/*! Write binary data. */
#define FWRITE(ptr, type, size, out) {					\
    if(fwrite(ptr, sizeof(type), size, out)!=size)			\
      ERRMSG("Error while writing!");					\
  }

/*! Initialize cache variables for interpolation. */
#define INTPOL_INIT						\
  double cw[4] = {0.0, 0.0, 0.0, 0.0}; int ci[3] = {0, 0, 0};

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
    intpol_met_space_3d(met, met->z, p, lon, lat, &z, ci, cw, 1);	\
    intpol_met_space_3d(met, met->t, p, lon, lat, &t, ci, cw, 0);	\
    intpol_met_space_3d(met, met->u, p, lon, lat, &u, ci, cw, 0);	\
    intpol_met_space_3d(met, met->v, p, lon, lat, &v, ci, cw, 0);	\
    intpol_met_space_3d(met, met->w, p, lon, lat, &w, ci, cw, 0);	\
    intpol_met_space_3d(met, met->pv, p, lon, lat, &pv, ci, cw, 0);	\
    intpol_met_space_3d(met, met->h2o, p, lon, lat, &h2o, ci, cw, 0);	\
    intpol_met_space_3d(met, met->o3, p, lon, lat, &o3, ci, cw, 0);	\
    intpol_met_space_3d(met, met->lwc, p, lon, lat, &lwc, ci, cw, 0);	\
    intpol_met_space_3d(met, met->iwc, p, lon, lat, &iwc, ci, cw, 0);	\
    intpol_met_space_3d(met, met->cc, p, lon, lat, &cc, ci, cw, 0);	\
    intpol_met_space_2d(met, met->ps, lon, lat, &ps, ci, cw, 0);	\
    intpol_met_space_2d(met, met->ts, lon, lat, &ts, ci, cw, 0);	\
    intpol_met_space_2d(met, met->zs, lon, lat, &zs, ci, cw, 0);	\
    intpol_met_space_2d(met, met->us, lon, lat, &us, ci, cw, 0);	\
    intpol_met_space_2d(met, met->vs, lon, lat, &vs, ci, cw, 0);	\
    intpol_met_space_2d(met, met->lsm, lon, lat, &lsm, ci, cw, 0);	\
    intpol_met_space_2d(met, met->sst, lon, lat, &sst, ci, cw, 0);	\
    intpol_met_space_2d(met, met->pbl, lon, lat, &pbl, ci, cw, 0);	\
    intpol_met_space_2d(met, met->pt, lon, lat, &pt, ci, cw, 0);	\
    intpol_met_space_2d(met, met->tt, lon, lat, &tt, ci, cw, 0);	\
    intpol_met_space_2d(met, met->zt, lon, lat, &zt, ci, cw, 0);	\
    intpol_met_space_2d(met, met->h2ot, lon, lat, &h2ot, ci, cw, 0);	\
    intpol_met_space_2d(met, met->pct, lon, lat, &pct, ci, cw, 0);	\
    intpol_met_space_2d(met, met->pcb, lon, lat, &pcb, ci, cw, 0);	\
    intpol_met_space_2d(met, met->cl, lon, lat, &cl, ci, cw, 0);	\
    intpol_met_space_2d(met, met->plcl, lon, lat, &plcl, ci, cw, 0);	\
    intpol_met_space_2d(met, met->plfc, lon, lat, &plfc, ci, cw, 0);	\
    intpol_met_space_2d(met, met->pel, lon, lat, &pel, ci, cw, 0);	\
    intpol_met_space_2d(met, met->cape, lon, lat, &cape, ci, cw, 0);	\
    intpol_met_space_2d(met, met->cin, lon, lat, &cin, ci, cw, 0);	\
    intpol_met_space_2d(met, met->o3c, lon, lat, &o3c, ci, cw, 0);	\
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
    intpol_met_time_3d(met0, met0->cc, met1, met1->cc, time, p, lon, lat, &cc, ci, cw, 0); \
    intpol_met_time_2d(met0, met0->ps, met1, met1->ps, time, lon, lat, &ps, ci, cw, 0); \
    intpol_met_time_2d(met0, met0->ts, met1, met1->ts, time, lon, lat, &ts, ci, cw, 0); \
    intpol_met_time_2d(met0, met0->zs, met1, met1->zs, time, lon, lat, &zs, ci, cw, 0); \
    intpol_met_time_2d(met0, met0->us, met1, met1->us, time, lon, lat, &us, ci, cw, 0); \
    intpol_met_time_2d(met0, met0->vs, met1, met1->vs, time, lon, lat, &vs, ci, cw, 0); \
    intpol_met_time_2d(met0, met0->lsm, met1, met1->lsm, time, lon, lat, &lsm, ci, cw, 0); \
    intpol_met_time_2d(met0, met0->sst, met1, met1->sst, time, lon, lat, &sst, ci, cw, 0); \
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
    intpol_met_time_2d(met0, met0->o3c, met1, met1->o3c, time, lon, lat, &o3c, ci, cw, 0); \
  }

/*! Calculate lapse rate between pressure levels. */
#define LAPSE(p1, t1, p2, t2)						\
  (1e3 * G0 / RA * ((t2) - (t1)) / ((t2) + (t1))			\
   * ((p2) + (p1)) / ((p2) - (p1)))

/*! Compute linear interpolation. */
#define LIN(x0, y0, x1, y1, x)			\
  ((y0)+((y1)-(y0))/((x1)-(x0))*((x)-(x0)))

/*! Write header for meteo data files. */
#define MET_HEADER							\
  fprintf(out,								\
	  "# $1 = time [s]\n"						\
	  "# $2 = altitude [km]\n"					\
	  "# $3 = longitude [deg]\n"					\
	  "# $4 = latitude [deg]\n"					\
	  "# $5 = pressure [hPa]\n"					\
	  "# $6 = temperature [K]\n"					\
	  "# $7 = zonal wind [m/s]\n"					\
	  "# $8 = meridional wind [m/s]\n"				\
	  "# $9 = vertical velocity [hPa/s]\n"				\
	  "# $10 = H2O volume mixing ratio [ppv]\n");			\
  fprintf(out,								\
	  "# $11 = O3 volume mixing ratio [ppv]\n"			\
	  "# $12 = geopotential height [km]\n"				\
	  "# $13 = potential vorticity [PVU]\n"				\
	  "# $14 = surface pressure [hPa]\n"				\
	  "# $15 = surface temperature [K]\n"				\
	  "# $16 = surface geopotential height [km]\n"			\
	  "# $17 = surface zonal wind [m/s]\n"				\
	  "# $18 = surface meridional wind [m/s]\n"			\
    	  "# $19 = land-sea mask [1]\n"					\
    	  "# $20 = sea surface temperature [K]\n");			\
  fprintf(out,								\
	  "# $21 = tropopause pressure [hPa]\n"				\
	  "# $22 = tropopause geopotential height [km]\n"		\
	  "# $23 = tropopause temperature [K]\n"			\
	  "# $24 = tropopause water vapor [ppv]\n"			\
	  "# $25 = cloud liquid water content [kg/kg]\n"		\
	  "# $26 = cloud ice water content [kg/kg]\n"			\
	  "# $27 = cloud cover [1]\n"					\
	  "# $28 = total column cloud water [kg/m^2]\n"			\
	  "# $29 = cloud top pressure [hPa]\n"				\
	  "# $30 = cloud bottom pressure [hPa]\n");			\
  fprintf(out,								\
	  "# $31 = pressure at lifted condensation level (LCL) [hPa]\n"	\
	  "# $32 = pressure at level of free convection (LFC) [hPa]\n"	\
	  "# $33 = pressure at equilibrium level (EL) [hPa]\n"		\
	  "# $34 = convective available potential energy (CAPE) [J/kg]\n" \
	  "# $35 = convective inhibition (CIN) [J/kg]\n"		\
	  "# $36 = relative humidity over water [%%]\n"			\
	  "# $37 = relative humidity over ice [%%]\n"			\
	  "# $38 = dew point temperature [K]\n"				\
	  "# $39 = frost point temperature [K]\n"			\
	  "# $40 = NAT temperature [K]\n");				\
  fprintf(out,								\
	  "# $41 = HNO3 volume mixing ratio [ppv]\n"			\
	  "# $42 = OH volume mixing ratio [ppv]\n"			\
	  "# $43 = H2O2 volume mixing ratio [ppv]\n"			\
	  "# $44 = HO2 volume mixing ratio [ppv]\n"			\
	  "# $45 = O(1D) volume mixing ratio [ppv]\n"			\
	  "# $46 = boundary layer pressure [hPa]\n"			\
	  "# $47 = total column ozone [DU]\n"				\
	  "# $48 = number of data points\n"				\
	  "# $49 = number of tropopause data points\n"			\
	  "# $50 = number of CAPE data points\n");			\

/*! Calculate molecular density of an ideal gas. */
#define MOLEC_DENS(p,t)			\
  (AVO * 1e-6 * ((p) * 100) / (RI * (t)))

/*! Execute netCDF library command and check result. */
#define NC(cmd) {				     \
    int nc_result=(cmd);			     \
    if(nc_result!=NC_NOERR)			     \
      ERRMSG("%s", nc_strerror(nc_result));	     \
  }

/*! Define netCDF variable. */
#define NC_DEF_VAR(varname, type, ndims, dims, long_name, units) {	\
    NC(nc_def_var(ncid, varname, type, ndims, dims, &varid));		\
    NC(nc_put_att_text(ncid, varid, "long_name", strlen(long_name), long_name)); \
    NC(nc_put_att_text(ncid, varid, "units", strlen(units), units));	\
  }

/*! Read netCDF double array. */
#define NC_GET_DOUBLE(varname, ptr, force) {			\
    if(force) {							\
      NC(nc_inq_varid(ncid, varname, &varid));			\
      NC(nc_get_var_double(ncid, varid, ptr));			\
    } else {							\
      if(nc_inq_varid(ncid, varname, &varid) == NC_NOERR) {	\
	NC(nc_get_var_double(ncid, varid, ptr));		\
      } else							\
	WARN("netCDF variable %s is missing!", varname);	\
    }								\
  }

/*! Read netCDF dimension. */
#define NC_INQ_DIM(dimname, ptr, min, max) {		\
    int dimid; size_t naux;				\
    NC(nc_inq_dimid(ncid, dimname, &dimid));		\
    NC(nc_inq_dimlen(ncid, dimid, &naux));		\
    *ptr = (int)naux;					\
    if ((*ptr) < (min) || (*ptr) > (max))		\
      ERRMSG("Dimension %s is out of range!", dimname);	\
  }

/*! Write netCDF double array. */
#define NC_PUT_DOUBLE(varname, ptr, hyperslab) {		\
    NC(nc_inq_varid(ncid, varname, &varid));			\
    if(hyperslab) {						\
      NC(nc_put_vara_double(ncid, varid, start, count, ptr));	\
    } else {							\
      NC(nc_put_var_double(ncid, varid, ptr));			\
    }								\
  }

/*! Write netCDF integer array. */
#define NC_PUT_INT(varname, ptr, hyperslab) {			\
    NC(nc_inq_varid(ncid, varname, &varid));			\
    if(hyperslab) {						\
      NC(nc_put_vara_int(ncid, varid, start, count, ptr));	\
    } else {							\
      NC(nc_put_var_int(ncid, varid, ptr));			\
    }								\
  }

/*! Set netCDF attribute. */
#define NC_PUT_ATT(varname, attname, text) {				\
    NC(nc_inq_varid(ncid, varname, &varid));				\
    NC(nc_put_att_text(ncid, varid, attname, strlen(text), text));	\
  }

/*! Set netCDF global attribute. */
#define NC_PUT_ATT_GLOBAL(attname, text)				\
  NC(nc_put_att_text(ncid, NC_GLOBAL, attname, strlen(text), text));

/*! Write netCDF float array. */
#define NC_PUT_FLOAT(varname, ptr, hyperslab) {			\
    NC(nc_inq_varid(ncid, varname, &varid));			\
    if(hyperslab) {						\
      NC(nc_put_vara_float(ncid, varid, start, count, ptr));	\
    } else {							\
      NC(nc_put_var_float(ncid, varid, ptr));			\
    }								\
  }

/*! Compute nearest neighbor interpolation. */
#define NN(x0, y0, x1, y1, x)				\
  (fabs((x) - (x0)) <= fabs((x) - (x1)) ? (y0) : (y1))

/*! Compute norm of a vector. */
#define NORM(a)					\
  sqrt(DOTP(a, a))

/*! Convert altitude to pressure. */
#define P(z)					\
  (P0 * exp(-(z) / H0))

/*! Compute saturation pressure over water (WMO, 2018). */
#define PSAT(t)							\
  (6.112 * exp(17.62 * ((t) - T0) / (243.12 + (t) - T0)))

/*! Compute saturation pressure over ice (WMO, 2018). */
#define PSICE(t)						\
  (6.112 * exp(22.46 * ((t) - T0) / (272.62 + (t) - T0)))

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

/*! Roeth approximation formula for photolysis reactions. */
#define ROETH_PHOTOL(a, b, c, sza)					\
  ((c)*(sza) < M_PI/2. ? (a) * exp((b) * (1 - 1/cos((c) * (sza)))) : 0)

/*! Set atmospheric quantity value. */
#define SET_ATM(qnt, val)			\
  if (ctl->qnt >= 0)				\
    atm->q[ctl->qnt][ip] = val;

/*! Set atmospheric quantity index. */
#define SET_QNT(qnt, name, longname, unit)		\
  if (strcasecmp(ctl->qnt_name[iq], name) == 0) {	\
    ctl->qnt = iq;					\
    sprintf(ctl->qnt_longname[iq], longname);		\
    sprintf(ctl->qnt_unit[iq], unit);			\
  } else

/*! Compute specific humidity from water vapor volume mixing ratio. */
#define SH(h2o)					\
  (EPS * GSL_MAX((h2o), 0.1e-6))

/*! Compute square of x. */
#define SQR(x)					\
  ((x)*(x))

/*! Swap macro. */
#define SWAP(x, y, type)			\
  do {type tmp = x; x = y; y = tmp;} while(0);

/*! Calculate dew point temperature (WMO, 2018). */
#define TDEW(p, h2o)				\
  (T0 + 243.12 * log(PW((p), (h2o)) / 6.112)	\
   / (17.62 - log(PW((p), (h2o)) / 6.112)))

/*! Calculate frost point temperature (WMO, 2018). */
#define TICE(p, h2o)				\
  (T0 + 272.62 * log(PW((p), (h2o)) / 6.112)	\
   / (22.46 - log(PW((p), (h2o)) / 6.112)))

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

  /*! Type of meteo data files (0=netCDF, 1=binary, 2=pack, 3=zfp, 4=zstd). */
  int atm_type_out;

  /*! Coupled use of pressure based modules and diabatic advection. 
     (0= no coupling, 1= coupling) */
  int cpl_zeta_and_press_modules;

  /*! Use predefined pressure levels or not. */
  int press_level_def;

  /*! Vertical coordinate of air parcels (0=pressure, 1=zeta). */
  int vert_coord_ap;

  /*! Vertical coordinate of input meteo data (0=automatic, 1=eta). */
  int vert_coord_met;

  /*! Read MPTRAC or CLaMS meteo data (0=MPTRAC, 1=CLaMS). */
  int clams_met_data;

  /*! Number of quantities. */
  int nq;

  /*! Quantity names. */
  char qnt_name[NQ][LEN];

  /*! Quantity long names. */
  char qnt_longname[NQ][LEN];

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

  /*! Quantity array index for land-sea mask. */
  int qnt_lsm;

  /*! Quantity array index for sea surface temperature. */
  int qnt_sst;

  /*! Quantity array index for boundary layer pressure. */
  int qnt_pbl;

  /*! Quantity array index for tropopause pressure. */
  int qnt_pt;

  /*! Quantity array index for tropopause temperature. */
  int qnt_tt;

  /*! Quantity array index for tropopause geopotential height. */
  int qnt_zt;

  /*! Quantity array index for tropopause water vapor volume mixing ratio. */
  int qnt_h2ot;

  /*! Quantity array index for geopotential height. */
  int qnt_zg;

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

  /*! Quantity array index for water vapor volume mixing ratio. */
  int qnt_h2o;

  /*! Quantity array index for ozone volume mixing ratio. */
  int qnt_o3;

  /*! Quantity array index for cloud liquid water content. */
  int qnt_lwc;

  /*! Quantity array index for cloud ice water content. */
  int qnt_iwc;

  /*! Quantity array index for cloud cover. */
  int qnt_cc;

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

  /*! Quantity array index for total column ozone. */
  int qnt_o3c;

  /*! Quantity array index for HNO3 volume mixing ratio (climatology). */
  int qnt_hno3;

  /*! Quantity array index for OH volume mixing ratio (climatology). */
  int qnt_oh;

  /*! Quantity array index for H2O2 volume mixing ratio (climatology). */
  int qnt_h2o2;

  /*! Quantity array index for HO2 volume mixing ratio (climatology). */
  int qnt_ho2;

  /*! Quantity array index for O(1D) volume mixing ratio (climatology). */
  int qnt_o1d;

  /*! Quantity array index for total mass loss due to OH chemistry. */
  int qnt_mloss_oh;

  /*! Quantity array index for total mass loss due to H2O2 chemistry. */
  int qnt_mloss_h2o2;

  /*! Quantity array index for total mass loss due to wet deposition. */
  int qnt_mloss_wet;

  /*! Quantity array index for total mass loss due to dry deposition. */
  int qnt_mloss_dry;

  /*! Quantity array index for total mass loss due to exponential decay. */
  int qnt_mloss_decay;

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

  /*! Quantity array index for diagnosed zeta vertical coordinate. */
  int qnt_zeta_d;

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

  /*! Quantity array index for trace species x volume mixing ratio (chemistry code). */
  int qnt_Cx;

  /*! Quantity array index for H2O volume mixing ratio (chemistry code). */
  int qnt_Ch2o;

  /*! Quantity array index for O3 volume mixing ratio (chemistry code). */
  int qnt_Co3;

  /*! Quantity array index for CO volume mixing ratio (chemistry code). */
  int qnt_Cco;

  /*! Quantity array index for OH volume mixing ratio (chemistry code). */
  int qnt_Coh;

  /*! Quantity array index for H volume mixing ratio (chemistry code). */
  int qnt_Ch;

  /*! Quantity array index for HO2 volume mixing ratio (chemistry code). */
  int qnt_Cho2;

  /*! Quantity array index for H2O2 volume mixing ratio (chemistry code). */
  int qnt_Ch2o2;

  /*! Quantity array index for O(1D) volume mixing ratio (chemistry code). */
  int qnt_Co1d;

  /*! Quantity array index for O(3P) volume mixing ratio (chemistry code). */
  int qnt_Co3p;

  /*! Quantity array index for CFC-10 volume mixing ratio (chemistry code). */
  int qnt_Cccl4;

  /*! Quantity array index for CFC-11 volume mixing ratio (chemistry code). */
  int qnt_Cccl3f;

  /*! Quantity array index for CFC-12 volume mixing ratio (chemistry code). */
  int qnt_Cccl2f2;

  /*! Quantity array index for N2O volume mixing ratio (chemistry code). */
  int qnt_Cn2o;

  /*! Quantity array index for SF6 volume mixing ratio (chemistry code). */
  int qnt_Csf6;

  /*! Quantity array index for age of air. */
  int qnt_aoa;

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

  /*! Meteo data layout (0=[lev, lat, lon], 1 = [lon, lat, lev]). */
  int met_convention;

  /*! Type of meteo data files (0=netCDF, 1=binary, 2=pack, 3=zfp, 4=zstd). */
  int met_type;

  /*! Check netCDF scaling factors (0=no, 1=yes). */
  int met_nc_scale;

  /*! ZFP compression precision (for all variables, except z and T). */
  int met_zfp_prec;

  /*! ZFP compression tolerance (for temperature). */
  double met_zfp_tol_t;

  /*! ZFP compression tolerance (for geopotential height). */
  double met_zfp_tol_z;

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

  /*! Try to read relative humidity (0=no, 1=yes). */
  int met_relhum;

  /*! Tropopause definition
     (0=none, 1=clim, 2=cold point, 3=WMO_1st, 4=WMO_2nd, 5=dynamical). */
  int met_tropo;

  /*! Dynamical tropopause potential vorticity threshold [PVU]. */
  double met_tropo_pv;

  /*! Dynamical tropopause potential temperature threshold [K]. */
  double met_tropo_theta;

  /*! Tropopause interpolation method (0=linear, 1=spline). */
  int met_tropo_spline;

  /*! Cloud data (0=none, 1=LWC+IWC, 2=RWC+SWC, 3=all). */
  int met_cloud;

  /*! Minimum cloud ice water content [kg/kg]. */
  double met_cloud_min;

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

  /*! Advection scheme (0=off, 1=Euler, 2=midpoint, 4=Runge-Kutta). */
  int advect;

  /*! Reflection of particles at top and bottom boundary (0=no, 1=yes). */
  int reflect;

  /*! Random number generator (0=GSL, 1=Squares, 2=cuRAND). */
  int rng_type;

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

  /*! Boundary conditions surface layer depth [hPa]. */
  double bound_dps;

  /*! Boundary conditions surface layer depth [km]. */
  double bound_dzs;

  /*! Boundary conditions surface layer zeta [K]. */
  double bound_zetas;

  /*! Boundary conditions planetary boundary layer (0=no, 1=yes). */
  int bound_pbl;

  /*! Species. */
  char species[LEN];

  /*! Molar mass [g/mol]. */
  double molmass;

  /*! Life time of particles in the troposphere [s]. */
  double tdec_trop;

  /*! Life time of particles in the stratosphere  [s]. */
  double tdec_strat;

  /*! Filename of photolysis rates climatology. */
  char clim_photo[LEN];

  /*! Filename of HNO3 climatology. */
  char clim_hno3_filename[LEN];

  /*! Filename of OH climatology. */
  char clim_oh_filename[LEN];

  /*! Filename of H2O2 climatology. */
  char clim_h2o2_filename[LEN];

  /*! Filename of HO2 climatology. */
  char clim_ho2_filename[LEN];

  /*! Filename of O(1D) climatology. */
  char clim_o1d_filename[LEN];

  /*! Filename of O3 climatology. */
  char clim_o3_filename[LEN];

  /*! Filename of CFC-10 time series. */
  char clim_ccl4_timeseries[LEN];

  /*! Filename of CFC-11 time series. */
  char clim_ccl3f_timeseries[LEN];

  /*! Filename of CFC-12 time series. */
  char clim_ccl2f2_timeseries[LEN];

  /*! Filename of N2O time series. */
  char clim_n2o_timeseries[LEN];

  /*! Filename of SF6 time series. */
  char clim_sf6_timeseries[LEN];

  /*! Time interval for mixing [s]. */
  double mixing_dt;

  /*! Interparcel exchange parameter for mixing in the troposphere. */
  double mixing_trop;

  /*! Interparcel exchange parameter for mixing in the stratosphere. */
  double mixing_strat;

  /*! Number of altitudes of mixing grid. */
  int mixing_nz;

  /*! Lower altitude of mixing grid [km]. */
  double mixing_z0;

  /*! Upper altitude of mixing grid [km]. */
  double mixing_z1;

  /*! Number of longitudes of mixing grid. */
  int mixing_nx;

  /*! Lower longitude of mixing grid [deg]. */
  double mixing_lon0;

  /*! Upper longitude of mixing grid [deg]. */
  double mixing_lon1;

  /*! Number of latitudes of mixing grid. */
  int mixing_ny;

  /*! Lower latitude of mixing grid [deg]. */
  double mixing_lat0;

  /*! Upper latitude of mixing grid [deg]. */
  double mixing_lat1;

  /*! Number of altitudes of chemistry grid. */
  int chemgrid_nz;

  /*! Lower altitude of chemistry grid [km]. */
  double chemgrid_z0;

  /*! Upper altitude of chemistry grid [km]. */
  double chemgrid_z1;

  /*! Number of longitudes of chemistry grid. */
  int chemgrid_nx;

  /*! Lower longitude of chemistry grid [deg]. */
  double chemgrid_lon0;

  /*! Upper longitude of chemistry grid [deg]. */
  double chemgrid_lon1;

  /*! Number of latitudes of chemistry grid. */
  int chemgrid_ny;

  /*! Lower latitude of chemistry grid [deg]. */
  double chemgrid_lat0;

  /*! Upper latitude of chemistry grid [deg]. */
  double chemgrid_lat1;

  /*! Reaction type for OH chemistry (0=none, 2=bimolecular, 3=termolecular). */
  int oh_chem_reaction;

  /*! Coefficients for OH reaction rate (A, E/R or k0, n, kinf, m). */
  double oh_chem[4];

  /*! Beta parameter for diurnal variablity of OH. */
  double oh_chem_beta;

  /*! Reaction type for H2O2 chemistry (0=none, 1=SO2). */
  int h2o2_chem_reaction;

  /*! Switch for KPP chemistry module (0=off, 1=on). */
  int kpp_chem;

  /*! Time step for KPP chemistry [s]. */
  double dt_kpp;

  /*! Switch for first order tracer chemistry module (0=off, 1=on). */
  int tracer_chem;

  /*! Coefficients for precipitation calculation. */
  double wet_depo_pre[2];

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

  /*! Dry deposition surface layer [hPa]. */
  double dry_depo_dp;

  /*! Dry deposition velocity [m/s]. */
  double dry_depo_vdep;

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

  /*! Type of atmospheric data files
     (0=ASCII, 1=binary, 2=netCDF, 3=CLaMS). */
  int atm_type;

  /*! Basename of CSI data files. */
  char csi_basename[LEN];

  /*! Kernel data file for CSI output. */
  char csi_kernel[LEN];

  /*! Time step for CSI output [s]. */
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

  /*! Basename of ensemble data file. */
  char ens_basename[LEN];

  /*! Time step for ensemble output [s]. */
  double ens_dt_out;

  /*! Basename of grid data files. */
  char grid_basename[LEN];

  /*! Kernel data file for grid output. */
  char grid_kernel[LEN];

  /*! Gnuplot file for gridded data. */
  char grid_gpfile[LEN];

  /*! Time step for gridded data output [s]. */
  double grid_dt_out;

  /*! Sparse output in grid data files (0=no, 1=yes). */
  int grid_sparse;

  /*! Include standard deviations in grid output (0=no, 1=yes). */
  int grid_stddev;

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

  /*! Type of grid data files (0=ASCII, 1=netCDF). */
  int grid_type;

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

  /*! Basename of sample data file. */
  char sample_basename[LEN];

  /*! Kernel data file for sample output. */
  char sample_kernel[LEN];

  /*! Observation data file for sample output. */
  char sample_obsfile[LEN];

  /*! Horizontal radius for sample output [km]. */
  double sample_dx;

  /*! Layer depth for sample output [km]. */
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

  /*! Basename of VTK data files. */
  char vtk_basename[LEN];

  /*! Time step for VTK data output [s]. */
  double vtk_dt_out;

  /*! Particle index stride for VTK data. */
  int vtk_stride;

  /*! Vertical scaling factor for VTK data. */
  double vtk_scale;

  /*! Vertical offset for VTK data [km]. */
  double vtk_offset;

  /*! Spherical projection for VTK data (0=no, 1=yes). */
  int vtk_sphere;

} ctl_t;

/*! Atmospheric data. */
typedef struct {

  /*! Number of air parcels. */
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

/*! Climatological data in form of photolysis rates. */
typedef struct {

  /*! Number of pressure levels. */
  int np;

  /*! Number of solar zenith angles. */
  int nsza;

  /*! Number of total ozone columns. */
  int no3c;

  /*! Pressure [hPa]. */
  double p[CP];

  /*! Solar zenith angle [rad]. */
  double sza[CSZA];

  /*! Total column ozone [DU]. */
  double o3c[CO3];

  /*! N2O photolysis rate [1/s]. */
  double n2o[CP][CSZA][CO3];

  /*! CCl4 photolysis rate [1/s]. */
  double ccl4[CP][CSZA][CO3];

  /*! CCl3F photolysis rate [1/s]. */
  double ccl3f[CP][CSZA][CO3];

  /*! CCl2F2 photolysis rate [1/s]. */
  double ccl2f2[CP][CSZA][CO3];

  /*! O2 photolysis rate [1/s]. */
  double o2[CP][CSZA][CO3];

  /*! O3 photolysis rate [1/s]. o3 + hv = o1d + o2      */
  double o3_1[CP][CSZA][CO3];

  /*! O3 photolysis rate [1/s]. o3 + hv = o3p + o2      */
  double o3_2[CP][CSZA][CO3];

  /*! H2O2 photolysis rate [1/s]. */
  double h2o2[CP][CSZA][CO3];

} clim_photo_t;

/*! Climatological data in form of time series. */
typedef struct {

  /*! Number of timesteps. */
  int ntime;

  /*! Time [s]. */
  double time[CTS];

  /*! Volume mixing ratio [ppv]. */
  double vmr[CTS];

} clim_ts_t;

/*! Climatological data in form of zonal means. */
typedef struct {

  /*! Number of timesteps. */
  int ntime;

  /*! Number of latitudes. */
  int nlat;

  /*! Number of pressure levels. */
  int np;

  /*! Time [s]. */
  double time[CT];

  /*! Latitude [deg]. */
  double lat[CY];

  /*! Pressure [hPa]. */
  double p[CP];

  /*! Volume mixing ratio [ppv]. */
  double vmr[CT][CP][CY];

} clim_zm_t;

/*! Climatological data. */
typedef struct {

  /*! Number of tropopause timesteps. */
  int tropo_ntime;

  /*! Number of tropopause latitudes. */
  int tropo_nlat;

  /*! Tropopause time steps [s]. */
  double tropo_time[12];

  /*! Tropopause latitudes [deg]. */
  double tropo_lat[73];

  /*! Tropopause pressure values [hPa]. */
  double tropo[12][73];

  /*! Photolysis rates. */
  clim_photo_t photo;

  /*! HNO3 zonal means. */
  clim_zm_t hno3;

  /*! OH zonal means. */
  clim_zm_t oh;

  /*! H2O2 zonal means. */
  clim_zm_t h2o2;

  /*! HO2 zonal means. */
  clim_zm_t ho2;

  /*! O(1D) zonal means. */
  clim_zm_t o1d;

  /*! CFC-10 time series. */
  clim_ts_t ccl4;

  /*! CFC-11 time series. */
  clim_ts_t ccl3f;

  /*! CFC-12 time series. */
  clim_ts_t ccl2f2;

  /*! N2O time series. */
  clim_ts_t n2o;

  /*! SF6 time series. */
  clim_ts_t sf6;

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

  /*! Number of model levels. */
  int npl;

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

  /*! Land-sea mask [1]. */
  float lsm[EX][EY];

  /*! Sea surface temperature [K]. */
  float sst[EX][EY];

  /*! Boundary layer pressure [hPa]. */
  float pbl[EX][EY];

  /*! Tropopause pressure [hPa]. */
  float pt[EX][EY];

  /*! Tropopause temperature [K]. */
  float tt[EX][EY];

  /*! Tropopause geopotential height [km]. */
  float zt[EX][EY];

  /*! Tropopause water vapor volume mixing ratio [ppv]. */
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

  /*! Total column ozone [DU]. */
  float o3c[EX][EY];

  /*! Geopotential height [km]. */
  float z[EX][EY][EP];

  /*! Temperature [K]. */
  float t[EX][EY][EP];

  /*! Zonal wind [m/s]. */
  float u[EX][EY][EP];

  /*! Zonal wind on model levels [m/s]. */
  float ul[EX][EY][EP];

  /*! Meridional wind [m/s]. */
  float v[EX][EY][EP];

  /*! Meridional wind on model levels [m/s]. */
  float vl[EX][EY][EP];

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

  /*! Cloud cover [1]. */
  float cc[EX][EY][EP];

  /*! Pressure on model levels [hPa]. */
  float pl[EX][EY][EP];

  /*! Zeta [K]. */
  float zeta[EX][EY][EP];

  /*! Vertical velocity [K/s]. */
  float zeta_dot[EX][EY][EP];

  /*! Zeta on model levels [K]. */
  float zetal[EX][EY][EP];

  /*! Vertical velocity on model levels [K/s]. */
  float zeta_dotl[EX][EY][EP];

  /*! Hybrid model levels */
  double hybrid[EP];

#ifdef UVW
  /*! Cache for wind data. */
  float uvw[EX][EY][EP][3];
#endif

} met_t;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Convert Cartesian coordinates to geolocation. */
void cart2geo(
  const double *x,
  double *z,
  double *lon,
  double *lat);

/*! Check if x is finite. */
#ifdef _OPENACC
#pragma acc routine (check_finite)
#endif
int check_finite(
  const double x);

/*! Climatology of OH volume mixing ratios. */
#ifdef _OPENACC
#pragma acc routine (clim_oh)
#endif
double clim_oh(
  const ctl_t * ctl,
  const clim_t * clim,
  const double t,
  const double lon,
  const double lat,
  const double p);

/*! Initialization function for OH climatology. */
void clim_oh_diurnal_correction(
  ctl_t * ctl,
  clim_t * clim);

/*! Interpolate photolysis rate data. */
#ifdef _OPENACC
#pragma acc routine (clim_photo)
#endif
double clim_photo(
  double rate[CP][CSZA][CO3],
  clim_photo_t * photo,
  double p,
  double sza,
  double o3c);

/*! Climatology of tropopause pressure. */
#ifdef _OPENACC
#pragma acc routine (clim_tropo)
#endif
double clim_tropo(
  const clim_t * clim,
  const double t,
  const double lat);

/*! Initialize tropopause climatology. */
void clim_tropo_init(
  clim_t * clim);

/*! Interpolate time series. */
#ifdef _OPENACC
#pragma acc routine (clim_ts)
#endif
double clim_ts(
  const clim_ts_t * ts,
  const double t);

/*! Interpolate zonal mean climatology. */
#ifdef _OPENACC
#pragma acc routine (clim_zm)
#endif
double clim_zm(
  const clim_zm_t * zm,
  const double t,
  const double lat,
  const double p);

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
  const int year,
  const int mon,
  const int day,
  int *doy);

/*! Get date from day of year. */
void doy2day(
  const int year,
  const int doy,
  int *mon,
  int *day);

/*! Convert geolocation to Cartesian coordinates. */
void geo2cart(
  const double z,
  const double lon,
  const double lat,
  double *x);

/*! Get meteo data for given time step. */
void get_met(
  ctl_t * ctl,
  clim_t * clim,
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

/*! Spatiotemporal interpolation of meteo data. !*/
#ifdef _OPENACC
#pragma acc routine (intpol_met_4d_coord)
#endif
void intpol_met_4d_coord(
  met_t * met0,
  float height0[EX][EY][EP],
  float array0[EX][EY][EP],
  met_t * met1,
  float height1[EX][EY][EP],
  float array1[EX][EY][EP],
  double ts,
  double height,
  double lon,
  double lat,
  double *var,
  int *ci,
  double *cw,
  int init);

#ifdef UVW
#ifdef _OPENACC
#pragma acc routine (intpol_met_4d_coord_uvw)
#endif
void intpol_met_4d_coord_uvw(
  met_t * met0,
  float heights0[EX][EY][EP],
  float uvw0[EX][EY][EP][3],
  met_t * met1,
  float heights1[EX][EY][EP],
  float uvw1[EX][EY][EP][3],
  double ts,
  double height,
  double lon,
  double lat,
  double *u,
  double *v,
  double *zeta_dot,
  int *ci,
  double *cw,
  int init);
#endif

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
#ifdef UVW
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
#endif

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
#ifdef UVW
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
#endif

/*! Convert seconds to date. */
void jsec2time(
  const double jsec,
  int *year,
  int *mon,
  int *day,
  int *hour,
  int *min,
  int *sec,
  double *remain);

/*! Get weighting factor from kernel function. */
#ifdef _OPENACC
#pragma acc routine (kernel_weight)
#endif
double kernel_weight(
  const double kz[EP],
  const double kw[EP],
  const int nk,
  const double p);

/*! Calculate moist adiabatic lapse rate. */
#ifdef _OPENACC
#pragma acc routine (lapse_rate)
#endif
double lapse_rate(
  const double t,
  const double h2o);

/*! Get predefined pressure levels. */
void level_definitions(
  ctl_t * ctl);

/*! Find array index for irregular grid. */
#ifdef _OPENACC
#pragma acc routine (locate_irr)
#endif
int locate_irr(
  const double *xx,
  const int n,
  const double x);

/*! Find array index for regular grid. */
#ifdef _OPENACC
#pragma acc routine (locate_reg)
#endif
int locate_reg(
  const double *xx,
  const int n,
  const double x);

/*! Locate the four vertical indizes of a box for a given height value. */
#ifdef _OPENACC
#pragma acc routine (locate_vert)
#endif
void locate_vert(
  float profiles[EX][EY][EP],
  int np,
  int lon_ap_ind,
  int lat_ap_ind,
  double alt_ap,
  int *ind);

/*! locate the index in a column of a three dimensional array. */
#ifdef _OPENACC
#pragma acc routine (locate_irr_3d)
#endif
int locate_irr_3d(
  float profiles[EX][EY][EP],
  int np,
  int ind_lon,
  int ind_lat,
  double x);

/*! Calculate NAT existence temperature. */
#ifdef _OPENACC
#pragma acc routine (nat_temperature)
#endif
double nat_temperature(
  const double p,
  const double h2o,
  const double hno3);

/*! Read atmospheric data. */
int read_atm(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

/*! Read atmospheric data in ASCII format. */
int read_atm_asc(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

/*! Read atmospheric data in binary format. */
int read_atm_bin(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

/*! Read atmospheric data in CLaMS format. */
int read_atm_clams(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

/*! Read atmospheric data in netCDF format. */
int read_atm_nc(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

/*! Read climatological data. */
void read_clim(
  ctl_t * ctl,
  clim_t * clim);

/*! Read climatological time series. */
int read_clim_ts(
  const char *filename,
  clim_ts_t * ts);

/*! Read climatological zonal means. */
void read_clim_zm(
  const char *filename,
  char *varname,
  clim_zm_t * zm);

/*! Read climatological photolysis rates. */
void read_clim_photo(
  const char *filename,
  clim_photo_t * photo);

/*! Read control parameters. */
void read_ctl(
  const char *filename,
  int argc,
  char *argv[],
  ctl_t * ctl);

/*! Read kernel data file. */
void read_kernel(
  const char *filename,
  double kz[EP],
  double kw[EP],
  int *nk);

/*! Read meteo data file. */
int read_met(
  const char *filename,
  ctl_t * ctl,
  clim_t * clim,
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
  char *varname);

/*! Calculate convective available potential energy. */
void read_met_cape(
  clim_t * clim,
  met_t * met);

/*! Calculate cloud properties. */
void read_met_cloud(
  ctl_t * ctl,
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
  const char *filename,
  int ncid,
  ctl_t * ctl,
  met_t * met);

/*! Read meteo data on vertical levels. */
void read_met_levels(
  int ncid,
  ctl_t * ctl,
  met_t * met);

/*! Smooth vertical zeta and pressure profiles. */
void read_met_monotonize(
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
  ctl_t * ctl,
  met_t * met,
  float dest[EX][EY],
  float scl,
  int init);

/*! Read and convert 3D variable from meteo data file. */
int read_met_nc_3d(
  int ncid,
  char *varname,
  char *varname2,
  ctl_t * ctl,
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

/*! Fix polar winds. */
void read_met_polar_winds(
  met_t * met);

/*! Calculate potential vorticity. */
void read_met_pv(
  met_t * met);

/*! Calculate total column ozone. */
void read_met_ozone(
  met_t * met);

/*! Downsampling of meteo data. */
void read_met_sample(
  ctl_t * ctl,
  met_t * met);

/*! Read surface data. */
void read_met_surface(
  int ncid,
  met_t * met,
  ctl_t * ctl);

/*! Calculate tropopause data. */
void read_met_tropo(
  ctl_t * ctl,
  clim_t * clim,
  met_t * met);

/*! Read observation data. */
void read_obs(
  const char *filename,
  double *rt,
  double *rz,
  double *rlon,
  double *rlat,
  double *robs,
  int *nobs);

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
  const double p,
  const double T,
  const double rp,
  const double rhop);

/*! Spline interpolation. */
void spline(
  const double *x,
  const double *y,
  const int n,
  const double *x2,
  double *y2,
  const int n2,
  const int method);

/*! Calculate standard deviation. */
#ifdef _OPENACC
#pragma acc routine (stddev)
#endif
float stddev(
  const float *data,
  const int n);

/*! Calculate solar zenith angle. */
#ifdef _OPENACC
#pragma acc routine (sza_calc)
#endif
double sza_calc(
  const double sec,
  const double lon,
  const double lat);

/*! Convert date to seconds. */
void time2jsec(
  const int year,
  const int mon,
  const int day,
  const int hour,
  const int min,
  const int sec,
  const double remain,
  double *jsec);

/*! Measure wall-clock time. */
void timer(
  const char *name,
  const char *group,
  int output);

/*! Extract time information from filename. */
double time_from_filename(
  const char *filename,
  int offset);

/*! Get weighting factor based on tropopause distance. */
#ifdef _OPENACC
#pragma acc routine (tropo_weight)
#endif
double tropo_weight(
  const clim_t * clim,
  const double t,
  const double lat,
  const double p);

/*! Write atmospheric data. */
void write_atm(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm,
  double t);

/*! Write atmospheric data in ASCII format. */
void write_atm_asc(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm,
  double t);

/*! Write atmospheric data in binary format. */
void write_atm_bin(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

/*! Write atmospheric data in CLaMS position file format. */
void write_atm_clams(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

/*! Write atmospheric data in CLaMS position file and trajectory format */
void write_atm_clams_traj(
  const char *dirname,
  ctl_t * ctl,
  atm_t * atm,
  double t);

/*! Write atmospheric data in netCDF format. */
void write_atm_nc(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

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

/*! Write gridded data in ASCII format. */
void write_grid_asc(
  const char *filename,
  ctl_t * ctl,
  double *cd,
  double *mean[NQ],
  double *stddev[NQ],
  double *vmr_impl,
  double t,
  double *z,
  double *lon,
  double *lat,
  double *area,
  double dz,
  int *np);

/*! Write gridded data in netCDF format. */
void write_grid_nc(
  const char *filename,
  ctl_t * ctl,
  double *cd,
  double *mean[NQ],
  double *stddev[NQ],
  double *vmr_impl,
  double t,
  double *z,
  double *lon,
  double *lat,
  double *area,
  double dz,
  int *np);

/*! Read meteo data file. */
int write_met(
  const char *filename,
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

/*! Write VTK data. */
void write_vtk(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm,
  double t);

#endif /* LIBTRAC_H */
