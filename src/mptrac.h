/*
  This file is part of MPTRAC.
  
  MPTRAC is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
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
  
  \section Introduction

  The source code of MPTRAC is available from the
  [git repository](https://github.com/slcs-jsc/mptrac). Please see the
  [README.md](https://github.com/slcs-jsc/mptrac/blob/master/README.md)
  in the git repository for introductory information. More information
  can be found in the [user manual](https://slcs-jsc.github.io/mptrac).
  
  This doxygen manual contains information about the algorithms and
  data structures used in the code. Please refer to the `mptrac.h'
  documentation for a first overview.
  
  \section References
  
  For citing the model in scientific publications, please see
  [CITATION.cff](https://github.com/slcs-jsc/mptrac/blob/master/CITATION.cff)
  and refer to the following papers:
  
  _Hoffmann, L., Baumeister, P. F., Cai, Z., Clemens, J., Griessbach,
  S., Günther, G., Heng, Y., Liu, M., Haghighi Mood, K., Stein, O.,
  Thomas, N., Vogel, B., Wu, X., and Zou, L.: Massive-Parallel
  Trajectory Calculations version 2.2 (MPTRAC-2.2): Lagrangian
  transport simulations on graphics processing units (GPUs),
  Geosci. Model Dev., 15, 2731–2762,
  https://doi.org/10.5194/gmd-15-2731-2022, 2022._
  
  _Hoffmann, L., T. Rößler, S. Griessbach, Y. Heng, and O. Stein,
  Lagrangian transport simulations of volcanic sulfur dioxide
  emissions: Impact of meteorological data products,
  J. Geophys. Res. Atmos., 121, 4651-4673,
  https://doi.org/10.1002/2015JD023749, 2016._

  Additional references are collected here:
  https://slcs-jsc.github.io/mptrac/references
  
  \section License
  
  MPTRAC is being develop at the Jülich Supercomputing Centre,
  Forschungszentrum Jülich, Germany.

  MPTRAC is distributed under the terms of the
  [GNU General Public License v3.0](https://github.com/slcs-jsc/mptrac/blob/master/COPYING).
  
  \section Contributing

  We are interested in supporting operational and research
  applications with MPTRAC.
  
  You can submit bug reports or feature requests on the
  [issue tracker](https://github.com/slcs-jsc/mptrac/issues).
  
  Proposed code changes and fixes can be submitted as
  [pull requests](https://github.com/slcs-jsc/mptrac/pulls).

  Please do not hesitate to contact us if you have any questions or
  need assistance.
  
  \section Contact
  
  Dr. Lars Hoffmann
  
  Jülich Supercomputing Centre, Forschungszentrum Jülich
  
  e-mail: <l.hoffmann@fz-juelich.de>
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
#endif

#ifdef CURAND
#include "curand.h"
#endif

#ifdef ZFP
#include "zfp.h"
#endif

#ifdef ZSTD
#include "zstd.h"
#endif

#ifdef CMS
#include "cmultiscale.h"
#endif

#ifdef KPP
#include "chem_Parameters.h"
#include "chem_Global.h"
#include "chem_Sparse.h"
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

/**
 * @brief Allocate memory for a pointer with error handling.
 *
 * This macro allocates memory for a pointer of a given type and size
 * using the `calloc` function.  It includes error handling to check
 * if memory allocation was successful.  If the code is being compiled
 * with OpenACC support (_OPENACC macro defined), it additionally
 * checks if the code is running on a GPU device, and if not, it
 * raises an error.
 *
 * @param ptr Pointer variable to be allocated.
 * @param type Data type of the pointer.
 * @param n Number of elements to allocate memory for.
 *
 * @note If the code is compiled without OpenACC support, the
 * conditional check for GPU device is skipped.
 *
 * @author Lars Hoffmann
 */
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

/**
 * @brief Macro for computing the linear index of a 2D array element.
 *
 * The `ARRAY_2D` macro computes the linear index of a 2D array
 * element based on the specified row index (`ix`), column index
 * (`iy`), and number of columns (`ny`).
 *
 * @param ix Integer representing the row index of the 2D array element.
 * @param iy Integer representing the column index of the 2D array element.
 * @param ny Integer representing the number of columns in the 2D array.
 * @return The computed linear index of the 2D array element.
 *
 * The macro computes the linear index using the formula: `(ix) * (ny) + (iy)`.
 * This formula assumes row-major storage, where elements of
 * each row are stored sequentially in memory.
 *
 * @author Lars Hoffmann
 */
#define ARRAY_2D(ix, iy, ny)			\
  ((ix) * (ny) + (iy))

/**
 * @brief Compute the linear index of a 3D array element.
 *
 * This macro computes the linear index of a 3D array element based on
 * the specified row index (`ix`), column index (`iy`), depth index
 * (`iz`), number of columns (`ny`), and number of depths (`nz`).
 *
 * @param ix Row index of the 3D array element.
 * @param iy Column index of the 3D array element.
 * @param ny Number of columns in the 3D array.
 * @param iz Depth index of the 3D array element.
 * @param nz Number of depths in the 3D array.
 * @return Linear index of the 3D array element.
 *
 * @author Lars Hoffmann
 */
#define ARRAY_3D(ix, iy, ny, iz, nz)		\
  (((ix)*(ny) + (iy)) * (nz) + (iz))

/**
 * @brief Calculate the Arrhenius rate constant.
 *
 * The Arrhenius equation is commonly used in chemical kinetics to
 * describe the temperature dependence of reaction rates. This macro
 * calculates the rate constant (k) based on the Arrhenius equation:
 *
 * \f[ k = a \times \exp( -b / T ), \f]
 *
 * where:
 *   - k is the rate constant.
 *   - a is the pre-exponential factor or frequency factor.
 *   - b is the activation energy.
 *   - T is the temperature in Kelvin.
 *
 * @param a Pre-exponential factor or frequency factor.
 * @param b Activation energy.
 * @param t Temperature in Kelvin.
 * @return Calculated rate constant based on the Arrhenius equation.
 *
 * @author Mingzhao Liu
 */
#define ARRHENIUS(a, b, t)			\
  ((a) * exp( -(b) / (t)))

/**
 * @brief Convert a longitude difference to a distance in the x-direction (east-west) at a specific latitude.
 *
 * This macro calculates the distance in the x-direction (east-west)
 * corresponding to a given longitude difference at a specific
 * latitude using the formula:
 *
 * \f[ dx = dlon \times \pi \times RE / 180 times \cos(lat), \f]
 *
 * where:
 *   - dx is the distance in the x-direction (east-west).
 *   - dlon is the difference in longitudes in degrees.
 *   - RE is the Earth's radius.
 *   - lat is the latitude in degrees.
 *
 * @param dlon Difference in longitudes in degrees.
 * @param lat Latitude in degrees.
 * @return Distance in the x-direction (east-west) corresponding to the given longitude difference at the specified latitude.
 *
 * @author Lars Hoffmann
 */
#define DEG2DX(dlon, lat)					\
  ((dlon) * M_PI * RE / 180. * cos((lat) / 180. * M_PI))

/**
 * @brief Convert a latitude difference to a distance in the y-direction (north-south).
 *
 * This macro calculates the distance in the y-direction (north-south)
 * corresponding to a given latitude difference using the formula:
 *
 * \f[ dy = dlat \times \pi \times RE / 180, \f]
 *
 * where:
 *   - dy is the distance in the y-direction (north-south).
 *   - dlat is the difference in latitudes in degrees.
 *   - RE is the Earth's radius.
 *
 * @param dlat Difference in latitudes in degrees.
 * @return Distance in the y-direction (north-south) corresponding to the given latitude difference.
 *
 * @author Lars Hoffmann
 */
#define DEG2DY(dlat)				\
  ((dlat) * M_PI * RE / 180.)

/**
 * @brief Convert a pressure difference to a height difference in the vertical direction.
 *
 * This macro calculates the change in height (altitude) corresponding
 * to a given pressure difference using the formula:
 *
 * \f[ dz = - (dp) \times H_0 / p \f]
 *
 * where:
 *   - dz is the change in height (altitude) in meters.
 *   - dp is the pressure difference in hPa.
 *   - H0 is a reference scale height in km.
 *   - p is the reference pressure in hPa.
 *
 * @param dp Pressure difference in hPa.
 * @param p Reference pressure in hPa.
 * @return Change in height (altitude) in kilometers corresponding to the given pressure difference.
 *
 * @author Lars Hoffmann
 */
#define DP2DZ(dp, p)				\
  (- (dp) * H0 / (p))

/**
 * @brief Convert a distance in kilometers to degrees longitude at a given latitude.
 *
 * This macro calculates the change in longitude in degrees
 * corresponding to a given distance in kilometers at a specified
 * latitude on the Earth's surface. It uses the formula:
 *
 * \f[ dlon = \frac{dx \times 180}{\pi \times RE \times \cos(lat)} \f]
 *
 * @param dx Distance in kilometers.
 * @param lat Latitude in degrees.
 * @return Change in longitude in degrees.
 *
 * @note The latitude must be in the range [-89.999, 89.999] degrees.
 * Otherwise, the macro return value will be zero. This avoids issues
 * with the singularities at the poles.
 *
 * @author Lars Hoffmann
 */
#define DX2DEG(dx, lat)						\
  (((lat) < -89.999 || (lat) > 89.999) ? 0			\
   : (dx) * 180. / (M_PI * RE * cos((lat) / 180. * M_PI)))

/**
 * @brief Convert a distance in kilometers to degrees latitude.
 *
 * This macro calculates the change in latitude in degrees
 * corresponding to a given distance in kilometers on the Earth's
 * surface. It uses the formula:
 *
 * \f[ dlat = \frac{dy \times 180}{\pi \times RE} \f]
 *
 * @param dy Distance in kilometers.
 * @return Change in latitude in degrees.
 *
 * @author Lars Hoffmann
 */
#define DY2DEG(dy)				\
  ((dy) * 180. / (M_PI * RE))

/**
 * @brief Convert a change in altitude to a change in pressure.
 *
 * This macro calculates the change in pressure corresponding to a
 * given change in altitude.  It uses the hydrostatic equation:
 *
 * \f[ dp = -\left(dz \times \frac{p}{H_0}\right) \f]
 *
 * @param dz Change in altitude in kilometers.
 * @param p Current pressure in hPa.
 * @return Change in pressure in hPa.
 *
 * @author Lars Hoffmann
 */
#define DZ2DP(dz, p)				\
  (-(dz) * (p) / H0)

/**
 * @brief Calculate the distance between two points in Cartesian coordinates.
 *
 * This macro calculates the Euclidean distance between two points in
 * Cartesian coordinates.  It uses the square root of the square of
 * the distance obtained from the DIST2 macro.
 *
 * @param a Coordinates of the first point as an array of doubles.
 * @param b Coordinates of the second point as an array of doubles.
 * @return The distance between the two points.
 *
 * @author Lars Hoffmann
 */
#define DIST(a, b)				\
  sqrt(DIST2(a, b))

/**
 * @brief Calculate the squared Euclidean distance between two points in Cartesian coordinates.
 *
 * This macro calculates the squared Euclidean distance between two
 * points in Cartesian coordinates.  It computes the sum of the
 * squares of the differences of corresponding coordinates.
 *
 * @param a Coordinates of the first point as an array of doubles.
 * @param b Coordinates of the second point as an array of doubles.
 * @return The squared distance between the two points.
 *
 * @author Lars Hoffmann
 */
#define DIST2(a, b)                                                     \
  ((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))

/**
 * @brief Calculate the dot product of two vectors.
 *
 * This macro computes the dot product of two vectors represented as
 * arrays of doubles.  It multiplies corresponding components of the
 * vectors and sums the results.
 *
 * @param a The first vector as an array of doubles.
 * @param b The second vector as an array of doubles.
 * @return The dot product of the two vectors.
 *
 * @author Lars Hoffmann
 */
#define DOTP(a, b)				\
  (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

/**
 * @brief Calculate the floating-point remainder of dividing x by y.
 *
 * This macro computes the floating-point remainder of dividing x by
 * y.  It calculates this remainder as x minus the integer part of
 * (x / y) times y.
 *
 * @param x The dividend.
 * @param y The divisor.
 * @return The floating-point remainder of x divided by y.
 *
 * @note Macro has been added as a substitute when a GPU version of
 * fmod() is missing.
 *
 * @author Lars Hoffmann
 */
#define FMOD(x, y)				\
  ((x) - (int) ((x) / (y)) * (y))

/**
 * @brief Read data from a file stream and store it in memory.
 *
 * This macro reads data of a specified type from the given input file
 * stream and stores it in the specified memory location.  It ensures
 * that the correct amount of data is read from the file stream, and
 * if not, it raises an error.
 *
 * @param ptr Pointer to the memory location where the data will be stored.
 * @param type Type of the data elements to be read.
 * @param size Number of elements to read.
 * @param in File stream from which to read the data.
 *
 * @author Lars Hoffmann
 */
#define FREAD(ptr, type, size, in) {					\
    if(fread(ptr, sizeof(type), size, in)!=size)			\
      ERRMSG("Error while reading!");					\
  }

/**
 * @brief Write data from memory to a file stream.
 *
 * This macro writes data of a specified type from the specified
 * memory location to the given output file stream.  It ensures that
 * the correct amount of data is written to the file stream, and if
 * not, it raises an error.
 *
 * @param ptr Pointer to the memory location containing the data to be written.
 * @param type Type of the data elements to be written.
 * @param size Number of elements to write.
 * @param out File stream to which the data will be written.
 *
 * @author Lars Hoffmann
 */
#define FWRITE(ptr, type, size, out) {					\
    if(fwrite(ptr, sizeof(type), size, out)!=size)			\
      ERRMSG("Error while writing!");					\
  }

/**
 * @brief Initialize arrays for interpolation.
 *
 * This macro initializes arrays used for interpolation. It sets the
 * weights `cw` and indices `ci` to zero.  These arrays are used
 * during interpolation to store the interpolation weights and
 * indices.
 *
 * @author Lars Hoffmann
 */
#define INTPOL_INIT						\
  double cw[4] = {0.0, 0.0, 0.0, 0.0}; int ci[3] = {0, 0, 0};

/**
 * @brief Perform 2D interpolation for a meteorological variable.
 *
 * This macro performs 2D interpolation for a given meteorological variable at a specific time and location.
 *
 * @param var The variable to interpolate.
 * @param init A flag indicating whether to initialize the interpolation arrays (`cw` and `ci`). Set to 1 for initialization, 0 otherwise.
 * @return The interpolated value of the variable `var`.
 *
 * @author Lars Hoffmann
 */
#define INTPOL_2D(var, init)						\
  intpol_met_time_2d(met0, met0->var, met1, met1->var,			\
		     atm->time[ip], atm->lon[ip], atm->lat[ip],		\
		     &var, ci, cw, init);

/**
 * @brief Perform 3D interpolation for a meteorological variable.
 *
 * This macro performs 3D interpolation for a given meteorological
 * variable at a specific time, pressure level, and location.
 *
 * @param var The variable to interpolate.
 * @param init A flag indicating whether to initialize the interpolation arrays (`cw` and `ci`). Set to 1 for initialization, 0 otherwise.
 * @return The interpolated value of the variable `var`.
 *
 * @author Lars Hoffmann
 */
#define INTPOL_3D(var, init)						\
  intpol_met_time_3d(met0, met0->var, met1, met1->var,			\
		     atm->time[ip], atm->p[ip],				\
		     atm->lon[ip], atm->lat[ip],			\
		     &var, ci, cw, init);

/**
 * @brief Interpolate multiple meteorological variables in space.
 *
 * This macro performs spatial interpolation for multiple
 * meteorological variables at a given pressure level, longitude, and
 * latitude.
 *
 * @param p The pressure level at which to interpolate the variables.
 * @param lon The longitude at which to interpolate the variables.
 * @param lat The latitude at which to interpolate the variables.
 *
 * @author Lars Hofmann
 */
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
    intpol_met_space_3d(met, met->rwc, p, lon, lat, &rwc, ci, cw, 0);	\
    intpol_met_space_3d(met, met->iwc, p, lon, lat, &iwc, ci, cw, 0);	\
    intpol_met_space_3d(met, met->swc, p, lon, lat, &swc, ci, cw, 0);	\
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

/**
 * @brief Interpolate multiple meteorological variables in time.
 *
 * This macro performs temporal interpolation for multiple
 * meteorological variables at a given time, pressure level,
 * longitude, and latitude.
 *
 * @param time The time at which to interpolate the variables.
 * @param p The pressure level at which to interpolate the variables.
 * @param lon The longitude at which to interpolate the variables.
 * @param lat The latitude at which to interpolate the variables.
 *
 * @author Lars Hoffmann
 */
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
    intpol_met_time_3d(met0, met0->rwc, met1, met1->rwc, time, p, lon, lat, &rwc, ci, cw, 0); \
    intpol_met_time_3d(met0, met0->iwc, met1, met1->iwc, time, p, lon, lat, &iwc, ci, cw, 0); \
    intpol_met_time_3d(met0, met0->swc, met1, met1->swc, time, p, lon, lat, &swc, ci, cw, 0); \
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

/**
 * @brief Calculate lapse rate.
 *
 * This macro calculates the lapse rate between two pressure levels
 * given their temperatures and pressures.
 * 
 * @param p1 Pressure at the first level (in hPa).
 * @param t1 Temperature at the first level (in K).
 * @param p2 Pressure at the second level (in hPa).
 * @param t2 Temperature at the second level (in K).
 * @return The lapse rate (in K/km).
 *
 * @author Lars Hoffmann
 */
#define LAPSE(p1, t1, p2, t2)						\
  (1e3 * G0 / RA * ((t2) - (t1)) / ((t2) + (t1))			\
   * ((p2) + (p1)) / ((p2) - (p1)))

/**
 * @brief Linear interpolation.
 *
 * This macro performs linear interpolation to estimate the value of y
 * at a given x based on two points (x0, y0) and (x1, y1).
 * 
 * @param x0 X-coordinate of the first point.
 * @param y0 Y-coordinate of the first point.
 * @param x1 X-coordinate of the second point.
 * @param y1 Y-coordinate of the second point.
 * @param x The x-coordinate at which to estimate the y-value.
 * @return The estimated y-value at the given x-coordinate.
 *
 * @author Lars Hoffmann
 */
#define LIN(x0, y0, x1, y1, x)			\
  ((y0)+((y1)-(y0))/((x1)-(x0))*((x)-(x0)))

/**
 * @brief Macro to determine the maximum of two values.
 *
 * This macro evaluates to the larger of its two arguments, `a` and
 * `b`.  It uses a ternary conditional operator to compare the values
 * of `a` and `b` and returns `a` if `a` is greater than `b`;
 * otherwise, it returns `b`.
 *
 * @param a The first value to compare. Can be of any type that supports comparison.
 * @param b The second value to compare. Can be of any type that supports comparison.
 *
 * @return The larger of the two values, `a` or `b`.
 *
 * @note Both `a` and `b` are evaluated twice. If `a` or `b` have side
 * effects (e.g., increment operators, function calls), the side
 * effects will occur more than once. This can lead to unexpected
 * behavior.
 *
 * @warning The macro does not perform type checking, so `a` and `b`
 * should be of compatible types to avoid potential issues with
 * comparison and return value.
 *
 * @author Lars Hoffmann
 */
#define MAX(a,b)				\
  (((a)>(b))?(a):(b))

/**
 * @brief Write header for meteorological data file.
 *
 * This macro writes a header to a meteorological data file, providing
 * information about the variables stored in the file and their
 * corresponding columns.
 *
 * @param out Pointer to the file stream where the header will be written.
 *
 * @author Lars Hoffmann
 */
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
    	  "# $26 = cloud rain water content [kg/kg]\n"			\
	  "# $27 = cloud ice water content [kg/kg]\n"			\
    	  "# $28 = cloud snow water content [kg/kg]\n"			\
	  "# $29 = cloud cover [1]\n"					\
	  "# $30 = total column cloud water [kg/m^2]\n");		\
  fprintf(out,								\
	  "# $31 = cloud top pressure [hPa]\n"				\
	  "# $32 = cloud bottom pressure [hPa]\n"			\
	  "# $33 = pressure at lifted condensation level (LCL) [hPa]\n"	\
	  "# $34 = pressure at level of free convection (LFC) [hPa]\n"	\
	  "# $35 = pressure at equilibrium level (EL) [hPa]\n"		\
	  "# $36 = convective available potential energy (CAPE) [J/kg]\n" \
	  "# $37 = convective inhibition (CIN) [J/kg]\n"		\
	  "# $38 = relative humidity over water [%%]\n"			\
	  "# $39 = relative humidity over ice [%%]\n"			\
	  "# $40 = dew point temperature [K]\n");			\
  fprintf(out,								\
	  "# $41 = frost point temperature [K]\n"			\
	  "# $42 = NAT temperature [K]\n"				\
	  "# $43 = HNO3 volume mixing ratio [ppv]\n"			\
	  "# $44 = OH volume mixing ratio [ppv]\n"			\
	  "# $45 = H2O2 volume mixing ratio [ppv]\n"			\
	  "# $46 = HO2 volume mixing ratio [ppv]\n"			\
	  "# $47 = O(1D) volume mixing ratio [ppv]\n"			\
	  "# $48 = boundary layer pressure [hPa]\n"			\
	  "# $49 = total column ozone [DU]\n"				\
	  "# $50 = number of data points\n");				\
  fprintf(out,								\
	  "# $51 = number of tropopause data points\n"			\
	  "# $52 = number of CAPE data points\n");

/**
 * @brief Macro to determine the minimum of two values.
 *
 * This macro evaluates to the smaller of its two arguments, `a` and
 * `b`.  It uses a ternary conditional operator to compare the values
 * of `a` and `b` and returns `a` if `a` is less than `b`; otherwise,
 * it returns `b`.
 *
 * @param a The first value to compare. Can be of any type that supports comparison.
 * @param b The second value to compare. Can be of any type that supports comparison.
 *
 * @return The smaller of the two values, `a` or `b`.
 *
 * @note Both `a` and `b` are evaluated twice. If `a` or `b` have side
 * effects (e.g., increment operators, function calls), the side
 * effects will occur more than once. This can lead to unexpected
 * behavior.
 *
 * @warning The macro does not perform type checking, so `a` and `b`
 * should be of compatible types to avoid potential issues with
 * comparison and return value.
 *
 * @author Lars Hoffmann
 */
#define MIN(a,b)				\
  (((a)<(b))?(a):(b))

/**
 * @brief Calculate the density of a gas molecule.
 *
 * This macro calculates the density of a gas molecule using the
 * provided pressure and temperature values.
 *
 * @param p Pressure of the gas in Pascals.
 * @param t Temperature of the gas in Kelvin.
 * @return Density of the gas molecule in kg/m^3.
 *
 * @author Lars Hoffmann
 */
#define MOLEC_DENS(p,t)			\
  (AVO * 1e-6 * ((p) * 100) / (RI * (t)))

/**
 * @brief Broadcasts a one-dimensional array to all processes in the MPI communicator.
 *
 * This macro wraps the `MPI_Bcast` function to simplify broadcasting
 * a one-dimensional array from the root process (rank 0) to all other
 * processes in the MPI communicator (MPI_COMM_WORLD).
 *
 * @param pointer Pointer to the beginning of the array to be broadcasted.
 * @param nx      Number of elements in the array.
 * @param type    MPI data type of the elements in the array.
 *
 * @note The root process is assumed to have the correct data in the
 *       array before calling this macro.  All other processes will
 *       receive the data from the root process.
 *
 * @warning Ensure that the array has enough space allocated on all
 *          processes to hold the broadcasted data.
 *
 * @author Lars Hoffmann
 */
#define MPI_BCAST_1D(pointer, nx, type)			\
  MPI_Bcast(pointer, nx, type, 0, MPI_COMM_WORLD);

/**
 * @brief Broadcasts a two-dimensional array to all processes in the MPI communicator.
 *
 * This macro wraps the `MPI_Bcast` function to broadcast a
 * two-dimensional array from the root process (rank 0) to all other
 * processes in the MPI communicator (MPI_COMM_WORLD).  It iterates
 * over the first dimension of the array and broadcasts each
 * sub-array.
 *
 * @param pointer Pointer to the beginning of the two-dimensional array to be broadcasted.
 *                It is assumed that the array is allocated as a contiguous block of memory.
 * @param nx      Number of elements in the first dimension of the array.
 * @param ny      Number of elements in the second dimension of the array.
 * @param type    MPI data type of the elements in the array.
 *
 * @note The root process is assumed to have the correct data in the
 *       array before calling this macro.  All other processes will
 *       receive the data from the root process.
 *
 * @warning Ensure that the array has enough space allocated on all
 *          processes to hold the broadcasted data.  The array should
 *          be allocated as a contiguous block of memory, or
 *          adjustments should be made to account for non-contiguous
 *          memory layouts.
 *
 * @author Lars Hoffmann
 */
#define MPI_BCAST_2D(pointer, nx, ny, type)				\
  for (int ix = 0; ix < nx; ix++)					\
    MPI_Bcast(pointer[ix], ny, type, 0, MPI_COMM_WORLD);

/**
 * @brief Broadcasts a three-dimensional array to all processes in the MPI communicator.
 *
 * This macro wraps the `MPI_Bcast` function to broadcast a
 * three-dimensional array from the root process (rank 0) to all other
 * processes in the MPI communicator (MPI_COMM_WORLD).  It iterates
 * over the first and second dimension of the array and broadcasts
 * each sub-array.
 *
 * @param pointer Pointer to the beginning of the three-dimensional array to be broadcasted.
 *                It is assumed that the array is allocated as a contiguous block of memory.
 * @param nx      Number of elements in the first dimension of the array.
 * @param ny      Number of elements in the second dimension of the array.
 * @param ny      Number of elements in the third dimension of the array.
 * @param type    MPI data type of the elements in the array.
 *
 * @note The root process is assumed to have the correct data in the
 *       array before calling this macro.  All other processes will
 *       receive the data from the root process.
 *
 * @warning Ensure that the array has enough space allocated on all
 *          processes to hold the broadcasted data.  The array should
 *          be allocated as a contiguous block of memory, or
 *          adjustments should be made to account for non-contiguous
 *          memory layouts.
 *
 * @author Lars Hoffmann
 */
#define MPI_BCAST_3D(pointer, nx, ny, nz, type)				\
  for (int ix = 0; ix < nx; ix++)					\
    for (int iy = 0; iy < ny; iy++)					\
      MPI_Bcast(pointer[ix][iy], nz, type, 0, MPI_COMM_WORLD);

/**
 * @brief Execute a NetCDF command and check for errors.
 *
 * This macro executes a NetCDF command and checks the result. If the
 * result indicates an error, it prints the error message using
 * ERRMSG.
 *
 * @param cmd NetCDF command to execute.
 *
 * @author Lars Hoffmann
 */
#define NC(cmd) {				     \
    int nc_result=(cmd);			     \
    if(nc_result!=NC_NOERR)			     \
      ERRMSG("%s", nc_strerror(nc_result));	     \
  }

/**
 * @brief Define a NetCDF variable with attributes.
 *
 * This macro defines a NetCDF variable with the specified name, data
 * type, dimensions, long name, and units. It also sets the
 * "long_name" and "units" attributes for the variable.
 *
 * @param varname Name of the variable.
 * @param type Data type of the variable.
 * @param ndims Number of dimensions for the variable.
 * @param dims Array of dimension IDs.
 * @param long_name Long name of the variable.
 * @param units Units of the variable.
 *
 * @author Lars Hoffmann
 */
#define NC_DEF_VAR(varname, type, ndims, dims, long_name, units) {	\
    NC(nc_def_var(ncid, varname, type, ndims, dims, &varid));		\
    NC(nc_put_att_text(ncid, varid, "long_name", strnlen(long_name, LEN), long_name)); \
    NC(nc_put_att_text(ncid, varid, "units", strnlen(units, LEN), units)); \
  }

/**
 * @brief Retrieve a double-precision variable from a NetCDF file.
 *
 * This macro retrieves a double-precision variable from a NetCDF
 * file. It first checks if the variable exists in the file and then
 * reads its data into the specified pointer. If the `force` parameter
 * is set to true, it forces the retrieval of the variable, raising an
 * error if the variable does not exist.  If `force` is false, it
 * retrieves the variable if it exists and issues a warning if it does
 * not.
 *
 * @param varname Name of the variable to retrieve.
 * @param ptr Pointer to the memory location where the data will be stored.
 * @param force Boolean flag indicating whether to force retrieval (true) or not (false).
 *
 * @author Lars Hoffmann
 */
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

/**
 * @brief Inquire the length of a dimension in a NetCDF file.
 *
 * This macro retrieves the length of a specified dimension from a
 * NetCDF file.  It checks if the length of the dimension is within a
 * specified range and assigns the length to the provided pointer. If
 * the length is outside the specified range, an error message is
 * raised.
 *
 * @param dimname Name of the dimension to inquire.
 * @param ptr Pointer to an integer where the dimension length will be stored.
 * @param min Minimum acceptable length for the dimension.
 * @param max Maximum acceptable length for the dimension.
 *
 * @author Lars Hoffmann
 */
#define NC_INQ_DIM(dimname, ptr, min, max) {		\
    int dimid; size_t naux;				\
    NC(nc_inq_dimid(ncid, dimname, &dimid));		\
    NC(nc_inq_dimlen(ncid, dimid, &naux));		\
    *ptr = (int)naux;					\
    if ((*ptr) < (min) || (*ptr) > (max))		\
      ERRMSG("Dimension %s is out of range!", dimname);	\
  }

/**
 * @brief Write double precision data to a NetCDF variable.
 *
 * This macro writes data to a specified NetCDF variable. It can
 * handle both full variable writes and hyperslab writes depending on
 * the `hyperslab` parameter. If `hyperslab` is true, the data is
 * written as a hyperslab; otherwise, the entire variable is written.
 *
 * @param varname Name of the NetCDF variable to write to.
 * @param ptr Pointer to the data to be written.
 * @param hyperslab Boolean indicating whether to write the data as a hyperslab.
 *
 * @author Lars Hoffmann
 */
#define NC_PUT_DOUBLE(varname, ptr, hyperslab) {		\
    NC(nc_inq_varid(ncid, varname, &varid));			\
    if(hyperslab) {						\
      NC(nc_put_vara_double(ncid, varid, start, count, ptr));	\
    } else {							\
      NC(nc_put_var_double(ncid, varid, ptr));			\
    }								\
  }

/**
 * @brief Write a float array to a NetCDF file.
 *
 * This macro writes a float array to a specified variable in a NetCDF
 * file.  Depending on the value of the hyperslab parameter, the data
 * can be written as a hyperslab or as a whole variable.
 *
 * @param varname Name of the variable to which the float array will be written.
 * @param ptr Pointer to the float array to be written.
 * @param hyperslab Boolean flag indicating if the data should be written as a hyperslab. 
 *        - If true, the data will be written as a hyperslab using the start and count arrays.
 *        - If false, the data will be written to the entire variable.
 *
 * @author Lars Hoffmann
 */
#define NC_PUT_FLOAT(varname, ptr, hyperslab) {			\
    NC(nc_inq_varid(ncid, varname, &varid));			\
    if(hyperslab) {						\
      NC(nc_put_vara_float(ncid, varid, start, count, ptr));	\
    } else {							\
      NC(nc_put_var_float(ncid, varid, ptr));			\
    }								\
  }

/**
 * @brief Write integer data to a NetCDF variable.
 *
 * This macro writes data to a specified NetCDF variable. It can
 * handle both full variable writes and hyperslab writes depending on
 * the `hyperslab` parameter. If `hyperslab` is true, the data is
 * written as a hyperslab; otherwise, the entire variable is written.
 *
 * @param varname Name of the NetCDF variable to write to.
 * @param ptr Pointer to the data to be written.
 * @param hyperslab Boolean indicating whether to write the data as a hyperslab.
 *
 * @author Lars Hoffmann
 */
#define NC_PUT_INT(varname, ptr, hyperslab) {			\
    NC(nc_inq_varid(ncid, varname, &varid));			\
    if(hyperslab) {						\
      NC(nc_put_vara_int(ncid, varid, start, count, ptr));	\
    } else {							\
      NC(nc_put_var_int(ncid, varid, ptr));			\
    }								\
  }

/**
 * @brief Add a text attribute to a NetCDF variable.
 *
 * This macro adds a text attribute to a specified NetCDF variable. It
 * first retrieves the variable ID using its name, then it attaches
 * the text attribute to the variable.
 *
 * @param varname Name of the NetCDF variable to which the attribute will be added.
 * @param attname Name of the attribute to be added.
 * @param text Text of the attribute to be added.
 *
 * @author Lars Hoffmann
 */
#define NC_PUT_ATT(varname, attname, text) {				\
    NC(nc_inq_varid(ncid, varname, &varid));				\
    NC(nc_put_att_text(ncid, varid, attname, strnlen(text, LEN), text)); \
  }

/**
 * @brief Add a global text attribute to a NetCDF file.
 *
 * This macro adds a text attribute to the global attributes of a
 * NetCDF file.  It directly attaches the attribute to the file,
 * rather than to a specific variable.
 *
 * @param attname Name of the global attribute to be added.
 * @param text Text of the attribute to be added.
 *
 * @author Lars Hoffmann
 */
#define NC_PUT_ATT_GLOBAL(attname, text)				\
  NC(nc_put_att_text(ncid, NC_GLOBAL, attname, strnlen(text, LEN), text));

/**
 * @brief Perform nearest-neighbor interpolation.
 *
 * This macro returns the value of the nearest neighbor (y0 or y1) for
 * a given x value.  It compares the distances between x and x0, and
 * between x and x1, and returns the y value corresponding to the
 * closer x value.
 *
 * @param x0 The x-coordinate of the first point.
 * @param y0 The y-coordinate of the first point.
 * @param x1 The x-coordinate of the second point.
 * @param y1 The y-coordinate of the second point.
 * @param x The x-coordinate for which the nearest neighbor is to be found.
 * @return The y-coordinate of the nearest neighbor (either y0 or y1).
 *
 * @author Lars Hoffmann
 */
#define NN(x0, y0, x1, y1, x)				\
  (fabs((x) - (x0)) <= fabs((x) - (x1)) ? (y0) : (y1))

/**
 * @brief Compute the norm (magnitude) of a vector.
 *
 * This macro computes the Euclidean norm (magnitude) of a vector `a`
 * using the dot product of the vector with itself. The vector is
 * assumed to have three components: a[0], a[1], and a[2].
 *
 * @param a The vector for which the norm is to be computed.
 * @return The Euclidean norm of the vector.
 *
 * @author Lars Hoffmann
 */
#define NORM(a)					\
  sqrt(DOTP(a, a))

/**
 * @brief Loop over particle indices with OpenACC acceleration.
 *
 * This macro defines a loop over particle indices from `ip0` to `ip1`
 * with optional checking of `dt`. If `_OPENACC` is defined, the loop
 * is accelerated using OpenACC directives. Otherwise, OpenMP
 * parallelization is used.
 * 
 * @param ip0 The starting index of the loop (inclusive).
 * @param ip1 The ending index of the loop (exclusive).
 * @param check_dt Flag indicating whether to check the array `dt` for non-zero values.
 * @param ... Optional pragma directives to be applied.
 *
 * @author Lars Hoffmann
 */
#ifdef _OPENACC
#define PARTICLE_LOOP(ip0, ip1, check_dt, ...)		\
  const int ip0_const = ip0;                            \
  const int ip1_const = ip1;                            \
  _Pragma(__VA_ARGS__)					\
  _Pragma("acc parallel loop independent gang vector")  \
  for (int ip = ip0_const; ip < ip1_const; ip++)        \
    if (!check_dt || dt[ip] != 0)
#else
#define PARTICLE_LOOP(ip0, ip1, check_dt, ...)		\
  const int ip0_const = ip0;                            \
  const int ip1_const = ip1;                            \
  _Pragma("omp parallel for default(shared)")           \
  for (int ip = ip0_const; ip < ip1_const; ip++)        \
    if (!check_dt || dt[ip] != 0)
#endif

/**
 * @brief Compute pressure at given altitude.
 *
 * This macro calculates the pressure at a given altitude using the
 * barometric formula.
 * 
 * @param z The altitude in kilometers.
 * @return The pressure in hPa at the given altitude.
 *
 * The barometric formula used for this calculation is:
 *
 * \f[ P(z) = P_0 \times e^{-(z / H_0)}, \f]
 *
 * where:
 * - \f$ P(z) \f$ is the pressure at altitude \f$ z \f$,
 * - \f$ P_0 \f$ is the standard pressure,
 * - \f$ H_0 \f$ is the scale height.
 *
 * Note: The constants \f$ P_0 \f$ and \f$ H_0 \f$ must be defined before using this macro.
 *
 * @author Lars Hoffmann
 */
#define P(z)					\
  (P0 * exp(-(z) / H0))

/**
 * @brief Compute saturation pressure over water.
 *
 * This macro calculates the saturation pressure over water based on
 * the WMO (2018) formula.
 * 
 * @param t The temperature in degrees Celsius.
 * @return The saturation pressure over water at the given temperature.
 *
 * The saturation pressure over water is calculated using the formula:
 *
 * \f[ P_{\textrm{sat}}(t) = 6.112 \times e^{17.62 \times \frac{(t - T_0)}{243.12 + (t - T_0)}}, \f]
 *
 * where:
 * - \f$ P_{\textrm{sat}}(t) \f$ is the saturation pressure over water at temperature \f$ t \f$,
 * - \f$ T_0 \f$ is the reference temperature (0°C).
 *
 * Note: The constants \f$ T_0 \f$ must be defined before using this macro.
 *
 * @author Lars Hoffmann
 */
#define PSAT(t)							\
  (6.112 * exp(17.62 * ((t) - T0) / (243.12 + (t) - T0)))

/**
 * @brief Compute saturation pressure over ice (WMO, 2018).
 *
 * This macro calculates the saturation pressure over ice based on the
 * WMO (2018) formula.
 * 
 * @param t The temperature in K.
 * @return The saturation pressure over ice at the given temperature.
 *
 * The saturation pressure over ice is calculated using the formula:
 *
 * \f[ P_{\textrm{ice}}(t) = 6.112 \times e^{22.46 \times \frac{(t - T_0)}{272.62 + (t - T_0)}}, \f]
 *
 * where:
 * - \f$ P_{\textrm{ice}}(t) \f$ is the saturation pressure over ice at temperature \f$ t \f$,
 * - \f$ T_0 \f$ is the reference temperature (0°C).
 *
 * Note: The constant \f$ T_0 \f$ must be defined before using this macro.
 *
 * @author Lars Hoffmann
 */
#define PSICE(t)						\
  (6.112 * exp(22.46 * ((t) - T0) / (272.62 + (t) - T0)))

/**
 * @brief Calculate partial water vapor pressure.
 *
 * This macro calculates the partial water vapor pressure using the
 * given total pressure and water vapor mixing ratio.
 * 
 * @param p The total pressure in hPa (hectopascals).
 * @param h2o The water vapor mixing ratio in ppv (parts per volume).
 * @return The partial water vapor pressure.
 *
 * The partial water vapor pressure is calculated using the formula:
 *
 * \f[ P_{\textrm{w}}(p, h_2o) = \frac{p \times \max(h_2o, 0.1 \times 10^{-6})}{1 + (1 - \epsilon) \times \max(h_2o, 0.1 \times 10^{-6})}, \f]
 *
 * where:
 * - \f$ P_{\textrm{w}}(p, h_2o) \f$ is the partial water vapor pressure,
 * - \f$ p \f$ is the total pressure in hPa,
 * - \f$ h_2o \f$ is the water vapor mixing ratio in ppv,
 * - \f$ \epsilon \f$ is the factor to account for saturation vapor pressure over water.
 *
 * Note: The constant \f$ \epsilon \f$ must be defined before using this macro.
 * 
 * @author Lars Hoffmann
 */
#define PW(p, h2o)							\
  ((p) * MAX((h2o), 0.1e-6) / (1. + (1. - EPS) * MAX((h2o), 0.1e-6)))

/**
 * @brief Compute relative humidity over water.
 *
 * This macro calculates the relative humidity over water using the
 * given total pressure, temperature, and water vapor mixing ratio.
 * 
 * @param p The total pressure in hPa.
 * @param t The temperature in K.
 * @param h2o The water vapor mixing ratio in ppv (parts per volume).
 * @return The relative humidity over water in percentage.
 *
 * The relative humidity over water is calculated using the formula:
 *
 * \f[ RH_{\textrm{w}}(p, t, h_2o) = \frac{P_{\textrm{w}}(p, h_2o)}{P_{\textrm{sat}}(t)} \times 100, \f]
 *
 * where:
 * - \f$ RH_{\textrm{w}}(p, t, h_2o) \f$ is the relative humidity over water,
 * - \f$ P_{\textrm{w}}(p, h_2o) \f$ is the partial water vapor pressure,
 * - \f$ P_{\textrm{sat}}(t) \f$ is the saturation pressure over water at the given temperature,
 * - \f$ p \f$ is the total pressure in hPa,
 * - \f$ t \f$ is the temperature in Kelvin,
 * - \f$ h_2o \f$ is the water vapor mixing ratio in ppv.
 *
 * Note: The macros PW() and PSAT() must be defined before using this macro.
 * 
 * @author Lars Hoffmann
 */
#define RH(p, t, h2o)				\
  (PW(p, h2o) / PSAT(t) * 100.)

/**
 * @brief Compute relative humidity over ice.
 *
 * This macro calculates the relative humidity over ice using the
 * given total pressure, temperature, and water vapor mixing ratio.
 * 
 * @param p The total pressure in hPa.
 * @param t The temperature in K.
 * @param h2o The water vapor mixing ratio in ppv (parts per volume).
 * @return The relative humidity over ice in percentage.
 *
 * The relative humidity over ice is calculated using the formula:
 *
 * \f[ RH_{\textrm{ice}}(p, t, h_2o) = \frac{P_{\textrm{w}}(p, h_2o)}{P_{\textrm{ice}}(t)} \times 100, \f]
 *
 * where:
 * - \f$ RH_{\textrm{ice}}(p, t, h_2o) \f$ is the relative humidity over ice,
 * - \f$ P_{\textrm{w}}(p, h_2o) \f$ is the partial water vapor pressure,
 * - \f$ P_{\textrm{ice}}(t) \f$ is the saturation pressure over ice at the given temperature,
 * - \f$ p \f$ is the total pressure in hPa,
 * - \f$ t \f$ is the temperature in Kelvin,
 * - \f$ h_2o \f$ is the water vapor mixing ratio in ppv.
 *
 * Note: The macros PW() and PSICE() must be defined before using this macro.
 * 
 * @author Lars Hoffmann
 */
#define RHICE(p, t, h2o)			\
  (PW(p, h2o) / PSICE(t) * 100.)

/**
 * @brief Compute density of air.
 *
 * This macro calculates the density of air using the given total
 * pressure and temperature.
 * 
 * @param p The total pressure in hPa.
 * @param t The temperature in K.
 * @return The density of air in kg/m^3.
 *
 * The density of air is calculated using the formula:
 *
 * \f[ \rho(p, t) = \frac{100 \times p}{R_a \times t}, \f]
 *
 * where:
 * - \f$ \rho(p, t) \f$ is the density of air,
 * - \f$ p \f$ is the total pressure in hPa,
 * - \f$ t \f$ is the temperature in Kelvin,
 * - \f$ R_a \f$ is the specific gas constant for dry air (287.05 J/(kg·K)).
 *
 * @author Lars Hoffmann
 */
#define RHO(p, t)				\
  (100. * (p) / (RA * (t)))

/**
 * @brief Roeth approximation formula for photolysis reactions.
 *
 * This macro calculates the rate of a photolysis reaction using the
 * Roeth approximation formula, which takes into account the solar
 * zenith angle (SZA).
 * 
 * @param a Coefficient 'a' in the Roeth formula.
 * @param b Coefficient 'b' in the Roeth formula.
 * @param c Coefficient 'c' in the Roeth formula.
 * @param sza The solar zenith angle in radians.
 * @return The rate of the photolysis reaction.
 *
 * The Roeth approximation formula for photolysis reactions is given by:
 *
 * \f[ ROETH\_PHOTOL(a, b, c, sza) = 
 *   \begin{cases} 
 *     a \times \exp\left(b \times \left(1 - \frac{1}{\cos(c \times sza)}\right)\right), \textrm{if } c \times sza < \frac{\pi}{2} \\
 *     0, \textrm{otherwise}
 *   \end{cases}
 * \f]
 * 
 * where:
 * - 'a', 'b', and 'c' are coefficients specific to the photolysis reaction.
 * - 'sza' is the solar zenith angle in radians.
 *
 * @author Mingzhao Liu
 */
#define ROETH_PHOTOL(a, b, c, sza)					\
  ((c)*(sza) < M_PI/2. ? (a) * exp((b) * (1 - 1/cos((c) * (sza)))) : 0)

/**
 * @brief Set atmospheric quantity value.
 *
 * This macro sets the value of a specific atmospheric quantity at a
 * given index 'ip'.  The macro first checks if the control index
 * 'ctl->qnt' is non-negative before assigning the value, ensuring
 * that the quantity index is valid.
 * 
 * @param qnt The index representing the atmospheric quantity to set.
 * @param val The value to set for the atmospheric quantity.
 *
 * @note The macro assumes the existence of structures 'ctl' and 'atm' containing the control indices
 *       and atmospheric data, respectively, and an index 'ip' representing the data point.
 *
 * @author Lars Hoffmann
 */
#define SET_ATM(qnt, val)			\
  if (ctl->qnt >= 0)				\
    atm->q[ctl->qnt][ip] = val;

/**
 * @brief Set atmospheric quantity index.
 *
 * This macro sets the index, long name, and unit of a specific
 * atmospheric quantity based on its name.  It compares the name
 * parameter with the name of the atmospheric quantity stored in
 * 'ctl->qnt_name'.  If a match is found, it assigns the index to
 * 'ctl->qnt', updates the long name, and updates the unit.
 * 
 * @param qnt The index representing the atmospheric quantity.
 * @param name The name of the atmospheric quantity.
 * @param longname The long name of the atmospheric quantity.
 * @param unit The unit of the atmospheric quantity.
 *
 * @note The macro assumes the existence of structures 'ctl' containing control information.
 *       It also assumes the presence of 'iq', representing the index of the atmospheric quantity.
 *
 * @author Lars Hoffmann
 */
#define SET_QNT(qnt, name, longname, unit)		\
  if (strcasecmp(ctl->qnt_name[iq], name) == 0) {	\
    ctl->qnt = iq;					\
    sprintf(ctl->qnt_longname[iq], longname);		\
    sprintf(ctl->qnt_unit[iq], unit);			\
  } else

/**
 * @brief Compute specific humidity from water vapor volume mixing ratio.
 *
 * This macro calculates the specific humidity from the water vapor
 * volume mixing ratio.  Specific humidity represents the ratio of the
 * mass of water vapor to the total mass of air and is dimensionless.
 * 
 * @param h2o The water vapor volume mixing ratio.
 * @return The specific humidity.
 *
 * @note The macro assumes that 'EPS' is defined and represents the ratio of the molecular weight of water vapor to the molecular weight of dry air.
 *
 * @author Lars Hoffmann
 */
#define SH(h2o)					\
  (EPS * MAX((h2o), 0.1e-6))

/**
 * @brief Compute the square of a value.
 *
 * This macro computes the square of the input value.
 * 
 * @param x The input value.
 * @return The square of the input value.
 *
 * @author Lars Hoffmann
 */
#define SQR(x)					\
  ((x)*(x))

/**
 * @brief Swap two values.
 *
 * This macro swaps the values of two variables of the specified type.
 * 
 * @param x The first variable to be swapped.
 * @param y The second variable to be swapped.
 * @param type The type of the variables.
 *
 * @author Lars Hoffmann
 */
#define SWAP(x, y, type)			\
  do {type tmp = x; x = y; y = tmp;} while(0);

/**
 * @brief Calculate dew point temperature.
 *
 * This macro computes the dew point temperature using the formula
 * provided by the World Meteorological Organization (WMO, 2018).
 * 
 * @param p The atmospheric pressure in hPa.
 * @param h2o The water vapor volume mixing ratio.
 * @return The dew point temperature in Kelvin.
 *
 * Formula:
 *
 * \f[ T_{\textrm{dew}} = T_0 + \frac{243.12 \times \ln\left(\frac{{P_W(p, h_{2}O)}}{6.112}\right)}{17.62 - \ln\left(\frac{{P_W(p, h_{2}O)}}{6.112}\right)} \f]
 *
 * where:
 * - \f$ T_{\textrm{dew}} \f$ is the dew point temperature.
 * - \f$ T_0 \f$ is the reference temperature in Kelvin (typically 273.15 K).
 * - \f$ P_W(p, h_{2}O) \f$ is the partial water vapor pressure.
 * 
 * @author Lars Hoffmann
 */
#define TDEW(p, h2o)				\
  (T0 + 243.12 * log(PW((p), (h2o)) / 6.112)	\
   / (17.62 - log(PW((p), (h2o)) / 6.112)))

/**
 * @brief Calculate frost point temperature (WMO, 2018).
 *
 * This macro computes the frost point temperature using the formula
 * provided by the World Meteorological Organization (WMO, 2018).
 * 
 * @param p The atmospheric pressure in hPa.
 * @param h2o The water vapor volume mixing ratio.
 * @return The frost point temperature in Kelvin.
 *
 * Formula:
 *
 * \f[ T_{\textrm{ice}} = T_0 + \frac{272.62 \times \ln\left(\frac{{P_W(p, h_{2}O)}}{6.112}\right)}{22.46 - \ln\left(\frac{{P_W(p, h_{2}O)}}{6.112}\right)} \f]
 *
 * where:
 * - \f$ T_{\textrm{ice}} \f$ is the frost point temperature.
 * - \f$ T_0 \f$ is the reference temperature in Kelvin (typically 273.15 K).
 * - \f$ P_W(p, h_{2}O) \f$ is the partial water vapor pressure.
 * 
 * @author Lars Hoffmann
 */
#define TICE(p, h2o)				\
  (T0 + 272.62 * log(PW((p), (h2o)) / 6.112)	\
   / (22.46 - log(PW((p), (h2o)) / 6.112)))

/**
 * @brief Compute potential temperature.
 *
 * This macro calculates the potential temperature of the atmosphere.
 * 
 * @param p The atmospheric pressure in hPa.
 * @param t The temperature in Kelvin.
 * @return The potential temperature in Kelvin.
 *
 * Formula:
 *
 * \f[ \theta = T \left( \frac{1000}{P} \right)^{0.286} \f]
 *
 * where:
 * - \f$ \theta \f$ is the potential temperature.
 * - \f$ T \f$ is the temperature in Kelvin.
 * - \f$ P \f$ is the atmospheric pressure in hPa.
 * 
 * @author Lars Hoffmann
 */
#define THETA(p, t)				\
  ((t) * pow(1000. / (p), 0.286))

/**
 * @brief Compute virtual potential temperature.
 *
 * This macro calculates the virtual potential temperature of the
 * atmosphere, which takes into account the effect of water vapor on
 * the atmosphere's buoyancy.
 * 
 * @param p The atmospheric pressure in hPa.
 * @param t The temperature in Kelvin.
 * @param h2o The water vapor volume mixing ratio (ppv).
 * @return The virtual potential temperature in Kelvin.
 *
 * Formula:
 *
 * The virtual potential temperature (\f$ \theta_v \f$) is computed as
 *
 * \f[ \theta_v = \theta \left( 1 + \frac{0.61 \times q}{\epsilon} \right), \f]
 *
 * where:
 * - \f$ \theta_v \f$ is the virtual potential temperature.
 * - \f$ \theta \f$ is the potential temperature.
 * - \f$ q \f$ is the specific humidity.
 * - \f$ \epsilon \f$ is the ratio of the molecular weight of water vapor to dry air.
 * 
 * @author Lars Hoffmann
 */
#define THETAVIRT(p, t, h2o)			\
  (TVIRT(THETA((p), (t)), MAX((h2o), 0.1e-6)))

/**
 * @brief Get string tokens.
 *
 * This macro extracts tokens from a given string, typically used for
 * parsing input lines.
 * 
 * @param line The input string containing tokens.
 * @param tok A pointer to the token string.
 * @param format The format string specifying the expected format of the token.
 * @param var The variable to store the parsed token value.
 *
 * The macro tokenizes the input line using space and tab characters
 * as delimiters.  It then parses each token according to the
 * specified format string and stores the parsed value in the provided
 * variable.
 * 
 * @author Lars Hoffmann
 */
#define TOK(line, tok, format, var) {					\
    if(((tok)=strtok((line), " \t"))) {					\
      if(sscanf(tok, format, &(var))!=1) continue;			\
    } else ERRMSG("Error while reading!");				\
  }

/**
 * @brief Compute virtual temperature.
 *
 * This macro calculates the virtual temperature of air given its
 * temperature and water vapor volume mixing ratio.
 * 
 * @param t The temperature of the air in Kelvin.
 * @param h2o The water vapor volume mixing ratio.
 * @return The virtual temperature of the air.
 *
 * The virtual temperature (T_v) is computed as the temperature (t)
 * multiplied by (1 + (1 - EPS) * max(h2o, 0.1e-6)), where EPS is the
 * ratio of the molar mass of water vapor to the molar mass of dry
 * air.
 * 
 * @note EPS is typically defined as 0.622.
 * 
 * @author Lars Hoffmann
 */
#define TVIRT(t, h2o)					\
  ((t) * (1. + (1. - EPS) * MAX((h2o), 0.1e-6)))

/**
 * @brief Convert pressure to altitude.
 *
 * This macro calculates the altitude from the given pressure using
 * the barometric formula.
 * 
 * @param p The pressure in hPa (hectopascal).
 * @return The altitude in kilometers (km).
 *
 * Formula:
 *
 * The altitude (z) is computed as H0 times the natural logarithm of the
 * ratio of the reference pressure (P0) to the given pressure (p), where H0
 * is the scale height and P0 is the reference pressure at sea level.
 * 
 * @note H0 and P0 are typically defined as constants specific to the atmosphere.
 * 
 * @author Lars Hoffmann
 */
#define Z(p)					\
  (H0 * log(P0 / (p)))

/**
 * @brief Calculate geopotential height difference.
 *
 * This macro calculates the geopotential height difference between two
 * pressure levels using the hypsometric equation.
 * 
 * @param lnp0 The natural logarithm of the pressure at the first level.
 * @param t0 The temperature at the first level in Kelvin (K).
 * @param h2o0 The water vapor volume mixing ratio at the first level.
 * @param lnp1 The natural logarithm of the pressure at the second level.
 * @param t1 The temperature at the second level in Kelvin (K).
 * @param h2o1 The water vapor volume mixing ratio at the second level.
 * @return The geopotential height difference in kilometers (km).
 *
 * Formula:
 * The geopotential height difference (dz) is computed as a function of the
 * difference in natural logarithm of pressure (lnp) between the two levels,
 * the average virtual temperature (ThetaVirt) of the two levels, the specific
 * gas constant for dry air (RI), and the acceleration due to gravity at the
 * surface of the Earth (G0).
 * 
 * @note
 * The specific gas constant for dry air (RI), the molar mass of dry air (MA),
 * and the acceleration due to gravity at the surface of the Earth (G0) are
 * typically defined as constants specific to the atmosphere.
 * 
 * @author Lars Hoffmann
 */
#define ZDIFF(lnp0, t0, h2o0, lnp1, t1, h2o1)				\
  (RI / MA / G0 * 0.5 * (TVIRT((t0), (h2o0)) + TVIRT((t1), (h2o1)))	\
   * ((lnp0) - (lnp1)))

/**
 * @brief Calculate potential vorticity using the Zeta approximation.
 *
 * This macro calculates the potential vorticity using the Zeta
 * approximation, which is a function of pressure, surface pressure,
 * and temperature.
 * 
 * @param ps Surface pressure in hPa.
 * @param p Pressure at the given level in hPa.
 * @param t Temperature at the given level in Kelvin (K).
 * @return The potential vorticity.
 *
 * Formula:
 * The potential vorticity (Zeta) is computed based on the given
 * pressure (p), surface pressure (ps), and temperature (t). It
 * involves the potential temperature (Theta) at the given pressure
 * level.
 * 
 * @note
 * The potential vorticity is approximated using the Zeta
 * approximation, which includes conditions based on pressure ratio
 * and a sinusoidal function.  Adjust the constants if the
 * approximation or conditions differ from the standard Zeta
 * approximation formula.
 * 
 * @author Lars Hoffmann
 */
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

/*!
 * \brief Print a log message with a specified logging level.
 *
 * This macro prints a formatted log message to the standard output if
 * the specified logging level meets certain conditions. The message
 * will be indented if the logging level is greater than or equal to
 * 2.
 * 
 * \param level The logging level of the message. This should be an integer value.
 * \param ... The formatted message string and its arguments, similar to printf.
 *
 * \details
 * The `LOG` macro provides a simple way to log messages with
 * different levels of importance. The message is only printed if the
 * specified `level` is less than or equal to the pre-defined `LOGLEV`
 * macro. If the `level` is greater than or equal to 2, the message is
 * preceded by two spaces for indentation.
 *
 * The macro expands to a block of code that:
 * - Checks if the `level` is greater than or equal to 2, and if so, prints two spaces.
 * - Checks if the `level` is less than or equal to `LOGLEV`, and if so, prints the
 *   formatted message followed by a newline.
 *
 * \note
 * The `LOGLEV` macro must be defined with an appropriate logging level
 * before using the `LOG` macro.
 * 
 * @author Lars Hoffmann
 */
#define LOG(level, ...) {						\
    if(level >= 2)							\
      printf("  ");							\
    if(level <= LOGLEV) {						\
      printf(__VA_ARGS__);						\
      printf("\n");							\
    }									\
  }

/*!
 * \brief Print a warning message with contextual information.
 *
 * This macro prints a formatted warning message to the standard
 * output, including the file name, function name, and line number
 * where the warning occurred. The message is then passed to the `LOG`
 * macro with a logging level of 0.
 * 
 * \param ... The formatted warning message string and its arguments, similar to printf.
 *
 * \details
 * The `WARN` macro is used to print warning messages with additional context
 * about where the warning was triggered. The message includes the following
 * contextual information:
 * - The name of the source file where the macro is called (`__FILE__`).
 * - The name of the function where the macro is called (`__func__`).
 * - The line number in the source file where the macro is called (`__LINE__`).
 *
 * After printing this contextual information, the macro uses the
 * `LOG` macro with a logging level of 0 to print the actual warning
 * message. This ensures that warning messages are always logged,
 * regardless of the value of `LOGLEV`.
 *
 * \note
 * The `LOG` macro must be defined before using the `WARN` macro.
 * 
 * @author Lars Hoffmann
 */
#define WARN(...) {							\
    printf("\nWarning (%s, %s, l%d): ", __FILE__, __func__, __LINE__);	\
    LOG(0, __VA_ARGS__);						\
  }

/*!
 * \brief Print an error message with contextual information and terminate the program.
 *
 * This macro prints a formatted error message to the standard output,
 * including the file name, function name, and line number where the
 * error occurred. After printing the message, the program is
 * terminated with an exit status indicating failure.
 * 
 * \param ... The formatted error message string and its arguments, similar to printf.
 *
 * \details
 * The `ERRMSG` macro is used to report critical errors that require the
 * program to terminate immediately. The message includes the following
 * contextual information:
 * - The name of the source file where the macro is called (`__FILE__`).
 * - The name of the function where the macro is called (`__func__`).
 * - The line number in the source file where the macro is called (`__LINE__`).
 *
 * After printing this contextual information, the macro uses the
 * `LOG` macro with a logging level of 0 to print the actual error
 * message. Finally, the program exits with a failure status
 * (`EXIT_FAILURE`).
 *
 * \note
 * The `LOG` macro must be defined before using the `ERRMSG` macro.
 * 
 * @author Lars Hoffmann
 */
#define ERRMSG(...) {							\
    printf("\nError (%s, %s, l%d): ", __FILE__, __func__, __LINE__);	\
    LOG(0, __VA_ARGS__);						\
    exit(EXIT_FAILURE);							\
  }

/*!
 * \brief Print the value of a variable with contextual information.
 *
 * This macro prints the value of a variable to the standard output,
 * including the file name, function name, and line number where the
 * macro is called. The output also includes the variable's name and
 * value in a formatted string.
 * 
 * \param format The format string used to print the variable's value, similar to printf.
 * \param var The variable to be printed.
 *
 * \details
 * The `PRINT` macro is used to output the value of a variable along with
 * additional context about where the macro is called. The message includes:
 * - The name of the source file where the macro is called (`__FILE__`).
 * - The name of the function where the macro is called (`__func__`).
 * - The line number in the source file where the macro is called (`__LINE__`).
 * - The name of the variable being printed (`#var`).
 * - The value of the variable, formatted according to the provided format string (`format`).
 *
 * This macro is particularly useful for debugging purposes, providing
 * a convenient way to trace variable values and their locations in
 * the code.
 *
 * \note
 * The format string must be compatible with the type of the variable being printed.
 * 
 * @author Lars Hoffmann
 */
#define PRINT(format, var)						\
  printf("Print (%s, %s, l%d): %s= "format"\n",				\
	 __FILE__, __func__, __LINE__, #var, var);

/* ------------------------------------------------------------
   Timers...
   ------------------------------------------------------------ */

/*! Maximum number of timers. */
#define NTIMER 100

/*!
 * \brief Print the current state of all timers.
 *
 * This macro calls the `timer` function with predefined arguments to
 * signify the end of the timer logging process. It is used to print
 * the results of all the timers that have been tracked.
 *
 * \note
 * The `timer` function must be defined elsewhere in the codebase for this
 * macro to function correctly.
 * 
 * @author Lars Hoffmann
 */
#define PRINT_TIMERS				\
  timer("END", "END", 1);

/*!
 * \brief Select and start a timer with specific attributes.
 *
 * This macro stops the current timer (if any) and starts a new timer
 * with the specified ID, group, and color. It uses the `NVTX_POP` and
 * `NVTX_PUSH` macros for managing timer events and the `timer`
 * function to log the timer start event.
 *
 * \param id The identifier for the timer.
 * \param group The group name associated with the timer.
 * \param color The color code associated with the timer for NVTX visualization.
 *
 * \note
 * The `NVTX_POP`, `NVTX_PUSH`, and `timer` functions/macros must be defined
 * elsewhere in the codebase for this macro to function correctly.
 * 
 * @author Lars Hoffmann
 */
#define SELECT_TIMER(id, group, color) {				\
    NVTX_POP;								\
    NVTX_PUSH(id, color);						\
    timer(id, group, 0);						\
  }

/*!
 * \brief Starts a timer for tracking.
 *
 * This macro initializes the timer tracking process by pushing a
 * start event onto the stack using the `NVTX_PUSH` macro with a
 * predefined ID ("START") and color (`NVTX_CPU`).
 *
 * \note
 * The `NVTX_PUSH` macro must be defined elsewhere in the codebase for this
 * macro to function correctly.
 * 
 * @author Lars Hoffmann
 */
#define START_TIMERS				\
  NVTX_PUSH("START", NVTX_CPU);

/*!
 * \brief Stop the current timer.
 *
 * This macro stops the current timer by popping the top event from
 * the stack using the `NVTX_POP` macro.
 *
 * \note
 * The `NVTX_POP` macro must be defined elsewhere in the codebase for this
 * macro to function correctly.
 * 
 * @author Lars Hoffmann
 */
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

/*! Light green color code (MPI receive). */
#define NVTX_RECV 0xFFCCFFCB

/*! Dark green color code (MPI send). */
#define NVTX_SEND 0xFF008B00

/*!
 * \brief Macro for calling `nvtxRangePushEx` to start a named and colored NVTX range.
 *
 * This macro initializes an `nvtxEventAttributes_t` structure with
 * the provided title and color, then calls `nvtxRangePushEx` to mark
 * the beginning of an NVTX range.
 *
 * \param range_title The title of the NVTX range, displayed in the NVTX visual profiler.
 * \param range_color The color of the NVTX range, specified as an ARGB value.
 *
 * \details
 * The macro sets up the `nvtxEventAttributes_t` structure with the following fields:
 * - `version`: Set to `NVTX_VERSION`.
 * - `size`: Set to `NVTX_EVENT_ATTRIB_STRUCT_SIZE`.
 * - `messageType`: Set to `NVTX_MESSAGE_TYPE_ASCII` to indicate the message is an ASCII string.
 * - `colorType`: Set to `NVTX_COLOR_ARGB` to specify the color format.
 * - `color`: Set to the value of `range_color`.
 * - `message.ascii`: Set to the value of `range_title`.
 *
 * It then calls `nvtxRangePushEx` with the initialized attributes to start the NVTX range.
 *
 * \note
 * The NVTX (NVIDIA Tools Extension) library must be included and
 * initialized in your project for this macro to function
 * correctly. If NVTX is not available, an empty definition is
 * provided.
 * 
 * @author Lars Hoffmann
 */
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

/*!
 * \brief Macro for calling `nvtxRangePop` to end the current NVTX range.
 *
 * This macro calls `nvtxRangePop` to mark the end of the most
 * recently started NVTX range.
 *
 * \note
 * The NVTX (NVIDIA Tools Extension) library must be included and initialized in your project for
 * this macro to function correctly. If NVTX is not available, an empty definition is provided.
 * 
 * @author Lars Hoffmann
 */
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

/*!
 * \brief Wrapper to Thrust sorting function.
 *
 * This function serves as a wrapper for a Thrust sorting operation,
 * sorting the array `c` and maintaining the correspondence with the
 * `index` array.
 *
 * \param c Pointer to the array of double values to be sorted.
 * \param n The number of elements in the array `c`.
 * \param index Pointer to the array of indices, which will be updated to reflect the sorted order.
 *
 * \details
 * The `thrustSortWrapper` function uses the Thrust library to sort the array `c` in ascending order.
 * The `index` array is updated to reflect the new order of elements in `c` after sorting.
 *
 * This function is particularly useful when the sorted order of
 * elements needs to be tracked by indices.
 *
 * \note
 * - The `c` and `index` arrays must be of the same length `n`.
 * - The function assumes that the Thrust library is properly included and configured in the project.
 *
 * @author Kaveh Haghighi Mood
 */
void thrustSortWrapper(
  double *__restrict__ c,
  int n,
  int *__restrict__ index);

/* ------------------------------------------------------------
   Structs...
   ------------------------------------------------------------ */

/**
 * @brief Control parameters.
 * 
 * This structure contains all control parameters used by the MPTRAC
 * model.  The struct is used to collect to easily pass the control
 * parameters on to the various functions.
 */
typedef struct {

  /*! Coupled use of pressure based modules and diabatic advection. 
     (0= no coupling, 1= coupling) */
  int cpl_zeta_and_press_modules;

  /*! Use predefined pressure levels or not. */
  int press_level_def;

  /*! Vertical coordinate of air parcels (0=pressure, 1=zeta, 2=eta). */
  int vert_coord_ap;

  /*! Vertical coordinate of input meteo data (0=pressure-level, 1=model-level). */
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

  /*! Quantity array index for cloud rain water content. */
  int qnt_rwc;

  /*! Quantity array index for cloud ice water content. */
  int qnt_iwc;

  /*! Quantity array index for cloud snow water content. */
  int qnt_swc;

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

  /*! Quantity array index for total mass loss due to KPP chemistry. */
  int qnt_mloss_kpp;

  /*! Quantity array index for total mass loss due to wet deposition. */
  int qnt_mloss_wet;

  /*! Quantity array index for total mass loss due to dry deposition. */
  int qnt_mloss_dry;

  /*! Quantity array index for total mass loss due to exponential decay. */
  int qnt_mloss_decay;

  /*! Quantity array index for total loss rate. */
  int qnt_loss_rate;

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

  /*! ZFP compression precision for all variables, except z and T. */
  int met_zfp_prec;

  /*! ZFP compression tolerance for temperature. */
  double met_zfp_tol_t;

  /*! ZFP compression tolerance for geopotential height. */
  double met_zfp_tol_z;

  /*! cmultiscale coarsening heuristics
     (0=default, 1=mean diff, 2=median diff, 3=max diff). */
  int met_cms_heur;

  /*! cmultiscale compression epsilon for geopotential height. */
  double met_cms_eps_z;

  /*! cmultiscale compression epsilon for temperature. */
  double met_cms_eps_t;

  /*! cmultiscale compression epsilon for zonal wind. */
  double met_cms_eps_u;

  /*! cmultiscale compression epsilon for meridional wind. */
  double met_cms_eps_v;

  /*! cmultiscale compression epsilon for vertical velocity. */
  double met_cms_eps_w;

  /*! cmultiscale compression epsilon for potential vorticity. */
  double met_cms_eps_pv;

  /*! cmultiscale compression epsilon for water vapor. */
  double met_cms_eps_h2o;

  /*! cmultiscale compression epsilon for ozone. */
  double met_cms_eps_o3;

  /*! cmultiscale compression epsilon for cloud liquid water content. */
  double met_cms_eps_lwc;

  /*! cmultiscale compression epsilon for cloud rain water content. */
  double met_cms_eps_rwc;

  /*! cmultiscale compression epsilon for cloud ice water content. */
  double met_cms_eps_iwc;

  /*! cmultiscale compression epsilon for cloud snow water content. */
  double met_cms_eps_swc;

  /*! cmultiscale compression epsilon for cloud cover. */
  double met_cms_eps_cc;

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

  /*! Time step for sampling of meteo data along trajectories [s]. */
  double met_dt_out;

  /*! Preload meteo data into disk cache (0=no, 1=yes). */
  int met_cache;

  /*! Use MPI to share meteo (0=no, 1=yes). */
  int met_mpi_share;

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

  /*! Type of atmospheric data files for output
     (-1=same as ATM_TYPE, 0=netCDF, 1=binary, 2=pack, 3=zfp, 4=zstd). */
  int atm_type_out;

  /*! Type of observation data files
     (0=ASCII, 1=netCDF). */
  int obs_type;

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

/**
 * @brief Air parcel data.
 * 
 * This structure contains information related to air parcel data,
 * including the number of air parcels, their respective times,
 * pressures, longitudes, latitudes, and various user-defined
 * attributes.
 */
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

/**
 * @brief Cache data structure.
 * 
 * This structure contains data related to cached isosurface variables
 * and wind perturbations for a given set of data points.
 */
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

/**
 * @brief Climatological data in the form of photolysis rates.
 * 
 * This structure contains climatological data related to photolysis
 * rates at various pressure levels, solar zenith angles, and total
 * ozone columns.
 */
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

  /*! O3 photolysis rate (O3 + hv = O1d + O2) [1/s]. */
  double o3_1[CP][CSZA][CO3];

  /*! O3 photolysis rate (O3 + hv = O3p + O2) [1/s]. */
  double o3_2[CP][CSZA][CO3];

  /*! H2O2 photolysis rate [1/s]. */
  double h2o2[CP][CSZA][CO3];

  /*! H2O photolysis rate [1/s]. */
  double h2o[CP][CSZA][CO3];

} clim_photo_t;

/**
 * @brief Climatological data in the form of time series.
 * 
 * This structure contains climatological data in the form of time
 * series, representing the evolution of volume mixing ratio over
 * time.
 */
typedef struct {

  /*! Number of timesteps. */
  int ntime;

  /*! Time [s]. */
  double time[CTS];

  /*! Volume mixing ratio [ppv]. */
  double vmr[CTS];

} clim_ts_t;

/**
 * @brief Climatological data in the form of zonal means.
 * 
 * This structure contains climatological data organized as zonal
 * means, representing the distribution of volume mixing ratio over
 * latitudes and pressure levels across multiple timesteps.
 */
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

/**
 * @brief Climatological data.
 * 
 * This structure represents climatological data containing various
 * atmospheric parameters organized in different formats such as zonal
 * means, time series, photolysis rates, and tropopause data.
 */
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

/*!
 * @brief Meteo data structure.
 *
 * This structure holds meteorological data such as time, dimensions,
 * coordinates, surface properties, atmospheric profiles, and derived
 * variables of a given meteorological model.
 */
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

  /*! Pressure levels [hPa]. */
  double p[EP];

  /*! Model hybrid levels. */
  double hybrid[EP];

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

  /*! Cloud rain water content [kg/kg]. */
  float rwc[EX][EY][EP];

  /*! Cloud ice water content [kg/kg]. */
  float iwc[EX][EY][EP];

  /*! Cloud snow water content [kg/kg]. */
  float swc[EX][EY][EP];

  /*! Cloud cover [1]. */
  float cc[EX][EY][EP];

  /*! Pressure on model levels [hPa]. */
  float pl[EX][EY][EP];

  /*! Zonal wind on model levels [m/s]. */
  float ul[EX][EY][EP];

  /*! Meridional wind on model levels [m/s]. */
  float vl[EX][EY][EP];

  /*! Vertical velocity on model levels [hPa/s]. */
  float wl[EX][EY][EP];

  /*! Zeta on model levels [K]. */
  float zetal[EX][EY][EP];

  /*! Vertical velocity on model levels [K/s]. */
  float zeta_dotl[EX][EY][EP];

} met_t;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/**
 * @brief Converts Cartesian coordinates to geographic coordinates.
 *
 * This function converts a point from Cartesian coordinates (x, y, z)
 * to geographic coordinates (longitude, latitude, and altitude).  It
 * uses the spherical Earth approximation for the conversion.
 *
 * @param x Pointer to an array containing the Cartesian coordinates (x, y, z) in kilometers.
 * @param z Pointer to a double where the computed altitude (above the reference ellipsoid) will be stored, in kilometers.
 * @param lon Pointer to a double where the computed longitude (in degrees) will be stored.
 * @param lat Pointer to a double where the computed latitude (in degrees) will be stored.
 *
 * @author Lars Hoffmann
 */
void cart2geo(
  const double *x,
  double *z,
  double *lon,
  double *lat);

/**
 * @brief Calculates the hydroxyl radical (OH) concentration from climatology data, 
 *        with an optional diurnal correction based on solar zenith angle.
 *
 * This function retrieves OH data from a given climatology and
 * applies a diurnal correction if the correction factor
 * (oh_chem_beta) is greater than zero. The diurnal correction
 * accounts for the variation in OH concentration due to changes in
 * the solar zenith angle.
 *
 * @param ctl  Pointer to the control structure containing configuration parameters.
 * @param clim Pointer to the climatology structure containing OH data.
 * @param t    Time at which the OH concentration is to be calculated.
 * @param lon  Longitude at which the OH concentration is to be calculated.
 * @param lat  Latitude at which the OH concentration is to be calculated.
 * @param p    Pressure level at which the OH concentration is to be calculated.
 * @return     The OH concentration at the specified time, location, and pressure,
 *             possibly adjusted by a diurnal correction.
 *
 * @authors Lars Hoffmann
 * @authors Mingzhao Liu
 */
double clim_oh(
  const ctl_t * ctl,
  const clim_t * clim,
  const double t,
  const double lon,
  const double lat,
  const double p);

/**
 * @brief Applies a diurnal correction to the hydroxyl radical (OH) concentration 
 *        in climatology data.
 *
 * This function iterates over the climatology data points for OH
 * concentration and integrates the day/night correction factor over
 * longitude. The correction factor is based on the solar zenith
 * angle, and it adjusts the OH data to account for diurnal
 * variations. The corrected OH data is scaled accordingly.
 *
 * @param ctl  Pointer to the control structure containing configuration parameters,
 *             including the correction factor (oh_chem_beta).
 * @param clim Pointer to the climatology structure containing OH data that will be
 *             corrected.
 *
 * @authors Lars Hoffmann
 * @authors Mingzhao Liu
 */
void clim_oh_diurnal_correction(
  ctl_t * ctl,
  clim_t * clim);

/**
 * @brief Calculates the photolysis rate for a given set of atmospheric conditions.
 *
 * This function computes the photolysis rate based on provided
 * climatology data and input parameters such as pressure, solar
 * zenith angle (SZA), and ozone column. It ensures that the input
 * parameters are within the valid range of the climatology data and
 * interpolates the photolysis rate accordingly.
 *
 * @param rate 3D array containing the photolysis rates for different combinations 
 *             of pressure, SZA, and ozone column.
 * @param photo Pointer to the climatology data structure containing arrays of valid 
 *              pressure levels, SZAs, and ozone columns.
 * @param p    Pressure at which the photolysis rate is to be calculated.
 * @param sza  Solar zenith angle at which the photolysis rate is to be calculated.
 * @param o3c  Ozone column at which the photolysis rate is to be calculated.
 * @return     The interpolated photolysis rate for the specified conditions. If the 
 *             calculated rate is negative, it returns 0.0.
 *
 * This function performs the following steps:
 * 1. Checks and adjusts the input parameters (pressure, SZA, and ozone column) to 
 *    ensure they are within the valid range.
 * 2. Determines the appropriate indices in the climatology data for interpolation.
 * 3. Performs trilinear interpolation to calculate the photolysis rate based on 
 *    the input parameters.
 *
 * @authors Lars Hoffmann
 * @authors Mingzhao Liu
 */
double clim_photo(
  double rate[CP][CSZA][CO3],
  clim_photo_t * photo,
  double p,
  double sza,
  double o3c);

/**
 * @brief Calculates the tropopause pressure based on climatological data.
 *
 * This function computes the tropopause pressure using climatological
 * data for different times and latitudes. It interpolates the
 * tropopause pressure based on the input time and latitude
 * parameters.
 *
 * @param clim Pointer to the climatology structure containing tropopause 
 *             pressure data.
 * @param t    Time for which the tropopause pressure is to be calculated, 
 *             in seconds since the beginning of the year.
 * @param lat  Latitude at which the tropopause pressure is to be calculated.
 * @return     The interpolated tropopause pressure for the specified time 
 *             and latitude.
 *
 * This function performs the following steps:
 * 1. Calculates the number of seconds since the beginning of the year.
 * 2. Determines the appropriate indices in the climatology data for 
 *    interpolation based on time and latitude.
 * 3. Interpolates the tropopause pressure using linear interpolation 
 *    based on latitude and time.
 *
 * @author Lars Hoffmann
 */
double clim_tropo(
  const clim_t * clim,
  const double t,
  const double lat);

/**
 * @brief Initializes the tropopause data in the climatology structure.
 *
 * This function initializes the tropopause data in the climatology
 * structure.  It sets the time steps, latitudes, and tropopause
 * pressure values based on predefined arrays.
 *
 * @param clim Pointer to the climatology structure to be initialized.
 *
 * This function performs the following steps:
 * 1. Sets the number of time steps and initializes the time array.
 * 2. Sets the number of latitudes and initializes the latitude array.
 * 3. Initializes the tropopause pressure values based on predefined arrays.
 * 4. Computes the range of tropopause pressure values.
 * 5. Logs information about the initialization process.
 *
 * @author Lars Hoffmann
 */
void clim_tropo_init(
  clim_t * clim);

/**
 * @brief Interpolates a time series of climatological variables.
 *
 * This function interpolates a time series of climatological
 * variables based on the input time and the provided data points.
 *
 * @param ts Pointer to the time series structure containing data points.
 * @param t Time at which to interpolate the climatological variable (in seconds).
 * @return Interpolated value of the climatological variable at the given time.
 *
 * This function performs linear interpolation between the closest data points
 * to the input time `t`. If `t` is outside the range of the provided time
 * series, the value at the nearest boundary is returned.
 *
 * @author Lars Hoffmann
 */
double clim_ts(
  const clim_ts_t * ts,
  const double t);

/**
 * @brief Interpolates monthly mean zonal mean climatological variables.
 *
 * This function interpolates climatological variables based on
 * pressure, latitude, and time. The climatological data is provided
 * in the form of monthly mean zonal mean values.
 *
 * @param zm Pointer to the climatological zonal mean structure containing data points.
 * @param t Time at which to interpolate the climatological variable (in seconds since the beginning of the year).
 * @param lat Latitude at which to interpolate the climatological variable (in degrees).
 * @param p Pressure at which to interpolate the climatological variable (in hPa).
 * @return Interpolated value of the climatological variable at the given pressure, latitude, and time.
 *
 * This function performs trilinear interpolation between the nearest
 * data points to the input time `t`, latitude `lat`, and pressure or
 * altitude `p`. If the input values are outside the range of the
 * provided data, the function extrapolates by using the nearest
 * boundary values.
 *
 * @author Lars Hoffmann
 */
double clim_zm(
  const clim_zm_t * zm,
  const double t,
  const double lat,
  const double p);

/**
 * @brief Compresses or decompresses a 3D array of floats using a custom multiscale compression algorithm.
 *
 * This function either compresses or decompresses a 3D array of
 * floats based on the value of the `decompress` parameter.  The
 * compression and decompression are performed using a custom
 * multiscale module.
 *
 * @param varname The name of the variable being processed.
 * @param array Pointer to the 3D array of floats to be compressed or decompressed.
 * @param nx The number of elements in the x-dimension of the array.
 * @param ny The number of elements in the y-dimension of the array.
 * @param np The number of elements in the p-dimension of the array.
 * @param decompress If non-zero, the function will decompress the data; otherwise, it will compress the data.
 * @param inout File pointer for input or output operations. It is used for reading compressed data during decompression
 * and writing compressed data during compression.
 *
 * The function performs the following steps:
 * - Determines grid properties from the input data dimensions.
 * - Initializes the multiscale module with the specified grid properties.
 * - Sets up longitude and latitude grids for the data.
 * - If decompressing:
 *   - Reads compressed data for each level and decompresses it.
 *   - Evaluates the decompressed data and stores it in the `array`.
 * - If compressing:
 *   - Copies data for each level into a temporary array.
 *   - Compresses the data and writes the compressed data to the file.
 *
 * The function logs the compression or decompression details and
 * frees allocated resources before returning.
 *
 * @note Ensure that the input `array` is already allocated and can hold the decompressed data.
 * 
 * @warning Ensure that the file pointer `inout` is correctly opened for reading or writing as required.
 * 
 * @see get_2d_grid_from_meteo_data, init_multiscale, read_sol, save_sol, eval, coarsening, delete_solution, delete_multiscale
 * 
 * @author 
 * Lars Hoffmann
 */
void compress_cms(
  ctl_t * ctl,
  char *varname,
  float *array,
  size_t nx,
  size_t ny,
  size_t np,
  int decompress,
  FILE * inout);

/**
 * @brief Compresses or decompresses a 3D array of floats.
 *
 * This function either compresses or decompresses a 3D array of
 * floats based on the value of the `decompress`
 * parameter. Compression reduces the storage size by converting float
 * values (4 bytes) to unsigned short values (2 bytes) with scaling
 * and offset. Decompression restores the original float values from
 * the compressed unsigned short representation.
 *
 * @param varname The name of the variable being processed.
 * @param array Pointer to the 3D array of floats to be compressed or decompressed.
 * @param nxy The number of elements in the first two dimensions of the array.
 * @param nz The number of elements in the third dimension of the array.
 * @param decompress If non-zero, the function will decompress the data; otherwise, it will compress the data.
 * @param inout File pointer for input or output operations. It is used for reading compressed data during decompression 
 * and writing compressed data during compression.
 *
 * The function performs the following steps:
 * - If decompressing:
 *   - Reads scaling factors, offsets, and compressed data from the file.
 *   - Decompresses the data and stores it in the `array`.
 * - If compressing:
 *   - Computes the minimum and maximum values for each slice in the third dimension.
 *   - Calculates scaling factors and offsets based on these values.
 *   - Compresses the data by converting floats to unsigned shorts using the scaling factors and offsets.
 *   - Writes the scaling factors, offsets, and compressed data to the file.
 *
 * The function allocates memory for the compressed data array and frees it before returning.
 *
 * @author Lars Hoffmann
 */
void compress_pack(
  char *varname,
  float *array,
  size_t nxy,
  size_t nz,
  int decompress,
  FILE * inout);

/**
 * @brief Compresses or decompresses a 3D array of floats using the ZFP library.
 *
 * This function either compresses or decompresses a 3D array of
 * floats based on the value of the `decompress`
 * parameter. Compression reduces the storage size using the ZFP
 * compression algorithm, which supports fixed-precision or
 * fixed-accuracy modes. Decompression restores the original float
 * values from the compressed representation.
 *
 * @param varname The name of the variable being processed.
 * @param array Pointer to the 3D array of floats to be compressed or decompressed.
 * @param nx The number of elements in the x-dimension of the array.
 * @param ny The number of elements in the y-dimension of the array.
 * @param nz The number of elements in the z-dimension of the array.
 * @param precision The precision parameter for ZFP compression. If greater than 0, it sets the fixed precision mode.
 * @param tolerance The tolerance parameter for ZFP compression. If greater than 0 and precision is 0, it sets the fixed accuracy mode.
 * @param decompress If non-zero, the function will decompress the data; otherwise, it will compress the data.
 * @param inout File pointer for input or output operations. It is used for reading compressed data during decompression 
 * and writing compressed data during compression.
 *
 * The function performs the following steps:
 * - Allocates metadata for the 3D array and the ZFP compressed stream.
 * - Sets the compression mode based on the precision or tolerance parameters.
 * - Allocates a buffer for the compressed data.
 * - Associates a bit stream with the allocated buffer and sets up the ZFP stream.
 * - If decompressing:
 *   - Reads the size of the compressed data and the compressed data itself from the file.
 *   - Decompresses the data and stores it in the `array`.
 * - If compressing:
 *   - Compresses the data and writes the compressed data size and the compressed data itself to the file.
 *
 * The function logs the compression or decompression details and frees allocated resources before returning.
 *
 * @note Ensure that either the precision or tolerance parameter is set to a value greater than 0.
 *
 * @author Lars Hoffmann
 */
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

/**
 * @brief Compresses or decompresses an array of floats using the Zstandard (ZSTD) library.
 *
 * This function either compresses or decompresses an array of floats
 * based on the value of the `decompress` parameter.  Compression
 * reduces the storage size using the ZSTD compression
 * algorithm. Decompression restores the original float values from
 * the compressed representation.
 *
 * @param varname The name of the variable being processed.
 * @param array Pointer to the array of floats to be compressed or decompressed.
 * @param n The number of elements in the array.
 * @param decompress If non-zero, the function will decompress the data; otherwise, it will compress the data.
 * @param inout File pointer for input or output operations. It is used for reading compressed data during decompression
 * and writing compressed data during compression.
 *
 * The function performs the following steps:
 * - Calculates the buffer sizes required for compression and decompression.
 * - Allocates memory for the compressed data buffer.
 * - If decompressing:
 *   - Reads the size of the compressed data and the compressed data itself from the file.
 *   - Decompresses the data and stores it in the `array`.
 * - If compressing:
 *   - Compresses the data and writes the compressed data size and the compressed data itself to the file.
 *
 * The function logs the compression or decompression details and frees allocated resources before returning.
 *
 * @note This function assumes that the input `array` is already allocated and can hold the decompressed data.
 * 
 * @warning Ensure that the file pointer `inout` is correctly opened for reading or writing as required.
 *
 * @see ZSTD_compress, ZSTD_decompress, ZSTD_isError
 * 
 * @author Lars Hoffmann
 */
void compress_zstd(
  char *varname,
  float *array,
  size_t n,
  int decompress,
  FILE * inout);

/*! Get day of year from date. */
/**
 * @brief Converts a given date to the day of the year (DOY).
 *
 * This function computes the day of the year (DOY) for a given date
 * specified by the year, month, and day. It takes into account
 * whether the given year is a leap year or not.
 *
 * @param year The year of the date.
 * @param mon The month of the date (1-12).
 * @param day The day of the month (1-31).
 * @param doy Pointer to an integer where the computed day of the year will be stored.
 *
 * The function uses two arrays, `d0` and `d0l`, which contain the
 * cumulative number of days at the start of each month for non-leap
 * years and leap years respectively. It checks if the year is a leap
 * year and calculates the day of the year accordingly.
 *
 * @note The function assumes that the input date is valid.
 * 
 * @author Lars Hoffmann
 */
void day2doy(
  const int year,
  const int mon,
  const int day,
  int *doy);

/**
 * @brief Converts a given day of the year (DOY) to a date (month and day).
 *
 * This function computes the month and day for a given day of the
 * year (DOY) and year. It accounts for whether the given year is a
 * leap year or not.
 *
 * @param year The year corresponding to the DOY.
 * @param doy The day of the year (1-365 or 1-366).
 * @param mon Pointer to an integer where the computed month will be stored.
 * @param day Pointer to an integer where the computed day of the month will be stored.
 *
 * The function uses two arrays, `d0` and `d0l`, which contain the
 * cumulative number of days at the start of each month for non-leap
 * years and leap years respectively. It checks if the year is a leap
 * year and calculates the month and day of the month accordingly.
 *
 * @note The function assumes that the input DOY is valid for the given year.
 * 
 * @author Lars Hoffmann
 */
void doy2day(
  const int year,
  const int doy,
  int *mon,
  int *day);

/**
 * @brief Computes the Fast Fourier Transform (FFT) of a complex sequence.
 *
 * This function calculates the FFT of a complex sequence represented
 * by separate arrays for the real and imaginary parts. The input
 * arrays `fcReal` and `fcImag` are modified in place to contain the
 * transformed data.
 *
 * @param fcReal Pointer to an array of doubles representing the real part of the input sequence. 
 *               The array should have at least `n` elements.
 * @param fcImag Pointer to an array of doubles representing the imaginary part of the input sequence. 
 *               The array should have at least `n` elements.
 * @param n The number of complex data points in the input sequence. This value should not exceed PMAX.
 *
 * @pre `fcReal` and `fcImag` must point to arrays of at least `n` elements.
 * @pre `n` must be less than or equal to PMAX.
 *
 * @post The arrays `fcReal` and `fcImag` will contain the real and imaginary parts of the FFT result, respectively.
 *
 * @note This function uses the GNU Scientific Library (GSL) for computing the FFT. Ensure that GSL is properly installed 
 *       and linked in your project.
 *
 * @warning If `n` exceeds PMAX, the function will trigger an error message and terminate.
 * 
 * @author Lars Hoffmann
 */
void fft_help(
  double *fcReal,
  double *fcImag,
  int n);

/**
 * @brief Converts geographic coordinates (longitude, latitude, altitude) to Cartesian coordinates.
 *
 * This function converts geographic coordinates specified by
 * longitude, latitude, and altitude into Cartesian coordinates. The
 * Earth is approximated as a sphere with radius defined by the
 * constant `RE`.
 *
 * @param z The altitude above the Earth's surface in kilometers.
 * @param lon The longitude in degrees.
 * @param lat The latitude in degrees.
 * @param x Pointer to an array of three doubles where the computed Cartesian coordinates (x, y, z) will be stored.
 *
 * The function computes the Cartesian coordinates using the given altitude, longitude, and latitude.
 * It assumes the Earth is a perfect sphere and uses the following formulas:
 * - \f$ x = (\textrm{radius}) \cos(\textrm{lat in radians}) \cos(\textrm{lon in radians}) \f$
 * - \f$ y = (\textrm{radius}) \cos(\textrm{lat in radians}) \sin(\textrm{lon in radians}) \f$
 * - \f$ z = (\textrm{radius}) \sin(\textrm{lat in radians}) \f$
 *
 * @note The constant `RE` is defined as the Earth's radius in kilometers.
 * @note Longitude and latitude should be in degrees.
 *
 * @see https://en.wikipedia.org/wiki/Geographic_coordinate_conversion
 *
 * @author Lars Hoffmann
 */
void geo2cart(
  const double z,
  const double lon,
  const double lat,
  double *x);

/**
 * @brief Retrieves meteorological data for the specified time.
 *
 * This function retrieves meteorological data for the given time `t`
 * and updates the provided pointers to the met0 and met1 structures
 * accordingly. It handles both the initialization and subsequent
 * updates of the meteorological data based on the direction of time
 * integration.
 *
 * @param ctl Pointer to the control structure containing configuration settings.
 * @param clim Pointer to the climate structure.
 * @param t The current time for which meteorological data is to be retrieved.
 * @param met0 Pointer to the pointer of the first meteorological data structure.
 * @param met1 Pointer to the pointer of the second meteorological data structure.
 *
 * The function performs the following steps:
 * - Initializes meteorological data on the first call or when the simulation restarts.
 * - Reads new meteorological data when advancing forward or backward in time.
 * - Swaps pointers to manage double buffering of the meteorological data.
 * - Performs caching to optimize subsequent data retrieval.
 * - Ensures consistency of the meteorological grids.
 *
 * @note This function utilizes GPU acceleration with OpenACC directives if enabled.
 * @note Ensure that `ctl`, `clim`, `met0`, and `met1` are properly initialized before calling this function.
 *
 * @see get_met_help
 * @see read_met
 * @see SELECT_TIMER
 * @see LOG
 * @see ERRMSG
 * @see WARN
 *
 * @author Lars Hoffmann
 */
void get_met(
  ctl_t * ctl,
  clim_t * clim,
  double t,
  met_t ** met0,
  met_t ** met1);

/**
 * @brief Helper function to generate the filename for meteorological data.
 *
 * This function generates the appropriate filename for the
 * meteorological data file based on the provided time `t`, direction
 * `direct`, and the base filename `metbase`.  The filename is
 * formatted according to the specified meteorological data type.
 *
 * @param ctl Pointer to the control structure containing configuration settings.
 * @param t The time for which the meteorological data filename is to be generated.
 * @param direct The direction of time integration (-1 for backward, 1 for forward).
 * @param metbase The base string for the meteorological data filenames.
 * @param dt_met The time interval between meteorological data files.
 * @param filename The generated filename for the meteorological data file.
 *
 * The function performs the following steps:
 * - Rounds the time to fixed intervals `dt_met` based on the direction.
 * - Decodes the time into year, month, day, hour, minute, and second.
 * - Constructs the filename based on the meteorological data type specified in `ctl`.
 * - Replaces placeholders (YYYY, MM, DD, HH) in the base filename with actual date and time values.
 *
 * @note Ensure that `ctl` and `filename` are properly initialized before calling this function.
 *
 * @see jsec2time
 * @see get_met_replace
 *
 * @author Lars Hoffmann
 */
void get_met_help(
  ctl_t * ctl,
  double t,
  int direct,
  char *metbase,
  double dt_met,
  char *filename);

/**
 * @brief Replaces occurrences of a substring in a string with another substring.
 *
 * This function replaces occurrences of the substring `search` in the
 * string `orig` with the substring `repl`. The replacement is
 * performed in-place.
 *
 * @param orig The original string where replacements are to be made.
 * @param search The substring to be replaced.
 * @param repl The substring to replace occurrences of `search`.
 *
 * The function iterates over the original string `orig` and replaces
 * each occurrence of the substring `search` with the substring
 * `repl`. It performs the replacement operation up to three times to
 * ensure multiple occurrences are replaced.
 *
 * @note We use this function to repace the strings `YYYY`, `MM`, and `DD` by year,
 *       month, and day in filenames.
 * @note Ensure that `orig`, `search`, and `repl` are properly initialized
 *       and have sufficient memory allocated before calling this function.
 *
 * @author Lars Hoffmann
 */
void get_met_replace(
  char *orig,
  char *search,
  char *repl);

/**
 * @brief Calculate tropopause data.
 *
 * This function reads and interpolates various meteorological
 * parameters such as tropopause pressure, temperature, and ozone
 * concentration at specified latitudes and longitudes. The
 * interpolated data is stored in the provided arrays.
 *
 * @param met_tropo An integer specifying the type of meteorological data to use.
 * @param ctl Pointer to a `ctl_t` structure that controls the meteorological data processing.
 * @param clim Pointer to a `clim_t` structure containing climatological data.
 * @param met Pointer to a `met_t` structure containing meteorological data.
 * @param lons Array of longitudes at which to interpolate data. The array should have `nx` elements.
 * @param nx Number of longitude points.
 * @param lats Array of latitudes at which to interpolate data. The array should have `ny` elements.
 * @param ny Number of latitude points.
 * @param pt Pointer to an array where the interpolated pressure values will be stored. The array should have `nx * ny` elements.
 * @param zt Pointer to an array where the interpolated height values will be stored. The array should have `nx * ny` elements.
 * @param tt Pointer to an array where the interpolated temperature values will be stored. The array should have `nx * ny` elements.
 * @param qt Pointer to an array where the interpolated specific humidity values will be stored. The array should have `nx * ny` elements.
 * @param o3t Pointer to an array where the interpolated ozone concentration values will be stored. The array should have `nx * ny` elements.
 * @param ps Pointer to an array where the interpolated surface pressure values will be stored. The array should have `nx * ny` elements.
 * @param zs Pointer to an array where the interpolated surface height values will be stored. The array should have `nx * ny` elements.
 *
 * @pre `lons` must have at least `nx` elements.
 * @pre `lats` must have at least `ny` elements.
 * @pre `pt`, `zt`, `tt`, `qt`, `o3t`, `ps`, and `zs` must have at least `nx * ny` elements.
 *
 * @post The arrays `pt`, `zt`, `tt`, `qt`, `o3t`, `ps`, and `zs` will contain the interpolated meteorological data.
 *
 * @note The function utilizes OpenMP for parallel processing of the interpolation tasks.
 *
 * @note The function uses the auxiliary functions `read_met_tropo`, `intpol_met_space_2d`, and `intpol_met_space_3d` for reading and interpolating the tropopause data.
 *
 * @author Lars Hoffmann
 */
void get_tropo(
  int met_tropo,
  ctl_t * ctl,
  clim_t * clim,
  met_t * met,
  double *lons,
  int nx,
  double *lats,
  int ny,
  double *pt,
  double *zt,
  double *tt,
  double *qt,
  double *o3t,
  double *ps,
  double *zs);

/**
 * @brief Interpolates meteorological variables to a given position and time.
 *
 * This function interpolates meteorological variables to a specified
 * position and time. It calculates the interpolated value based on
 * the values provided at two time steps and performs interpolation in
 * time, longitude, latitude, and altitude dimensions.
 *
 * @param met0 Pointer to the meteorological data at the first time step.
 * @param height0 Array containing heights at the first time step.
 * @param array0 Array containing meteorological variable values at the first time step.
 * @param met1 Pointer to the meteorological data at the second time step.
 * @param height1 Array containing heights at the second time step.
 * @param array1 Array containing meteorological variable values at the second time step.
 * @param ts Interpolation time (fractional time between met0 and met1).
 * @param height Altitude at which to interpolate.
 * @param lon Longitude at which to interpolate.
 * @param lat Latitude at which to interpolate.
 * @param var Pointer to store the interpolated variable value.
 * @param ci Array to store the calculated indices.
 * @param cw Array to store the weighting factors.
 * @param init Flag indicating if it's the first call (1) or not (0).
 *
 * The function first restricts the longitude within the range [0,
 * 360) degrees.  It then calculates the horizontal indices (`ci[0]`
 * and `ci[1]`) based on the provided longitude and latitude. Next, it
 * locates the vertical indices for each edge of the column based on
 * the provided height.
 *
 * The function then calculates the weighting factors for time,
 * longitude, latitude, and altitude. It iterates over the
 * interpolation process to determine the altitude weighting
 * factor. After initializing the interpolation parameters, it
 * calculates the interpolated variable value and stores it in the
 * memory location pointed to by `var`.
 *
 * @note Ensure that all arrays (`height0`, `array0`, `height1`, `array1`, `ci`, `cw`)
 *       have sufficient memory allocated before calling this function.
 *
 * @author Jan Clemens
 */
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

/**
 * @brief Interpolates meteorological variables in 3D space.
 *
 * This function interpolates meteorological variables at a specified
 * pressure level and geographic position. It calculates the
 * interpolated value based on the values provided at neighboring grid
 * points and performs interpolation in pressure, longitude, and
 * latitude dimensions.
 *
 * @param met Pointer to the meteorological data.
 * @param array Array containing meteorological variable values.
 * @param p Pressure level at which to interpolate.
 * @param lon Longitude at which to interpolate.
 * @param lat Latitude at which to interpolate.
 * @param var Pointer to store the interpolated variable value.
 * @param ci Array to store the calculated indices.
 * @param cw Array to store the weighting factors.
 * @param init Flag indicating if it's the first call (1) or not (0).
 *
 * The function first checks the longitude and adjusts it if necessary
 * to ensure it falls within the valid range. It then calculates the
 * interpolation indices based on the provided pressure level,
 * longitude, and latitude. Next, it computes the interpolation
 * weights for pressure, longitude, and latitude.
 *
 * The function interpolates vertically first and then
 * horizontally. The interpolated value is stored in the memory
 * location pointed to by `var`.
 *
 * @note Ensure that the `array`, `ci`, and `cw` arrays have sufficient memory allocated
 *       before calling this function.
 *
 * @author Lars Hoffmann
 */
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

/**
 * @brief Interpolates a meteorological variable in 3D space (longitude, latitude, pressure).
 *
 * This function performs trilinear interpolation of a meteorological variable
 * based on the provided longitude, latitude, and pressure coordinates. The
 * meteorological data is given in a 3D array, and the function calculates the
 * interpolated value and stores it in the variable pointed to by `var`.
 * The function operates on model level data.
 *
 * @param[in]  met   Pointer to a `met_t` structure containing the meteorological data.
 * @param[in]  array 3D array of meteorological data with dimensions [EX][EY][EP].
 * @param[in]  p     Pressure coordinate at which to interpolate.
 * @param[in]  lon   Longitude coordinate at which to interpolate.
 * @param[in]  lat   Latitude coordinate at which to interpolate.
 * @param[out] var   Pointer to a double where the interpolated value will be stored.
 *
 * @author Lars Hoffmann
 */
void intpol_met_space_3d_ml(
  met_t * met,
  float array[EX][EY][EP],
  double p,
  double lon,
  double lat,
  double *var);

/**
 * @brief Interpolates meteorological variables in 2D space.
 *
 * This function interpolates meteorological variables at a specified
 * geographic position. It calculates the interpolated value based on
 * the values provided at neighboring grid points and performs
 * interpolation in longitude and latitude dimensions.
 *
 * @param met Pointer to the meteorological data.
 * @param array Array containing meteorological variable values.
 * @param lon Longitude at which to interpolate.
 * @param lat Latitude at which to interpolate.
 * @param var Pointer to store the interpolated variable value.
 * @param ci Array to store the calculated indices.
 * @param cw Array to store the weighting factors.
 * @param init Flag indicating if it's the first call (1) or not (0).
 *
 * The function first checks the longitude and adjusts it if necessary
 * to ensure it falls within the valid range. It then calculates the
 * interpolation indices based on the provided longitude and
 * latitude. Next, it computes the interpolation weights for longitude
 * and latitude.
 *
 * The function interpolates horizontally and stores the interpolated
 * value in the memory location pointed to by `var`. If any of the
 * data values used in interpolation are not finite, the function
 * handles this situation by choosing a valid value or performing a
 * simple interpolation.
 *
 * @note Ensure that the `array`, `ci`, and `cw` arrays have sufficient memory allocated
 *       before calling this function.
 *
 * @author Lars Hoffmann
 */
void intpol_met_space_2d(
  met_t * met,
  float array[EX][EY],
  double lon,
  double lat,
  double *var,
  int *ci,
  double *cw,
  int init);

/**
 * @brief Interpolates meteorological data in 3D space and time.
 *
 * This function interpolates meteorological data in three dimensions
 * (longitude, latitude, and pressure) and time. It calculates the
 * interpolated value based on the values provided at neighboring grid
 * points and performs interpolation both spatially and temporally.
 *
 * @param met0 Pointer to the meteorological data at time t0.
 * @param array0 3D array of meteorological data at time t0.
 * @param met1 Pointer to the meteorological data at time t1.
 * @param array1 3D array of meteorological data at time t1.
 * @param ts Time stamp at which to interpolate.
 * @param p Pressure level at which to interpolate.
 * @param lon Longitude at which to interpolate.
 * @param lat Latitude at which to interpolate.
 * @param var Pointer to store the interpolated value.
 * @param ci Array to store the calculated indices.
 * @param cw Array to store the weighting factors.
 * @param init Flag indicating if it's the first call (1) or not (0).
 *
 * The function first performs spatial interpolation for both time
 * instances (t0 and t1) using the `intpol_met_space_3d` function. It
 * then calculates the weighting factor `wt` based on the time stamp
 * `ts`. Finally, it performs temporal interpolation using the
 * interpolated values at t0 and t1 along with the weighting factor to
 * compute the final interpolated value stored in `var`.
 *
 * @note Ensure that the `ci` and `cw` arrays have sufficient memory allocated
 *       before calling this function.
 *
 * @author Lars Hoffmann
 */
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

/**
 * @brief Interpolates a meteorological variable in time and 3D space (longitude, latitude, pressure).
 *
 * This function performs spatiotemporal interpolation of a meteorological
 * variable based on the provided longitude, latitude, pressure, and timestamp. 
 * The meteorological data is given in two 3D arrays corresponding to two 
 * different time steps, and the function calculates the interpolated value 
 * and stores it in the variable pointed to by `var`.
 * The function operates on model level data.
 *
 * @param[in]  met0   Pointer to a `met_t` structure containing the meteorological data for the first time step.
 * @param[in]  array0 3D array of meteorological data for the first time step with dimensions [EX][EY][EP].
 * @param[in]  met1   Pointer to a `met_t` structure containing the meteorological data for the second time step.
 * @param[in]  array1 3D array of meteorological data for the second time step with dimensions [EX][EY][EP].
 * @param[in]  ts     Timestamp at which to interpolate.
 * @param[in]  p      Pressure coordinate at which to interpolate.
 * @param[in]  lon    Longitude coordinate at which to interpolate.
 * @param[in]  lat    Latitude coordinate at which to interpolate.
 * @param[out] var    Pointer to a double where the interpolated value will be stored.
 *
 * @author Lars Hoffmann
 */
void intpol_met_time_3d_ml(
  met_t * met0,
  float array0[EX][EY][EP],
  met_t * met1,
  float array1[EX][EY][EP],
  double ts,
  double p,
  double lon,
  double lat,
  double *var);

/**
 * @brief Interpolates meteorological data in 2D space and time.
 *
 * This function interpolates meteorological data in two dimensions
 * (longitude and latitude) and time. It calculates the interpolated
 * value based on the values provided at neighboring grid points and
 * performs interpolation both spatially and temporally.
 *
 * @param met0 Pointer to the meteorological data at time t0.
 * @param array0 2D array of meteorological data at time t0.
 * @param met1 Pointer to the meteorological data at time t1.
 * @param array1 2D array of meteorological data at time t1.
 * @param ts Time stamp at which to interpolate.
 * @param lon Longitude at which to interpolate.
 * @param lat Latitude at which to interpolate.
 * @param var Pointer to store the interpolated value.
 * @param ci Array to store the calculated indices.
 * @param cw Array to store the weighting factors.
 * @param init Flag indicating if it's the first call (1) or not (0).
 *
 * The function first performs spatial interpolation for both time
 * instances (t0 and t1) using the `intpol_met_space_2d` function. It
 * then calculates the weighting factor `wt` based on the time stamp
 * `ts`. Finally, it performs temporal interpolation using the
 * interpolated values at t0 and t1 along with the weighting factor to
 * compute the final interpolated value stored in `var`.  If one of
 * the interpolated values is not finite, it selects the valid value
 * based on the weighting factor `wt`.
 *
 * @note Ensure that the `ci` and `cw` arrays have sufficient memory allocated
 *       before calling this function.
 *
 * @author Lars Hoffmann
 */
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

/**
 * @brief Interpolates tropopause data in 3D (latitude, longitude, and time).
 *
 * This function performs interpolation of tropopause data at a given
 * latitude, longitude, and time.  The interpolation can be performed
 * using either linear interpolation or nearest neighbor
 * interpolation.  The standard deviation of the data points used in
 * the interpolation is also computed.
 *
 * @param time0 Time corresponding to the first data array `array0`.
 * @param array0 A 2D array of tropopause data at `time0`. The dimensions are `EX` by `EY`.
 * @param time1 Time corresponding to the second data array `array1`.
 * @param array1 A 2D array of tropopause data at `time1`. The dimensions are `EX` by `EY`.
 * @param lons Array of longitudes with `EX` elements.
 * @param lats Array of latitudes with `EY` elements.
 * @param nlon Number of longitudes.
 * @param nlat Number of latitudes.
 * @param time The specific time at which to interpolate the data.
 * @param lon The specific longitude at which to interpolate the data.
 * @param lat The specific latitude at which to interpolate the data.
 * @param method Interpolation method: `1` for linear interpolation, otherwise nearest neighbor interpolation is used.
 * @param var Pointer to the variable where the interpolated value will be stored.
 * @param sigma Pointer to the variable where the standard deviation of the data points will be stored.
 *
 * @pre `array0` and `array1` must be 2D arrays of size `EX` by `EY`.
 * @pre `lons` must have at least `nlon` elements and `lats` must have at least `nlat` elements.
 *
 * @post `var` will contain the interpolated value.
 * @post `sigma` will contain the standard deviation of the data points used in the interpolation.
 *
 * @note The function adjusts the longitude to ensure it is within the range defined by `lons`.
 * @note This function uses the auxiliary functions `locate_reg`, `LIN`, and `NN` for locating indices and performing interpolation.
 *
 * @warning Ensure that `EX` and `EY` are defined appropriately to match the dimensions of `array0` and `array1`.
 *
 * @author Lars Hoffmann
 */
void intpol_tropo_3d(
  double time0,
  float array0[EX][EY],
  double time1,
  float array1[EX][EY],
  double lons[EX],
  double lats[EY],
  int nlon,
  int nlat,
  double time,
  double lon,
  double lat,
  int method,
  double *var,
  double *sigma);

/**
 * @brief Converts Julian seconds to calendar date and time components.
 *
 * This function converts Julian seconds to calendar date and time
 * components, including year, month, day, hour, minute, and
 * second. It also calculates the fractional part of the seconds.
 *
 * @param jsec Julian seconds to convert.
 * @param year Pointer to store the year.
 * @param mon Pointer to store the month.
 * @param day Pointer to store the day.
 * @param hour Pointer to store the hour.
 * @param min Pointer to store the minute.
 * @param sec Pointer to store the second.
 * @param remain Pointer to store the fractional part of seconds.
 *
 * The function initializes a time structure `t0` with a fixed
 * starting date and time. It then converts the Julian seconds to a
 * time_t type by adding the seconds to the epoch time. Next, it
 * converts the time_t value to a UTC time structure `t1`. Finally, it
 * extracts the year, month, day, hour, minute, and second components
 * from `t1` and calculates the fractional part of seconds, which is
 * stored in `remain`.
 *
 * @author Lars Hoffmann
 */
void jsec2time(
  const double jsec,
  int *year,
  int *mon,
  int *day,
  int *hour,
  int *min,
  int *sec,
  double *remain);

/**
 * @brief Calculates the kernel weight based on altitude and given kernel data.
 *
 * This function calculates the kernel weight based on altitude and
 * given kernel data. It takes arrays of altitudes (`kz`) and
 * corresponding weights (`kw`), the number of data points (`nk`), and
 * the current altitude (`p`) as input.
 *
 * @param kz Array of altitudes.
 * @param kw Array of corresponding weights.
 * @param nk Number of data points.
 * @param p Current altitude.
 * @return The calculated kernel weight.
 *
 * If the number of data points is less than 2 (`nk < 2`), the
 * function returns a default weight of 1.0.
 *
 * The function first computes the altitude `z` based on the current
 * altitude `p`.  Then it checks whether `z` is outside the range of
 * altitudes in the kernel data.  If so, it returns the corresponding
 * weight at the nearest altitude boundary.  Otherwise, it
 * interpolates linearly between the two closest altitudes in the
 * kernel data to determine the weight at altitude `z`.
 *
 * @author Lars Hoffmann
 */
double kernel_weight(
  const double kz[EP],
  const double kw[EP],
  const int nk,
  const double p);

/**
 * @brief Calculates the moist adiabatic lapse rate in Kelvin per kilometer.
 *
 * This function calculates the moist adiabatic lapse rate in Kelvin
 * per kilometer from the given temperature (`t`) in Kelvin and water
 * vapor volume mixing ratio (`h2o`).
 *
 * @param t Temperature in Kelvin.
 * @param h2o Water vapor volume mixing ratio.
 * @return The moist adiabatic lapse rate in Kelvin per kilometer.
 *
 * The moist adiabatic lapse rate is calculated using the formula:
 *
 * \f[ \Gamma = \frac{{1000 \times g \times \left(a + L_v \times r \times T\right)}} {{C_{pd} \times a + L_v^2 \times r \times \epsilon}} \f]
 *
 * where:
 * - \f$ \Gamma \f$ is the lapse rate in Kelvin per kilometer.
 * - \f$ g \f$ is the acceleration due to gravity (constant).
 * - \f$ a = R_a \times T^2 \f$ is a term based on the gas constant for dry air and temperature squared.
 * - \f$ R_a \f$ is the gas constant for dry air.
 * - \f$ T \f$ is the temperature in Kelvin.
 * - \f$ L_v \f$ is the latent heat of vaporization.
 * - \f$ r = \frac{{S_h(h_2o)}}{{1 - S_h(h_2o)}} \f$ is a term based on the water vapor mixing ratio.
 * - \f$ S_h(h_2o) \f$ is the saturation vapor pressure relative to the pressure at saturation.
 * - \f$ C_{pd} \f$ is the specific heat of dry air at constant pressure.
 * - \f$ \epsilon \f$ is the ratio of the gas constants for dry air and water vapor.
 *
 * The constants used in the calculation are defined externally:
 * - \f$ g \f$: Acceleration due to gravity (constant).
 * - \f$ R_a \f$: Gas constant for dry air.
 * - \f$ L_v \f$: Latent heat of vaporization.
 * - \f$ C_{pd} \f$: Specific heat of dry air at constant pressure.
 * - \f$ \epsilon \f$: Ratio of the gas constants for dry air and water vapor.
 *
 * @see [Wikipedia - Lapse rate](https://en.wikipedia.org/wiki/Lapse_rate)
 *
 * @author Lars Hoffmann
 */
double lapse_rate(
  const double t,
  const double h2o);

/**
 * @brief Defines pressure levels for meteorological data.
 *
 * This function defines pressure levels for meteorological data based
 * on the given control structure (`ctl`).  Pressure levels are
 * defined differently based on the value of `press_level_def` in
 * `ctl`.
 *
 * @param ctl Control structure containing information about pressure level definitions.
 *
 * The function determines the number of pressure levels (`met_np`)
 * and the corresponding pressure values (`met_p`) based on the value
 * of `press_level_def` in the control structure `ctl`. It initializes
 * the `met_np` and `met_p` fields accordingly.
 *
 * @note Valid values for `press_level_def` are:
 * - 0: Define 138 pressure levels.
 * - 1: Define 92 pressure levels.
 * - 2: Define 60 pressure levels.
 * - 3: Define 147 pressure levels.
 * - 4: Define 101 pressure levels.
 * - 5: Define 62 pressure levels.
 * - 6: Define 137 pressure levels.
 * - 7: Define 59 pressure levels.
 * Any other value for `press_level_def` will result in an error message.
 *
 * @author Jan Clemens
 */
void level_definitions(
  ctl_t * ctl);

/**
 * @brief Locate the index of the interval containing a given value in a sorted array.
 *
 * This function locates the index of the interval containing a given
 * value in a sorted array.  It uses a binary search algorithm to
 * efficiently find the interval.
 *
 * @param xx Pointer to the sorted array.
 * @param n Size of the array.
 * @param x Value to be located.
 * @return Index of the interval containing the value `x`.
 *
 * The function assumes that the array `xx` is sorted in ascending
 * order.  It returns the index of the interval where the value `x` is
 * located.  If the value `x` is outside the range of the array, the
 * function returns the index of the closest interval.
 *
 * @author Lars Hoffmann
 */
int locate_irr(
  const double *xx,
  const int n,
  const double x);

/**
 * @brief Locate the index of the interval containing a given value in an irregularly spaced array.
 *
 * This function performs a binary search to locate the interval in the array `xx` such that 
 * `xx[ig] <= x < xx[ig + 1]`. If the value `x` lies within the interval specified by the initial 
 * guess index `ig`, the function returns `ig`. Otherwise, it searches the array to find the correct interval.
 *
 * @param xx Pointer to the array of floats representing the irregularly spaced intervals. 
 *           The array must be of size `n`.
 * @param n The number of elements in the array `xx`.
 * @param x The value to locate within the intervals of the array `xx`.
 * @param ig The initial guess index. If the interval `[xx[ig], xx[ig+1])` contains `x`, 
 *           the function returns `ig` directly.
 * 
 * @return The index `i` such that `xx[i] <= x < xx[i + 1]`. If `x` is out of bounds, 
 *         it returns the index of the closest interval.
 *
 * @note The function assumes that the array `xx` contains at least two elements.
 * @note The function can handle both increasing and decreasing sequences in the array `xx`.
 *
 * @warning The behavior is undefined if the array `xx` is not sorted in either increasing or 
 *          decreasing order, or if it contains less than two elements.
 *
 * @author Lars Hoffmann
 */
int locate_irr_float(
  const float *xx,
  const int n,
  const double x,
  const int ig);

/**
 * @brief Locate the index of the interval containing a given value in a 3D irregular grid.
 *
 * This function locates the index of the interval containing a given
 * value in a 3D irregular grid.  It searches for the interval in the
 * specified longitude and latitude indices of the grid.
 *
 * @param profiles 3D array representing the irregular grid.
 * @param np Size of the profile (number of data points).
 * @param ind_lon Index of the longitude.
 * @param ind_lat Index of the latitude.
 * @param x Value to be located.
 * @return Index of the interval containing the value `x`.
 *
 * The function assumes that the array `profiles` represents a 3D
 * irregular grid.  It calculates the index of the interval where the
 * value `x` is located based on the data in the specified longitude
 * and latitude indices.  If the value `x` is outside the range of the
 * profile, the function returns the index of the closest interval.
 *
 * @author Jan Clemens
 */
int locate_irr_3d(
  float profiles[EX][EY][EP],
  int np,
  int ind_lon,
  int ind_lat,
  double x);

/**
 * @brief Locate the index of the interval containing a given value in a regular grid.
 *
 * This function locates the index of the interval containing a given
 * value in a regular grid.  It calculates the index based on the
 * spacing between grid points and the value to be located.
 *
 * @param xx Pointer to the array representing the regular grid.
 * @param n Size of the grid (number of grid points).
 * @param x Value to be located.
 * @return Index of the interval containing the value `x`.
 *
 * The function assumes that the array `xx` represents a regular grid
 * with equally spaced points.  It calculates the index of the
 * interval where the value `x` is located based on the spacing
 * between grid points.  If the value `x` is outside the range of the
 * grid, the function returns the index of the closest interval.
 *
 * @author Lars Hoffmann
 */
int locate_reg(
  const double *xx,
  const int n,
  const double x);

/**
 * @brief Locate the four vertical indizes of a box for a given height value.
 *
 * This function locates the vertical indices corresponding to a given
 * height in a 3D irregular grid.  It calculates the indices based on
 * the specified longitude and latitude indices of the grid.
 *
 * @param profiles 3D array representing the irregular grid.
 * @param np Size of the profile (number of data points).
 * @param lon_ap_ind Index of the longitude.
 * @param lat_ap_ind Index of the latitude.
 * @param alt_ap Height value.
 * @param ind Pointer to an array to store the resulting indices.
 *
 * The function calculates the indices corresponding to the specified
 * height in the 3D irregular grid. It stores the resulting indices
 * in the array pointed to by `ind`. The indices are calculated based
 * on the specified longitude and latitude indices of the grid.
 *
 * @author Lars Hoffmann
 */
void locate_vert(
  float profiles[EX][EY][EP],
  int np,
  int lon_ap_ind,
  int lat_ap_ind,
  double alt_ap,
  int *ind);

/**
 * @brief Performs the advection of atmospheric particles using meteorological data.
 * 
 * This function advects particles in the atmosphere using
 * meteorological data from two time steps, updating their positions
 * and interpolating necessary data.  It supports both pressure and
 * zeta vertical coordinate systems.
 * 
 * @param ctl   Pointer to the control structure containing configuration flags.
 * @param met0  Pointer to the initial meteorological data structure.
 * @param met1  Pointer to the final meteorological data structure.
 * @param atm   Pointer to the air parcel data structure.
 * @param dt    Array of time step values for each particle.
 * 
 * @details
 * The function performs the following operations:
 * - Sets up a timer labeled "MODULE_ADVECTION" within the "PHYSICS"
     category for GPU profiling using NVTX.
 * - Depending on the vertical coordinate system (pressure or zeta),
 *   it loops over each particle in the atmosphere (atm->np),
 *   initializing and updating their positions and velocities using
 *   meteorological data interpolation.
 * 
 * ### Pressure Coordinate System (ctl->vert_coord_ap == 0)
 * - Initializes particle positions and velocities.
 * - Performs integration over a specified number of nodes (ctl->advect) to update particle positions.
 * - Interpolates meteorological data for velocity components (u, v, w).
 * - Computes mean velocities and updates particle positions in longitude, latitude, and pressure.
 * 
 * ### Zeta Coordinate System (ctl->vert_coord_ap == 1)
 * - Translates pressure into zeta coordinate if other modules have modified the pressure.
 * - Initializes particle positions and velocities in zeta coordinates.
 * - Performs integration over a specified number of nodes (ctl->advect) to update particle positions.
 * - Interpolates meteorological data for velocity components (u, v) and zeta_dot.
 * - Computes mean velocities and updates particle positions in longitude, latitude, and zeta.
 * - Checks and corrects if zeta is below zero.
 * - Translates updated zeta back into pressure coordinates.
 * 
 * @note The function assumes that the atmospheric data structure
 * (atm) contains arrays for time, longitude, latitude, and pressure
 * for each point.  The specific quantification of zeta is stored in
 * atm->q[ctl.qnt_zeta].
 *
 * @authors Lars Hoffmann
 * @authors Jan Clemens
 */
void module_advect(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/**
 * @brief Initializes the advection module by setting up pressure fields.
 * 
 * This function initializes the advection module, setting up the air parcel pressure
 * to be consistent with the given zeta vertical coordinate.
 * It utilizes meteorological data from two time steps and interpolates
 * the pressure values accordingly.
 * 
 * @param ctl   Pointer to the control structure containing configuration flags.
 * @param met0  Pointer to the initial meteorological data structure.
 * @param met1  Pointer to the final meteorological data structure.
 * @param atm   Pointer to the air parcel data structure.
 * 
 * @details
 * The function performs the following operations:
 * - Sets up a timer labeled "MODULE_ADVECTION" within the "PHYSICS" category.
 * - If the zeta vertical coordinate system is specified (ctl->vert_coord_ap == 1), it initializes
 *   the pressure fields to be consistent with the zeta coordinate using 4D interpolation.
 * 
 * @author Jan Clemens
 */
void module_advect_init(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm);

/**
 * @brief Apply boundary conditions to particles based on meteorological and climatological data.
 *
 * This function applies boundary conditions to particles based on
 * specified criteria, including latitude, pressure, surface layer
 * parameters, and climatological data. It loops over each particle
 * and checks whether it satisfies the specified boundary
 * conditions. If a particle satisfies the conditions, its properties
 * such as mass, volume mixing ratio, and age of air are updated
 * accordingly.
 *
 * It checks for quantity flags to determine which properties need to
 * be updated. If the latitude or pressure of a particle falls outside
 * the specified ranges, it skips the particle. It also considers
 * surface layer parameters such as surface pressure, height, zeta
 * range, and planetary boundary layer. If a particle is within the
 * specified surface layer boundaries, its properties are updated
 * accordingly.
 *
 * The function updates properties such as mass and volume mixing
 * ratio if the corresponding flags are set. It retrieves volume
 * mixing ratio values for various trace gases (e.g., CFC-10, CFC-11,
 * N2O, SF6) from climatological time series data and updates the
 * particle properties accordingly.  Additionally, it updates the age
 * of air for each particle based on the current simulation time.
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param clim Pointer to the climatological data structure containing time series data.
 * @param met0 Pointer to the meteorological data structure at the initial time step.
 * @param met1 Pointer to the meteorological data structure at the next time step.
 * @param atm Pointer to the atmospheric data structure containing particle information.
 * @param dt Pointer to the time step value.
 *
 * @authors Lars Hoffmann
 * @authors Mingzhao Liu
 */
void module_bound_cond(
  ctl_t * ctl,
  clim_t * clim,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/**
 * @brief Calculate grid data for chemistry modules.
 *
 * This function initializes and updates chemical grid quantities
 * based on atmospheric data and interpolation of meteorological
 * variables.
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param clim Pointer to the climate data structure containing climatological data.
 * @param met0 Pointer to the first meteorological data structure.
 * @param met1 Pointer to the second meteorological data structure.
 * @param atm Pointer to the atmospheric data structure containing particle information.
 * @param t Time for which chemical grid is updated.
 *
 * @authors Mingzhao Liu
 * @authors Lars Hoffmann
 */
void module_chemgrid(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double t);

/**
 * @brief Initializes the chemistry modules by setting atmospheric composition.
 * 
 * This function initializes various chemical components of the
 * atmosphere using meteorological data and climatological
 * information. It interpolates and sets values for water vapor (H2O),
 * ozone (O3), and several radical species such as OH, HO2, H2O2, and
 * O1D for each air parcel.
 * 
 * @param ctl   Pointer to the control structure containing quantity flags.
 * @param clim  Pointer to the climatology structure containing climatological data.
 * @param met0  Pointer to the initial meteorological data structure.
 * @param met1  Pointer to the final meteorological data structure.
 * @param atm   Pointer to the air parcel data structure.
 * 
 * @details
 * The function uses OpenMP for parallel processing, iterating over
 * each point in the atmosphere (atm->np) to initialize chemical
 * species concentrations.  It performs the following steps:
 * - Interpolates H2O and O3 data from meteorological input if the
 *   respective quantity flags (ctl->qnt_Ch2o and ctl->qnt_Co3) are set.
 * - Sets the concentrations of OH, HO2, H2O2, and O1D using
 *   climatological data if the respective quantity flags are set.
 *
 * @authors Mingzhao Liu
 */
void module_chem_init(
  ctl_t * ctl,
  clim_t * clim,
  met_t * met0,
  met_t * met1,
  atm_t * atm);

/**
 * @brief Simulate convective processes for atmospheric particles.
 *
 * This function simulates convective processes for atmospheric
 * particles based on convective available potential energy (CAPE) and
 * convective inhibition (CIN).  It loops over each particle and
 * checks whether the CAPE exceeds a specified threshold. If CAPE is
 * above the threshold and meets the conditions, it performs
 * convective mixing to update the particle's pressure based on random
 * numbers generated for each particle.
 *
 * The function interpolates CAPE and CIN from meteorological data and
 * checks whether the particle is above the cloud top level. If the
 * conditions are met, it determines the pressure range for vertical
 * mixing and calculates the density range based on temperature at
 * different pressure levels. It then calculates the new density based
 * on random numbers and updates the particle's pressure accordingly.
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param met0 Pointer to the meteorological data structure at the initial time step.
 * @param met1 Pointer to the meteorological data structure at the next time step.
 * @param atm Pointer to the atmospheric data structure containing particle information.
 * @param dt Pointer to the time step value.
 * @param rs Pointer to an array of random numbers generated for each particle.
 *
 * @author Lars Hoffmann
 */
void module_convection(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt,
  double *rs);

/**
 * @brief Simulate exponential decay processes for atmospheric particles.
 *
 * This function simulates decay processes for atmospheric particles
 * based on their mass or volume mixing ratio. It loops over each
 * particle and calculates the decay rate using weighting factors for
 * tropospheric and stratospheric lifetimes.  Exponential decay is
 * then calculated, and the mass or volume mixing ratio of particles
 * is updated accordingly. Loss rates can also be calculated and
 * updated based on the decay process.
 *
 * The function checks for quantity flags to ensure that mass or
 * volume mixing ratio data is available. It then calculates the
 * weighting factor based on the particle's location in the atmosphere
 * and sets the lifetime accordingly. Exponential decay is calculated
 * using the time step and the lifetime, and the particle's mass or
 * volume mixing ratio is updated. Loss rates can also be updated
 * based on the decay process.
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param clim Pointer to the climate data structure containing atmospheric data.
 * @param atm Pointer to the atmospheric data structure containing particle information.
 * @param dt Pointer to the time step value.
 *
 * @author Lars Hoffmann
 */
void module_decay(
  ctl_t * ctl,
  clim_t * clim,
  atm_t * atm,
  double *dt);

/**
 * @brief Simulate mesoscale diffusion for atmospheric particles.
 *
 * This function simulates mesoscale diffusion for atmospheric
 * particles, including horizontal and vertical wind fluctuations. It
 * calculates standard deviations of local wind data and temporal
 * correlations for mesoscale fluctuations. Mesoscale wind
 * fluctuations are then calculated based on the provided random
 * numbers and turbulence parameters. The particle positions are
 * updated accordingly.
 *
 * The function loops over each particle and calculates indices for
 * interpolation of wind data. It then computes standard deviations of
 * local wind data and temporal correlations for mesoscale
 * fluctuations. Based on the turbulence parameters and provided
 * random numbers, it calculates horizontal and vertical mesoscale
 * wind fluctuations. Finally, it updates the particle positions based
 * on the calculated wind fluctuations.
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param met0 Pointer to the meteorological data structure at the current time step.
 * @param met1 Pointer to the meteorological data structure at the next time step.
 * @param atm Pointer to the atmospheric data structure containing particle information.
 * @param cache Pointer to the cache structure for temporary storage.
 * @param dt Pointer to the time step value.
 * @param rs Pointer to the array of random numbers.
 *
 * @author Lars Hoffmann
 */
void module_diffusion_meso(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  cache_t * cache,
  double *dt,
  double *rs);

/**
 * @brief Simulate turbulent diffusion for atmospheric particles.
 *
 * This function simulates turbulent diffusion for atmospheric
 * particles, including horizontal and vertical diffusion. It
 * calculates diffusivity based on the provided weighting factor and
 * turbulence parameters. The diffusion coefficients are then used to
 * calculate horizontal and vertical displacements of particles based
 * on the provided random numbers and time step.
 *
 * The function loops over each particle and calculates the weighting
 * factor based on the atmospheric properties. It then computes
 * diffusivity for horizontal and vertical diffusion. Based on the
 * diffusivity values and provided random numbers, it calculates
 * horizontal and vertical displacements of particles and updates
 * their positions accordingly.
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param clim Pointer to the climatological data structure containing atmospheric properties.
 * @param atm Pointer to the atmospheric data structure containing particle information.
 * @param dt Pointer to the time step value.
 * @param rs Pointer to the array of random numbers.
 *
 * @author Lars Hoffmann
 */
void module_diffusion_turb(
  ctl_t * ctl,
  clim_t * clim,
  atm_t * atm,
  double *dt,
  double *rs);

/**
 * @brief Simulate dry deposition of atmospheric particles.
 *
 * This function simulates the dry deposition of atmospheric
 * particles, including both particulate matter and gases. It
 * calculates the sedimentation velocity for particles based on the
 * atmospheric properties and applies it to determine the loss of mass
 * or volume mixing ratio due to deposition. The function loops over
 * each particle and calculates the loss of mass or volume mixing
 * ratio based on the deposition velocity and time step.
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param met0 Pointer to the meteorological data structure at the current time step.
 * @param met1 Pointer to the meteorological data structure at the next time step.
 * @param atm Pointer to the atmospheric data structure containing particle information.
 * @param dt Pointer to the time step value.
 *
 * @author Lars Hoffmann
 */
void module_dry_deposition(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/**
 * @brief Perform chemical reactions involving H2O2 within cloud particles.
 *
 * This function simulates chemical reactions involving hydrogen
 * peroxide (H2O2) within cloud particles. It calculates the change in
 * H2O2 concentration over time due to chemical reactions. The
 * reaction rates are determined based on temperature and cloud
 * properties such as liquid water content.
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param clim Pointer to the climatological data structure.
 * @param met0 Pointer to the first meteorological data structure.
 * @param met1 Pointer to the second meteorological data structure.
 * @param atm Pointer to the atmospheric data structure containing particle information.
 * @param dt Pointer to an array containing the time step for each particle.
 *
 * @note The function assumes that the necessary control structure (ctl), climatological
 *       data structure (clim), meteorological data structures (met0, met1), and atmospheric
 *       data structure (atm) have been initialized and are accessible.
 * @note Chemical reactions involving H2O2 are simulated for particles within clouds, as
 *       indicated by a positive liquid water content (LWC).
 * @note The function calculates reaction rates based on temperature and cloud properties,
 *       including the liquid water content (LWC) and the concentration of SO2.
 * @note The exponential decay of H2O2 concentration due to chemical reactions is calculated
 *       using the reaction rate coefficient and the time step (dt) for each particle.
 * @note If the particle has a quantity flag for either mass (ctl->qnt_m) or volume mixing
 *       ratio (ctl->qnt_vmr), the function updates the quantity based on the exponential decay.
 * @note If the particle has a loss rate quantity flag (ctl->qnt_loss_rate), the function
 *       accumulates the reaction rate coefficient to quantify the loss rate.
 *
 * @author Mingzhao Liu
 */
void module_h2o2_chem(
  ctl_t * ctl,
  clim_t * clim,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/**
 * @brief Initialize the isosurface module based on atmospheric data.
 *
 * This function initializes the isosurface module based on the
 * atmospheric data provided. It calculates the necessary variables
 * required for generating the isosurface, such as pressure, density,
 * or potential temperature. Additionally, it can read balloon
 * pressure data from a file if specified in the control
 * structure. The initialized data is stored in the cache for later
 * use.
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param met0 Pointer to the meteorological data structure at the current time step.
 * @param met1 Pointer to the meteorological data structure at the next time step.
 * @param atm Pointer to the atmospheric data structure containing particle information.
 * @param cache Pointer to the cache structure for storing initialized data.
 *
 * @author Lars Hoffmann
 */
void module_isosurf_init(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  cache_t * cache);

/**
 * @brief Apply the isosurface module to adjust atmospheric properties.
 *
 * This function applies the isosurface module to adjust atmospheric
 * properties based on the initialized data stored in the cache. It
 * interpolates and restores atmospheric pressure, density, or
 * potential temperature according to the specified isosurface mode in
 * the control structure.
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param met0 Pointer to the meteorological data structure at the current time step.
 * @param met1 Pointer to the meteorological data structure at the next time step.
 * @param atm Pointer to the atmospheric data structure containing particle information.
 * @param cache Pointer to the cache structure containing initialized data.
 * @param dt Array of time step values for each particle.
 *
 * @author Lars Hoffmann
 */
void module_isosurf(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  cache_t * cache,
  double *dt);

/*! KPP chemistry module. */
/**
 * @brief Simulate chemical reactions using the Kinetic PreProcessor (KPP) integration scheme.
 *
 * This function simulates chemical reactions using the Kinetic
 * PreProcessor (KPP) integration scheme for atmospheric particles. It
 * loops over each particle in the atmospheric data structure and
 * integrates chemical reactions over a specified time step using the
 * KPP algorithm.
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param clim Pointer to the climatological data structure.
 * @param met0 Pointer to the first meteorological data structure.
 * @param met1 Pointer to the second meteorological data structure.
 * @param atm Pointer to the atmospheric data structure containing particle information.
 * @param dt Pointer to an array containing the time step for each particle.
 *
 * @note The function initializes a timer to measure the execution time of the chemical simulation.
 * @note Chemical integration using KPP is performed for particles with a positive time step (dt > 0).
 * @note For each particle, the function allocates memory for variable (VAR) and fixed (FIX) arrays,
 *       sets the range of time steps (STEPMIN and STEPMAX), and defines relative and absolute tolerances.
 * @note The chemical system is initialized for each particle using the kpp_chem_initialize function.
 * @note Chemical integration is performed over a specified time step (ctl->dt_kpp) using the INTEGRATE
 *       macro, which is part of the KPP integration scheme.
 * @note The function outputs the integrated chemical concentrations back to the atmospheric data structure
 *       using the kpp_chem_output2atm function.
 * @note Memory allocated for the variable (VAR) and fixed (FIX) arrays is freed after the integration is
 *       completed for each particle.
 *
 * @author Mingzhao Liu
 */
void module_kpp_chem(
  ctl_t * ctl,
  clim_t * clim,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/**
 * @brief Update atmospheric properties using meteorological data.
 *
 * This function updates atmospheric properties based on
 * meteorological data interpolated between two time steps. It
 * calculates various atmospheric quantities such as pressure,
 * temperature, wind speed, humidity, etc., and updates the
 * corresponding fields in the atmospheric data structure.
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param clim Pointer to the climate data structure containing climatological data.
 * @param met0 Pointer to the meteorological data structure at the current time step.
 * @param met1 Pointer to the meteorological data structure at the next time step.
 * @param atm Pointer to the atmospheric data structure containing particle information.
 * @param dt Array of time step values for each particle.
 *
 * @author Lars Hoffmann
 */
void module_meteo(
  ctl_t * ctl,
  clim_t * clim,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/**
 * @brief Update atmospheric properties through interparcel mixing.
 *
 * This function updates atmospheric properties by performing
 * interparcel mixing based on the given meteorological and
 * climatological data. It calculates the indices of grid boxes and
 * performs mixing for various quantities such as mass, volume mixing
 * ratio, and other chemical species concentrations.
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param clim Pointer to the climate data structure containing climatological data.
 * @param atm Pointer to the atmospheric data structure containing particle information.
 * @param t Time at which mixing is performed.
 *
 * @authors Mingzhao Liu
 * @authors Lars Hoffmann
 */
void module_mixing(
  ctl_t * ctl,
  clim_t * clim,
  atm_t * atm,
  double t);

/**
 * @brief Perform interparcel mixing for a specific quantity.
 *
 * This function performs interparcel mixing for a specific quantity
 * based on the given indices of grid boxes. It calculates the mean
 * concentration within each grid box and adjusts the quantity for
 * each particle accordingly.
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param clim Pointer to the climate data structure containing climatological data.
 * @param atm Pointer to the atmospheric data structure containing particle information.
 * @param ixs Pointer to the array of grid box indices along the longitude direction.
 * @param iys Pointer to the array of grid box indices along the latitude direction.
 * @param izs Pointer to the array of grid box indices along the vertical direction.
 * @param qnt_idx Index of the quantity for which mixing is performed.
 *
 * @authors Mingzhao Liu
 * @authors Lars Hoffmann
 */
void module_mixing_help(
  ctl_t * ctl,
  clim_t * clim,
  atm_t * atm,
  const int *ixs,
  const int *iys,
  const int *izs,
  int qnt_idx);

/**
 * @brief Perform hydroxyl chemistry calculations for atmospheric particles.
 *
 * This function calculates the OH chemistry for each atmospheric
 * particle based on the specified reaction mechanism and updates the
 * particle quantities accordingly.  The OH chemistry includes
 * bimolecular and termolecular reactions, and the reaction rates are
 * calculated based on the provided climatological data and
 * atmospheric conditions. The function supports both mass and volume
 * mixing ratio quantities for the particles.
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param clim Pointer to the climate data structure containing climatological data.
 * @param met0 Pointer to the first meteorological data structure.
 * @param met1 Pointer to the second meteorological data structure.
 * @param atm Pointer to the atmospheric data structure containing particle information.
 * @param dt Array of time steps for each particle.
 *
 * @note The function assumes that the necessary meteorological and climatological
 *       data structures have been initialized and are accessible via the pointers
 *       met0, met1, and clim, respectively.
 * @note The reaction rates are calculated based on the provided reaction mechanism
 *       and atmospheric conditions, including temperature, pressure, and the
 *       concentrations of relevant species.
 * @note The function updates the particle quantities based on the calculated
 *       reaction rates and the specified time steps. The update can include both
 *       mass and volume mixing ratio quantities, as determined by the control
 *       structure (ctl).
 *
 * @authors Lars Hoffmann
 * @authors Mingzhao Liu
 */
void module_oh_chem(
  ctl_t * ctl,
  clim_t * clim,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/**
 * @brief Update the positions and pressure levels of atmospheric particles.
 *
 * This function updates the positions and pressure levels of atmospheric particles based on the
 * meteorological data and the specified time step. It loops over each particle in the atmospheric
 * data structure and performs the following operations:
 * - Initializes variables required for interpolation.
 * - Calculates modulo for longitude and latitude to ensure they remain within valid ranges.
 * - Adjusts latitude if it exceeds the range [-90, 90] degrees.
 * - Adjusts longitude if it exceeds the range [-180, 180] degrees.
 * - Checks and adjusts pressure levels:
 *   - Reflects pressure levels if they are below the minimum pressure in meteorological data.
 *   - Clamps pressure levels to the maximum pressure in meteorological data if they exceed a
 *     predefined threshold (300 hPa).
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param met0 Pointer to the first meteorological data structure.
 * @param met1 Pointer to the second meteorological data structure.
 * @param atm Pointer to the atmospheric data structure containing particle information.
 * @param dt Pointer to an array containing the time step for each particle.
 *
 * @note The function initializes a timer to measure the execution time of the position update process.
 * @note Position and pressure updates are performed for each particle using linear interpolation.
 * @note Longitude and latitude are adjusted to ensure they remain within valid ranges.
 * @note Pressure levels are adjusted based on meteorological data and a predefined threshold.
 *
 * @author Lars Hoffmann
 */
void module_position(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/**
 * @brief Initialize random number generators for parallel tasks.
 *
 * This function initializes random number generators for parallel
 * tasks using both GSL (GNU Scientific Library) and cuRAND (NVIDIA
 * CUDA Random Number Generation Library) if available. It sets up GSL
 * random number generators for each OpenMP thread and initializes
 * them with unique seeds. For cuRAND, it creates a pseudo-random
 * number generator and sets its seed. The initialization ensures that
 * each task or thread has its own independent random number generator
 * to prevent interference between parallel executions.
 *
 * @param ntask The number of tasks or parallel threads for which random number generators are initialized.
 *
 * @note This function must be called before using any random number generation functions to ensure proper 
 * initialization of random number generators.
 * @note GSL random number generators are initialized for each OpenMP thread, while cuRAND is initialized 
 * for the entire task set.
 * @note If cuRAND is not available (CURAND macro not defined), the cuRAND initialization section is skipped.
 * @note Random number generators are allocated and seeded uniquely for each task or thread to ensure 
 * independence and avoid interference between parallel executions.
 *
 * @author Lars Hoffmann
 */
void module_rng_init(
  int ntask);

/**
 * @brief Generate random numbers using various methods and distributions.
 *
 * This function generates random numbers using different methods and
 * distributions based on the specified method and random number
 * generator type. It supports uniform and normal distributions and
 * can utilize GSL, Squares (Widynski, 2022), or cuRAND random number
 * generators.
 *
 * @param ctl Pointer to the control structure containing parameters and settings.
 * @param rs Pointer to the array where the generated random numbers will be stored.
 * @param n The number of random numbers to generate.
 * @param method The method for generating random numbers: 
 *               - 0: Uniform distribution
 *               - 1: Normal distribution
 *
 * @note The function selects the appropriate random number generator based on the specified method and 
 * the random number generator type defined in the control structure (ctl->rng_type).
 * @note For uniform distribution, the generated random numbers are in the range [0, 1).
 * @note For normal distribution, the Box-Muller transform is used to generate pairs of random numbers 
 * and transform them into a normal distribution.
 * @note If cuRAND is not available (CURAND macro not defined), the function returns an error message.
 *
 * @author Lars Hoffmann
 */
void module_rng(
  ctl_t * ctl,
  double *rs,
  size_t n,
  int method);

/**
 * @brief Simulate sedimentation of particles in the atmosphere.
 *
 * This function calculates the sedimentation velocity of particles
 * based on atmospheric pressure, temperature, and particle properties
 * such as radius and density. It then updates the pressure of each
 * particle based on the sedimentation velocity and the specified time
 * step.
 *
 * @param ctl Pointer to the control structure containing parameters and settings.
 * @param met0 Pointer to the meteorological data at the current time step.
 * @param met1 Pointer to the meteorological data at the next time step.
 * @param atm Pointer to the atmospheric data containing particle information.
 * @param dt Pointer to the array of time steps for each particle.
 *
 * @note The sedimentation velocity is calculated using the `sedi` function, which takes atmospheric pressure, 
 * temperature, particle radius, and particle density as inputs.
 * @note The pressure change for each particle is calculated based on the sedimentation velocity and the 
 * specified time step using the `DZ2DP` function.
 *
 * @author Lars Hoffmann
 */
void module_sedi(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/*!  */
/**
 * @brief Sort particles according to box index.
 *
 * This function sorts particles within the atmosphere data structure
 * based on their geographical coordinates (longitude and latitude)
 * and pressure level. It allocates temporary arrays to store indices
 * and auxiliary data for sorting, then performs the sorting
 * operation. After sorting, it updates the order of particles in the
 * atmosphere data structure.
 *
 * @param ctl Pointer to the control structure containing parameters and settings.
 * @param met0 Pointer to the meteorological data at the current time step.
 * @param atm Pointer to the atmospheric data containing particle information.
 *
 * @note The function utilizes the `locate_reg` and `locate_irr` functions to determine the 
 * appropriate index for sorting particles based on their longitude, latitude, and pressure level.
 * @note Particle sorting is performed using either the Thrust library (if compiled with Thrust support) 
 * or a custom sorting algorithm. If compiled without Thrust support, an error message is displayed.
 * @note After sorting, the function updates the order of particle-related data arrays in the atmosphere 
 * data structure to maintain consistency.
 *
 * @author Lars Hoffmann
 */
void module_sort(
  ctl_t * ctl,
  met_t * met0,
  atm_t * atm);

/**
 * @brief Reorder an array based on a given permutation.
 *
 * This function reorders the elements of a given array based on a
 * specified permutation array.  It allocates temporary memory to
 * store the reordered elements, performs the reordering operation,
 * and then updates the original array with the reordered elements.
 *
 * @param a Pointer to the array to be reordered.
 * @param p Pointer to the permutation array defining the order of elements.
 * @param np The number of elements in the array.
 *
 * @note The function utilizes temporary memory to store the reordered elements before updating 
 * the original array to prevent data loss or corruption.
 * @note Reordering is performed based on the permutation array `p`, which defines the new order 
 * of elements in the array `a`.
 *
 * @author Lars Hoffmann
 */
void module_sort_help(
  double *a,
  const int *p,
  const int np);

/**
 * @brief Calculate time steps for air parcels based on specified conditions.
 *
 * This function calculates the time steps for air parcels based on
 * specified conditions, including the direction of simulation, start
 * and stop times, and a given target time.  It adjusts the time step
 * for each air parcel accordingly and checks for horizontal boundary
 * conditions of local meteorological data.
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param met0 Pointer to the initial meteorological data structure.
 * @param atm Pointer to the atmospheric data structure containing air parcel information.
 * @param dt Pointer to the array storing the calculated time steps for air parcels.
 * @param t The target time for which time steps are calculated.
 *
 * @note The function sets the time step for each air parcel based on its current time 
 * relative to the start and stop times of the simulation, as well as the specified 
 * target time `t`.
 * @note It also checks for horizontal boundaries of local meteorological data and adjusts 
 * the time step accordingly if necessary.
 *
 * @author Lars Hoffmann
 */
void module_timesteps(
  ctl_t * ctl,
  met_t * met0,
  atm_t * atm,
  double *dt,
  double t);

/**
 * @brief Initialize start time and time interval for time-stepping.
 *
 * This function initializes the start time and time interval for
 * time-stepping based on the direction of simulation and the provided
 * atmospheric data. It sets the start time according to the minimum
 * or maximum time in the atmospheric data, depending on the
 * simulation direction. Additionally, it checks the time interval and
 * adjusts the start time accordingly for rounding purposes.
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param atm Pointer to the atmospheric data structure containing air parcel information.
 *
 * @note The function sets the start time based on the direction of simulation and the 
 * minimum or maximum time in the atmospheric data.
 * @note It checks the time interval to ensure that there is a valid time range for 
 * simulation and adjusts the start time for rounding purposes.
 *
 * @author Lars Hoffmann
 */
void module_timesteps_init(
  ctl_t * ctl,
  atm_t * atm);

/**
 * @brief Simulate chemical reactions involving long-lived atmospheric tracers.
 *
 * This function simulates chemical reactions involving atmospheric
 * tracers, such as CFC-10, CFC-11, CFC-12, and N2O. It calculates the
 * change in tracer concentrations over time based on reaction rates
 * and environmental factors such as temperature, ozone concentration,
 * solar zenith angle, and O(1D) volume mixing ratio.
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param clim Pointer to the climatological data structure.
 * @param met0 Pointer to the first meteorological data structure.
 * @param met1 Pointer to the second meteorological data structure.
 * @param atm Pointer to the atmospheric data structure containing particle information.
 * @param dt Pointer to an array containing the time step for each particle.
 *
 * @note The function assumes that the necessary control structure (ctl), climatological
 *       data structure (clim), meteorological data structures (met0, met1), and atmospheric
 *       data structure (atm) have been initialized and are accessible.
 * @note Chemical reactions involving CFC-10, CFC-11, CFC-12, and N2O are simulated for each
 *       particle in the atmospheric data structure.
 * @note The function calculates reaction rates based on temperature, solar zenith angle,
 *       total column ozone, and the volume mixing ratio of O(1D).
 * @note The exponential decay of tracer concentrations due to chemical reactions is calculated
 *       using reaction rate coefficients and the time step (dt) for each particle.
 * @note If the particle has a quantity flag for the tracer species (e.g., ctl->qnt_Cccl4,
 *       ctl->qnt_Cccl3f, ctl->qnt_Cccl2f2, ctl->qnt_Cn2o), the function updates the
 *       concentration of the tracer based on the exponential decay.
 *
 * @authors Mingzhao Liu
 * @authors Lars Hoffmann
 */
void module_tracer_chem(
  ctl_t * ctl,
  clim_t * clim,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/**
 * @brief Perform wet deposition calculations for air parcels.
 *
 * This function calculates the wet deposition process for each air
 * parcel based on provided atmospheric and meteorological data. It
 * estimates the precipitation rate and scavenging coefficients for
 * particles and gases inside and below cloud layers.  The scavenging
 * coefficients are used to calculate the exponential decay of mass or
 * volume mixing ratio over time due to wet deposition.
 *
 * @param ctl Pointer to the control structure containing simulation parameters.
 * @param met0 Pointer to the initial meteorological data structure.
 * @param met1 Pointer to the updated meteorological data structure.
 * @param atm Pointer to the atmospheric data structure containing air parcel information.
 * @param dt Array containing the time step for each air parcel.
 *
 * @note The function calculates the wet deposition process for particles and gases 
 * based on precipitation rate and scavenging coefficients inside and below cloud layers.
 * @note It estimates the exponential decay of mass or volume mixing ratio over time 
 * due to wet deposition.
 *
 * @authors Lars Hoffmann
 * @authors Mingzhao Liu
 */
void module_wet_deposition(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/**
 * @brief Calculates the nitric acid trihydrate (NAT) temperature.
 *
 * This function computes the temperature at which nitric acid
 * trihydrate (NAT) can form given the partial pressures of water
 * vapor and nitric acid in the atmosphere.
 *
 * @param p     The total atmospheric pressure (in hPa).
 * @param h2o   The volume mixing ratio of water vapor (H2O).
 * @param hno3  The volume mixing ratio of nitric acid (HNO3).
 * @return      The NAT temperature (in Kelvin).
 *
 * This function follows these steps:
 * - Ensures the water vapor volume mixing ratio is above a minimum threshold.
 * - Converts the volume mixing ratios of H2O and HNO3 to partial pressures.
 * - Uses these partial pressures to compute coefficients for the quadratic equation
 *   that determines the NAT temperature.
 * - Solves the quadratic equation to find the NAT temperature.
 *
 * The calculations are based on empirical relationships involving logarithms of the 
 * partial pressures of H2O and HNO3.
 *
 * @note The constants and formulae used are specific to the context of atmospheric
 * chemistry and the formation of NAT.
 *
 * @author Lars Hoffmann
 */
double nat_temperature(
  const double p,
  const double h2o,
  const double hno3);

/**
 * @brief Reads air parcel data from a specified file into the given atmospheric structure.
 *
 * This function reads air parcel data from a file and populates the
 * provided `atm_t` structure based on the type of data specified in
 * the `ctl_t` control structure. It supports various data formats
 * including ASCII, binary, netCDF, and CLaMS.
 *
 * @param filename The name of the file containing the atmospheric data.
 * @param ctl      A pointer to the control structure (`ctl_t`) that specifies the type of data.
 * @param atm      A pointer to the atmospheric structure (`atm_t`) that will be populated with the data.
 * @return         Returns 1 on success, and 0 on failure.
 *
 * This function performs the following steps:
 * - Sets a timer for performance measurement.
 * - Initializes the atmospheric structure.
 * - Logs the file being read.
 * - Reads data from the file based on the specified type (`ctl->atm_type`):
 *   - `0` for ASCII data
 *   - `1` for binary data
 *   - `2` for netCDF data
 *   - `3` or `4` for CLaMS data
 * - Handles errors if the data type is not supported.
 * - Checks the result of the data reading function and ensures data was read successfully.
 * - Logs information about the number of air parcels and the ranges of various parameters (time, altitude, pressure, longitude, latitude, and other quantities).
 *
 * The function utilizes several helper functions and macros:
 * - `SELECT_TIMER` for setting the timer.
 * - `LOG` for logging information.
 * - `ERRMSG` for handling error messages.
 * - `gsl_stats_minmax` for calculating minimum and maximum values.
 * - `Z` for converting altitude.
 *
 * @author Lars Hoffmann
 */
int read_atm(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

/**
 * @brief Reads air parcel data from an ASCII file and populates the given atmospheric structure.
 *
 * This function reads air parcel data from an ASCII file and stores
 * the data in the provided `atm_t` structure. It reads each line of
 * the file, extracts the necessary data fields, and converts the
 * altitude to pressure.
 *
 * @param filename The name of the ASCII file containing the atmospheric data.
 * @param ctl      A pointer to the control structure (`ctl_t`) that specifies the number of quantities.
 * @param atm      A pointer to the atmospheric structure (`atm_t`) that will be populated with the data.
 * @return         Returns 1 on success, and 0 on failure.
 *
 * This function performs the following steps:
 * - Attempts to open the specified file for reading.
 * - Logs a warning and returns 0 if the file cannot be opened.
 * - Reads each line of the file and extracts data values for time, altitude, longitude, latitude, 
 *   and other specified quantities.
 * - Converts the altitude to pressure.
 * - Increments the data point counter.
 * - Checks if the number of data points exceeds the maximum allowed (`NP`) and logs an error message if so.
 * - Closes the file after reading all data.
 * - Returns 1 to indicate successful data reading.
 *
 * The function utilizes several macros and helper functions:
 * - `WARN` for logging warnings.
 * - `ERRMSG` for handling error messages.
 * - `TOK` for tokenizing and reading values from the line.
 * - `P` for converting altitude to pressure.
 *
 * @author Lars Hoffmann
 */
int read_atm_asc(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

/**
 * @brief Reads air parcel data from a binary file and populates the given atmospheric structure.
 *
 * This function reads air parcel data from a binary file and stores
 * the data in the provided `atm_t` structure. It checks the version
 * of the binary data, reads the data values, and verifies the
 * integrity of the data read.
 *
 * @param filename The name of the binary file containing the atmospheric data.
 * @param ctl      A pointer to the control structure (`ctl_t`) that specifies the number of quantities.
 * @param atm      A pointer to the atmospheric structure (`atm_t`) that will be populated with the data.
 * @return         Returns 1 on success, and 0 on failure.
 *
 * This function performs the following steps:
 * - Attempts to open the specified file for reading.
 * - Returns 0 if the file cannot be opened.
 * - Checks the version of the binary data and logs an error message if the version is incorrect.
 * - Reads the number of data points (`np`).
 * - Reads the data arrays for time, pressure, longitude, latitude, and other specified quantities.
 * - Checks a final flag to ensure the data was read correctly.
 * - Logs an error message if the final flag is incorrect.
 * - Closes the file after reading all data.
 * - Returns 1 to indicate successful data reading.
 *
 * The function utilizes several macros and helper functions:
 * - `ERRMSG` for handling error messages.
 * - `FREAD` for reading data from the binary file.
 *
 * @author Lars Hoffmann
 */
int read_atm_bin(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

/**
 * @brief Reads air parcel data from a CLaMS netCDF file and populates the given atmospheric structure.
 *
 * This function reads air parcel data from a CLaMS (Chemical
 * Lagrangian Model of the Stratosphere) netCDF file and stores the
 * data in the provided `atm_t` structure. It handles various
 * coordinate systems and ensures the necessary data is correctly read
 * and stored.
 *
 * @param filename The name of the netCDF file containing the atmospheric data.
 * @param ctl      A pointer to the control structure (`ctl_t`) that specifies the vertical coordinate system and quantities.
 * @param atm      A pointer to the atmospheric structure (`atm_t`) that will be populated with the data.
 * @return         Returns 1 on success, and 0 on failure.
 *
 * This function performs the following steps:
 * - Attempts to open the specified netCDF file for reading.
 * - Returns 0 if the file cannot be opened.
 * - Retrieves the number of particles (`np`) from the "NPARTS" dimension.
 * - Reads the initial time (`TIME_INIT`) if available, otherwise uses the "time" variable.
 * - Reads the vertical coordinate data based on the specified system (`vert_coord_ap`):
 *   - If `vert_coord_ap` is 1, reads "ZETA" and optionally "PRESS".
 *   - Otherwise, reads "PRESS_INIT" if available, otherwise uses "PRESS".
 * - Reads the longitude ("LON") and latitude ("LAT") data.
 * - Closes the netCDF file.
 * - Returns 1 to indicate successful data reading.
 *
 * The function utilizes several macros and helper functions:
 * - `NC_INQ_DIM` for inquiring about dimensions in the netCDF file.
 * - `NC_GET_DOUBLE` for reading double values from the netCDF file.
 * - `NC` for checking netCDF function return values.
 * - `WARN` for logging warnings.
 *
 * @author Jan Clemens
 */
int read_atm_clams(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

/**
 * @brief Reads air parcel data from a generic netCDF file and populates the given atmospheric structure.
 *
 * This function reads air parcel data from a netCDF file and stores
 * the data in the provided `atm_t` structure.  It retrieves the
 * dimensions, geolocations (time, pressure, longitude, latitude), and
 * specified variables from the file.
 *
 * @param filename The name of the netCDF file containing the atmospheric data.
 * @param ctl      A pointer to the control structure (`ctl_t`) that specifies the number of quantities and their names.
 * @param atm      A pointer to the atmospheric structure (`atm_t`) that will be populated with the data.
 * @return         Returns 1 on success, and 0 on failure.
 *
 * This function performs the following steps:
 * - Attempts to open the specified netCDF file for reading.
 * - Returns 0 if the file cannot be opened.
 * - Retrieves the number of observations (`np`) from the "obs" dimension.
 * - Reads the geolocation data arrays for time, pressure, longitude, and latitude.
 * - Reads the specified variables into the corresponding arrays in the `atm_t` structure.
 * - Closes the netCDF file after reading all data.
 * - Returns 1 to indicate successful data reading.
 *
 * The function utilizes several macros and helper functions:
 * - `NC_INQ_DIM` for inquiring about dimensions in the netCDF file.
 * - `NC_GET_DOUBLE` for reading double values from the netCDF file.
 * - `NC` for checking netCDF function return values.
 *
 * @author Lars Hoffmann
 */
int read_atm_nc(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

/**
 * @brief Reads various climatological data and populates the given climatology structure.
 *
 * This function reads a range of climatological datasets based on the
 * specified control settings and stores the data in the provided
 * `clim_t` structure. It handles initialization of tropopause
 * climatology, photolysis rates, and multiple gas species'
 * climatologies and time series.
 *
 * @param ctl   A pointer to the control structure (`ctl_t`) that specifies file names and parameters for climatology data.
 * @param clim  A pointer to the climatology structure (`clim_t`) that will be populated with the data.
 *
 * This function performs the following steps:
 * - Sets a timer for reading climatology data.
 * - Initializes the tropopause climatology.
 * - Reads photolysis rates if specified in `ctl`.
 * - Reads HNO3 climatology if specified in `ctl`.
 * - Reads OH climatology if specified in `ctl` and applies a diurnal correction if specified.
 * - Reads H2O2, HO2, O(1D) climatologies if specified in `ctl`.
 * - Reads time series data for various gases (CFC-10, CFC-11, CFC-12, N2O, SF6) if specified in `ctl`.
 *
 * The function utilizes several helper functions:
 * - `clim_tropo_init` for initializing tropopause climatology.
 * - `read_clim_photo` for reading photolysis rates.
 * - `read_clim_zm` for reading zonal mean climatologies.
 * - `clim_oh_diurnal_correction` for applying diurnal correction to OH climatology.
 * - `read_clim_ts` for reading time series data.
 *
 * @authors Lars Hoffmann
 * @authors Mingzhao Liu
 */
void read_clim(
  ctl_t * ctl,
  clim_t * clim);

/**
 * @brief Reads photolysis rates from a NetCDF file and populates the given photolysis structure.
 *
 * This function opens a NetCDF file specified by the filename, reads
 * various dimensions and data related to photolysis rates, and stores
 * this data in the provided `clim_photo_t` structure.  It includes
 * checks for data consistency and logs detailed information about the
 * loaded data.
 *
 * @param filename  A string containing the path to the NetCDF file containing photolysis rate data.
 * @param photo     A pointer to the photolysis structure (`clim_photo_t`) that will be populated with the data.
 *
 * The function performs the following steps:
 * - Logs the initiation of reading photolysis rates.
 * - Opens the NetCDF file in read-only mode.
 * - Reads pressure data and checks for descending order.
 * - Reads total column ozone data and checks for ascending order.
 * - Reads solar zenith angle data and checks for ascending order.
 * - Allocates memory for temporary arrays to hold the data.
 * - Reads various photolysis rates (e.g., J_N2O, J_CCl4, J_CFC-11, J_CFC-12, J_O2, J_O3b, J_O3a, J_H2O2, J_H2O) 
 *   and stores them in the `clim_photo_t` structure.
 * - Frees the allocated memory for temporary arrays.
 * - Closes the NetCDF file.
 * - Logs detailed information about the loaded data, including pressure levels, solar zenith angles, 
 *   and photolysis rates.
 *
 * @author Mingzhao Liu
 */
void read_clim_photo(
  const char *filename,
  clim_photo_t * photo);

/**
 * @brief Reads a climatological time series from a file and populates the given time series structure.
 *
 * This function reads time and volume mixing ratio (VMR) data from a
 * specified file, processes the data, and stores it in the provided
 * `clim_ts_t` structure. It also includes checks for data consistency
 * and logs detailed information about the loaded data.
 *
 * @param filename  A string containing the path to the file containing the climatological time series data.
 * @param ts        A pointer to the time series structure (`clim_ts_t`) that will be populated with the data.
 * @return          Returns 1 on success, and 0 on failure (e.g., if the file cannot be opened or data is invalid).
 *
 * The function performs the following steps:
 * - Logs the initiation of reading the climatological time series.
 * - Opens the file for reading.
 * - Reads time and VMR data from the file, converting years to seconds.
 * - Checks for ascending order of time data and ensures the number of data points does not exceed the limit.
 * - Closes the file after reading.
 * - Checks if there are enough data points.
 * - Logs detailed information about the loaded data, including the number of time steps and the range of VMR values.
 *
 * @author Lars Hoffmann
 */
int read_clim_ts(
  const char *filename,
  clim_ts_t * ts);

/**
 * @brief Reads zonally averaged climatological data from a netCDF file and populates the given structure.
 *
 * This function reads data from a specified netCDF file, including
 * pressure levels, latitudes, and volume mixing ratios (VMR) for a
 * specified variable. It performs necessary checks and logs detailed
 * information about the loaded data.
 *
 * @param filename  A string containing the path to the netCDF file.
 * @param varname   A string containing the name of the variable to be read from the netCDF file.
 * @param zm        A pointer to the structure (`clim_zm_t`) that will be populated with the data.
 *
 * The function performs the following steps:
 * - Logs the initiation of reading the specified data.
 * - Opens the netCDF file for reading.
 * - Reads pressure level data and checks for descending order.
 * - Reads latitude data and checks for ascending order.
 * - Sets the time data for monthly means.
 * - Checks the number of time steps.
 * - Reads the specified variable data from the file.
 * - Fixes any gaps in the data by interpolating from valid values.
 * - Logs detailed information about the loaded data, including the number of time steps, pressure levels, 
 *   latitude values, and the range of the variable's volume mixing ratios.
 *
 * @author Lars Hoffmann
 */
void read_clim_zm(
  const char *filename,
  char *varname,
  clim_zm_t * zm);

/**
 * @brief Reads control parameters from a configuration file and populates the given structure.
 *
 * This function reads control parameters from a specified
 * configuration file and command line arguments, populating the
 * provided `ctl_t` structure with the parsed data. It handles a wide
 * range of parameters, performing necessary checks and providing
 * default values where applicable.
 *
 * @param filename  A string containing the path to the configuration file.
 * @param argc      An integer representing the number of command line arguments.
 * @param argv      An array of strings containing the command line arguments.
 * @param ctl       A pointer to the structure (`ctl_t`) that will be populated with the control parameters.
 *
 * The function performs the following steps:
 * - Sets a timer for reading the control file.
 * - Logs information about the MPTRAC executable version and compilation details.
 * - Initializes quantity indices.
 * - Reads and sets various control parameters such as quantities, vertical coordinates, time steps,
 *   meteorological data, sorting options, isosurface parameters, random number generator type,
 *   advection parameters, diffusion parameters, convection parameters, boundary conditions,
 *   species parameters, molar mass, OH chemistry parameters, H2O2 chemistry parameters,
 *   KPP chemistry parameters, first order tracer chemistry parameters, wet deposition parameters,
 *   dry deposition parameters, climatological data paths, mixing parameters, chemistry grid parameters,
 *   exponential decay parameters, PSC analysis parameters, output parameters for atmospheric data,
 *   CSI data, ensemble data, grid data, profile data, sample data, station data, and VTK data.
 *
 * @author Lars Hoffmann
 */
void read_ctl(
  const char *filename,
  int argc,
  char *argv[],
  ctl_t * ctl);

/**
 * @brief Reads kernel function data from a file and populates the provided arrays.
 *
 * This function reads kernel function data from a specified file,
 * populating the provided arrays `kz` and `kw` with the parsed
 * data. It also updates the variable pointed to by `nk` with the
 * number of data points read. The function ensures that the height
 * levels are in ascending order and performs checks for the number of
 * height levels read.
 *
 * @param filename  A string containing the path to the file containing kernel function data.
 * @param kz        A double array to store the height levels of the kernel function.
 * @param kw        A double array to store the weights corresponding to the height levels.
 * @param nk        A pointer to an integer variable representing the number of data points read.
 *
 * The function performs the following steps:
 * - Logs information indicating the kernel function file being read.
 * - Attempts to open the specified file for reading.
 * - Reads data from the file line by line, parsing height levels and weights.
 * - Checks that the height levels are in ascending order and that the number of data points does
 *   not exceed the defined maximum.
 * - Closes the file after reading.
 * - Updates the value of `nk` with the number of data points read.
 * - Normalizes the kernel function weights by dividing each weight by the maximum weight.
 *
 * @author Lars Hoffmann
 */
void read_kernel(
  const char *filename,
  double kz[EP],
  double kw[EP],
  int *nk);

/**
 * @brief Reads meteorological data from a file and populates the provided structures.
 *
 * This function reads meteorological data from a specified file,
 * populating the provided structures `ctl`, `clim`, and `met` with
 * the parsed data. The function handles both netCDF and binary data
 * formats. For netCDF data, it performs various data processing steps
 * such as reading coordinates, surface data, atmospheric levels, and
 * calculating derived quantities. For binary data, it reads the data
 * directly into the appropriate arrays and structures.
 *
 * @param filename  A string containing the path to the file containing meteorological data.
 * @param ctl       A pointer to a structure containing control parameters.
 * @param clim      A pointer to a structure containing climatological data.
 * @param met       A pointer to a structure to store the meteorological data.
 * @return          Returns 1 on success and 0 on failure.
 *
 * The function performs the following steps:
 * - Logs information indicating the meteorological data file being read.
 * - Determines the type of meteorological data (netCDF or binary).
 * - For netCDF data:
 *   - Opens the netCDF file and reads various data components including grid coordinates,
 *     atmospheric levels, surface data, and derived quantities.
 *   - Performs data processing steps such as extrapolation, boundary condition handling,
 *     downsampling, and calculations of geopotential heights, potential vorticity, boundary
 *     layer data, tropopause data, cloud properties, convective available potential energy,
 *     total column ozone, and detrending.
 *   - Checks and adjusts meteorological data as necessary, including monotonicity of
 *     vertical coordinates.
 *   - Closes the netCDF file after reading.
 * - For binary data:
 *   - Opens the binary file and reads metadata including the type and version of the data,
 *     time information, dimensions, and grid coordinates.
 *   - Reads surface data, level data, and final flag from the binary file.
 *   - Closes the binary file after reading.
 * - Copies wind data to a cache for subsequent use.
 * - Returns 1 on successful reading and processing of meteorological data.
 *
 * @note The function handles different types of meteorological data formats and performs
 *       appropriate processing and error checking for each format.
 *
 * @author Lars Hoffmann
 */
int read_met(
  const char *filename,
  ctl_t * ctl,
  clim_t * clim,
  met_t * met);

/**
 * @brief Reads a 2-dimensional meteorological variable from a binary file and stores it in the provided array.
 *
 * This function reads a 2-dimensional meteorological variable from a
 * binary file, which is assumed to be uncompressed, and stores it in
 * the provided 2-dimensional array `var`. The variable name is used
 * for logging purposes to identify the data being read.
 *
 * @param in        A pointer to the FILE structure representing the binary file to read from.
 * @param met       A pointer to a structure containing meteorological data.
 * @param var       A 2-dimensional array to store the read variable.
 * @param varname   A string containing the name of the variable being read.
 *
 * The function performs the following steps:
 * - Allocates memory for a temporary buffer to hold the uncompressed data.
 * - Logs information about the variable being read from the file.
 * - Reads the uncompressed data from the file into the temporary buffer.
 * - Copies the data from the temporary buffer to the provided 2-dimensional array.
 * - Frees the memory allocated for the temporary buffer.
 *
 * @note The function assumes that the binary file contains uncompressed data and reads the data
 *       directly into the provided array without any additional processing.
 *
 * @author Lars Hoffmann
 */
void read_met_bin_2d(
  FILE * in,
  met_t * met,
  float var[EX][EY],
  char *varname);

/**
 * @brief Reads a 3-dimensional meteorological variable from a binary file and stores it in the provided array.
 *
 * This function reads a 3-dimensional meteorological variable from a
 * binary file, which can be uncompressed or compressed using
 * different methods, and stores it in the provided 3-dimensional
 * array `var`. The variable name is used for logging purposes to
 * identify the data being read.
 *
 * @param in        A pointer to the FILE structure representing the binary file to read from.
 * @param ctl       A pointer to a structure containing control parameters.
 * @param met       A pointer to a structure containing meteorological data.
 * @param var       A 3-dimensional array to store the read variable.
 * @param varname   A string containing the name of the variable being read.
 *
 * The function performs the following steps based on the compression type specified in `ctl->met_type`:
 * - If uncompressed data, it allocates memory for a temporary buffer, reads the data from the file,
 *   and then copies it to the provided array.
 * - If the data is compressed using packing, zfp, zstd, or cmultiscale compression, it calls corresponding
 *   compression functions to decompress the data, then copies the decompressed data to the provided array.
 * - The data is copied to the provided array using OpenMP parallelization for improved performance.
 * - Finally, the memory allocated for the temporary buffer is freed.
 *
 * @note This function supports different compression methods based on the value of `ctl->met_type`.
 *       If a particular compression method is not supported or the necessary libraries are not available,
 *       an error message is generated.
 *
 * @author Lars Hoffmann
 */
void read_met_bin_3d(
  FILE * in,
  ctl_t * ctl,
  met_t * met,
  float var[EX][EY][EP],
  char *varname);

/**
 * @brief Calculates Convective Available Potential Energy (CAPE) for each grid point.
 *
 * This function calculates the Convective Available Potential Energy
 * (CAPE) at each grid point based on the provided meteorological
 * data. CAPE is a measure of the energy available for deep
 * convection, which is essential for severe weather development.
 *
 * @param clim A pointer to a structure containing climatological data.
 * @param met  A pointer to a structure containing meteorological data.
 *
 * The function performs the following steps:
 * - Sets up a timer to monitor the calculation time.
 * - Initializes variables and constants required for the computation, such as vertical spacing and pressure levels.
 * - Iterates over each grid point in parallel using OpenMP.
 * - Calculates CAPE by integrating the difference in virtual temperatures between the environment and the parcel,
 *   up to the level of free convection (LFC).
 * - Determines the lifted condensation level (LCL), level of free convection (LFC), equilibrium level (EL),
 *   and Convective Inhibition (CIN) for each grid point.
 * - Checks the results and updates the corresponding fields in the meteorological data structure.
 *
 * @note The function utilizes OpenMP for parallelization to enhance performance by distributing
 *       the computation across multiple threads.
 *
 * @author Lars Hoffmann
 */
void read_met_cape(
  clim_t * clim,
  met_t * met);

/**
 * @brief Calculates cloud-related variables for each grid point.
 *
 * This function calculates cloud-related variables, such as cloud
 * cover, cloud top pressure, cloud bottom pressure, and total cloud
 * water content, based on the provided meteorological data.
 *
 * @param met A pointer to a structure containing meteorological data.
 *
 * The function performs the following steps:
 * - Sets up a timer to monitor the calculation time.
 * - Initializes variables and constants required for the computation.
 * - Iterates over each grid point in parallel using OpenMP.
 * - Determines cloud-related variables based on thresholds for liquid water content (LWC), rain water content (RWC), ice water content (IWC) and snow water content (SWC).
 * - Calculates cloud cover, cloud top pressure, cloud bottom pressure, and total cloud water content for each grid point.
 * - Updates the corresponding fields in the meteorological data structure.
 *
 * @note The function utilizes OpenMP for parallelization to enhance performance by distributing
 *       the computation across multiple threads.
 *
 * @author Lars Hoffmann
 */
void read_met_cloud(
  met_t * met);

/**
 * @brief Detrends meteorological data.
 *
 * This function detrends meteorological data by removing spatially
 * varying backgrounds from each grid point.  Detrending helps in
 * removing systematic biases and trends from the data, enabling
 * better analysis and modeling.
 *
 * @param ctl A pointer to a structure containing control parameters.
 * @param met A pointer to a structure containing meteorological data.
 *
 * The function performs the following steps:
 * - Checks if detrending is enabled based on the control parameters.
 * - Sets up a timer to monitor the detrending time.
 * - Allocates memory for a temporary meteorological data structure.
 * - Calculates the standard deviation and box size for detrending.
 * - Calculates the detrended data by subtracting spatially varying backgrounds.
 * - Updates the original meteorological data with the detrended values.
 * - Frees the allocated memory.
 *
 * @note Detrending is performed by subtracting spatially varying backgrounds calculated from neighboring grid points.
 * @note OpenMP is utilized for parallelization to enhance performance by distributing the computation across multiple threads.
 *
 * @author Lars Hoffmann
 */
void read_met_detrend(
  ctl_t * ctl,
  met_t * met);

/**
 * @brief Extrapolates meteorological data.
 *
 * This function extrapolates meteorological data by filling missing
 * or invalid data points with values from the nearest valid point
 * above.  Extrapolation is performed column-wise, ensuring that
 * missing data points are replaced with valid values.
 *
 * @param met A pointer to a structure containing meteorological data.
 *
 * The function performs the following steps:
 * - Sets up a timer to monitor the extrapolation time.
 * - Loops over each grid column in parallel.
 * - Finds the lowest valid data point within each column.
 * - Extrapolates missing or invalid data points by copying values from the nearest valid point above.
 * - Updates the meteorological data structure with the extrapolated values.
 *
 * @note Extrapolation is performed by copying values from the nearest valid point above to fill missing or invalid data points.
 *       OpenMP is utilized for parallelization to enhance performance by distributing the computation across multiple threads.
 *
 * @author Lars Hoffmann
 */
void read_met_extrapolate(
  met_t * met);

/**
 * @brief Calculates geopotential heights from meteorological data.
 *
 * This function calculates geopotential heights from provided
 * meteorological data using the hydrostatic equation.  Geopotential
 * heights are computed column-wise for each grid point based on the
 * temperature, pressure, and surface height information.  Optionally,
 * the calculated geopotential heights can be smoothed horizontally
 * using a weighted averaging scheme.
 *
 * @param ctl A pointer to a structure containing control parameters.
 * @param met A pointer to a structure containing meteorological data.
 *
 * The function performs the following steps:
 * - Sets up a timer to monitor the geopotential height calculation time.
 * - Calculates the logarithm of pressure levels for efficient computation.
 * - Applies the hydrostatic equation to determine geopotential heights based on temperature, pressure, and height information.
 * - Optionally, performs horizontal smoothing on the calculated geopotential heights.
 * - Updates the meteorological data structure with the computed geopotential heights.
 *
 * @note The hydrostatic equation is utilized to calculate geopotential heights, ensuring consistency with atmospheric conditions.
 *       Optionally, horizontal smoothing can be applied to the calculated geopotential heights to reduce spatial variability.
 *       OpenMP is utilized for parallelization to enhance performance by distributing the computation across multiple threads.
 *
 * @author Lars Hoffmann
 */
void read_met_geopot(
  ctl_t * ctl,
  met_t * met);

/**
 * @brief Reads meteorological grid information from a NetCDF file.
 *
 * This function reads meteorological grid information from a NetCDF
 * file, including time, spatial dimensions, and pressure levels.  It
 * also extracts longitudes, latitudes, and pressure levels from the
 * NetCDF file based on the specified control parameters.  The
 * function determines the time information either from the filename
 * or from the data file, depending on the file type.
 *
 * @param filename The filename of the NetCDF file.
 * @param ncid The NetCDF file identifier.
 * @param ctl A pointer to a structure containing control parameters.
 * @param met A pointer to a structure to store meteorological data.
 *
 * The function performs the following steps:
 * - Sets up a timer to monitor the reading time for meteorological grid information.
 * - Determines the time information from either the filename or the data file based on the file type.
 * - Checks the validity of the time information.
 * - Retrieves grid dimensions (longitude, latitude, and vertical levels) from the NetCDF file.
 * - Reads longitudes, latitudes, and pressure levels from the NetCDF file.
 * - Converts pressure levels to hPa if necessary.
 * - Logs the retrieved grid information for verification and debugging purposes.
 *
 * @note This function supports reading meteorological grid information from different types of NetCDF files, including MPTRAC and CLaMS.
 *       The time information is extracted either from the filename or from the data file, depending on the file type and control parameters.
 *       Spatial dimensions (longitude, latitude, and vertical levels) and pressure levels are retrieved from the NetCDF file.
 *
 * @authors Lars Hoffmann
 * @authors Jan Clemens
 */
void read_met_grid(
  const char *filename,
  int ncid,
  ctl_t * ctl,
  met_t * met);

/**
 * @brief Reads meteorological variables at different vertical levels from a NetCDF file.
 *
 * This function reads meteorological variables such as temperature,
 * wind components, specific humidity, ozone data, cloud parameters,
 * and cloud cover at various vertical levels from a NetCDF file.  The
 * function supports reading meteorological data from both MPTRAC and
 * CLaMS formats.  Depending on the file format, it reads specific
 * variables and performs necessary conversions or interpolations.
 *
 * @param ncid The NetCDF file identifier.
 * @param ctl A pointer to a structure containing control parameters.
 * @param met A pointer to a structure to store meteorological data.
 *
 * The function performs the following steps:
 * - Sets up a timer to monitor the reading time for meteorological level data.
 * - Reads meteorological variables from the NetCDF file based on the specified control parameters and file format.
 * - Handles specific variables differently depending on the file format, such as reading temperature, wind components, humidity, ozone data, and cloud parameters.
 * - Performs conversions or interpolations if necessary, such as converting specific humidity and ozone data from mixing ratio to volume mixing ratio.
 * - Transfers velocity components to model levels for diabatic advection if applicable.
 * - Reads pressure on model levels if specified in the control parameters.
 * - Performs vertical interpolation from model levels to pressure levels if needed.
 * - Checks the ordering of pressure levels to ensure they are in descending order.
 *
 * @note This function supports reading meteorological variables from NetCDF files in MPTRAC or CLaMS formats and handles specific variables differently based on the file format and control parameters.
 *       It performs necessary conversions or interpolations and ensures the correctness of pressure levels.
 *
 * @authors Lars Hoffmann
 * @authors Jan Clemens
 */
void read_met_levels(
  int ncid,
  ctl_t * ctl,
  met_t * met);

/**
 * @brief Interpolates meteorological data to specified pressure levels.
 *
 * This function interpolates meteorological data from model levels to pressure levels.
 * The interpolation is performed in parallel over the spatial grid defined in the 
 * meteorological data structure.
 *
 * @param[in] ctl A pointer to a control structure containing the number of pressure levels (`met_np`)
 *                and the pressure levels themselves (`met_p`).
 * @param[in] met A pointer to a meteorological data structure containing the grid dimensions 
 *                (`nx`, `ny`) and the pressure profile (`pl`).
 * @param[in, out] var A 3D array containing the meteorological variable to be interpolated.
 *                      On output, this array will contain the interpolated values at the specified
 *                      pressure levels.
 * @param[in] varname A string representing the name of the meteorological variable being interpolated.
 *
 * This function performs the following steps:
 * - Sets a timer for the operation.
 * - Logs the start of the interpolation process with the variable name.
 * - Loops over the spatial columns (grid points).
 * - For each column, copies the pressure profile.
 * - Interpolates the meteorological variable to the specified pressure levels.
 * - Copies the interpolated data back into the `var` array.
 *
 * @note The interpolation is performed in parallel using OpenMP.
 *
 * @author Lars Hoffmann
 */
void read_met_ml2pl(
  ctl_t * ctl,
  met_t * met,
  float var[EX][EY][EP],
  char *varname);

/**
 * @brief Makes zeta and pressure profiles monotone.
 *
 * This function ensures that zeta and pressure profiles are monotone
 * increasing and decreasing with altitude.  It iterates over each
 * grid point and each level to identify inversions and linearly
 * interpolate between them to maintain monotonicity.  The
 * interpolation is performed for both zeta and pressure profiles.
 *
 * @param met A pointer to a structure containing meteorological data.
 *
 * The function performs the following steps:
 * - Sets up a timer to monitor the processing time.
 * - Iterates over each grid point in parallel using OpenMP.
 * - Identifies inversions in both zeta and pressure profiles and interpolates linearly between them to make the profiles monotone increasing.
 *
 * @note This function is crucial for maintaining the physical consistency of meteorological profiles, ensuring accurate atmospheric simulations.
 *
 * @author Jan Clemens
 */
void read_met_monotonize(
  met_t * met);

/**
 * @brief Reads a 2-dimensional meteorological variable from a NetCDF file.
 *
 * This function reads a 2-dimensional meteorological variable from a
 * NetCDF file and stores it in a specified destination array.  It
 * supports both packed and unpacked data formats and handles missing
 * values and scaling factors accordingly.  The function also checks
 * the meteorological data layout to ensure correct data copying.
 *
 * @param ncid The NetCDF file ID.
 * @param varname The name of the variable to read.
 * @param varname2 An alternative name of the variable to read (in case varname is not found).
 * @param varname3 An alternative name of the variable to read (in case varname2 is not found).
 * @param varname4 An alternative name of the variable to read (in case varname3 is not found).
 * @param ctl A pointer to a structure containing control parameters.
 * @param met A pointer to a structure containing meteorological data.
 * @param dest The destination array to store the read data.
 * @param scl A scaling factor to apply to the read data.
 * @param init Flag indicating whether to initialize the destination array before reading.
 * @return Returns 1 on success, 0 on failure.
 *
 * The function performs the following steps:
 * - Checks if the specified variable exists in the NetCDF file.
 * - Reads packed data if scaling factors are available, otherwise reads unpacked data.
 * - Handles missing values and scaling factors appropriately.
 * - Copies the data to the destination array, applying the scaling factor if provided.
 *
 * @author Lars Hoffmann
 */
int read_met_nc_2d(
  int ncid,
  char *varname,
  char *varname2,
  char *varname3,
  char *varname4,
  ctl_t * ctl,
  met_t * met,
  float dest[EX][EY],
  float scl,
  int init);

/**
 * @brief Reads a 3-dimensional meteorological variable from a NetCDF file.
 *
 * This function reads a 3-dimensional meteorological variable from a
 * NetCDF file and stores it in a specified destination array. It
 * supports both packed and unpacked data formats and handles missing
 * values and scaling factors accordingly. The function also checks
 * the meteorological data layout to ensure correct data copying.
 *
 * @param ncid The NetCDF file ID.
 * @param varname The name of the variable to read.
 * @param varname2 An alternative name of the variable to read (in case varname is not found).
 * @param varname3 An alternative name of the variable to read (in case varname2 is not found).
 * @param varname4 An alternative name of the variable to read (in case varname3 is not found).
 * @param ctl A pointer to a structure containing control parameters.
 * @param met A pointer to a structure containing meteorological data.
 * @param dest The destination array to store the read data.
 * @param scl A scaling factor to apply to the read data.
 * @return Returns 1 on success, 0 on failure.
 *
 * The function performs the following steps:
 * - Checks if the specified variable exists in the NetCDF file.
 * - Reads packed data if scaling factors are available, otherwise reads unpacked data.
 * - Handles missing values and scaling factors appropriately.
 * - Copies the data to the destination array, applying the scaling factor if provided.
 *
 * @author Lars Hoffmann
 */
int read_met_nc_3d(
  int ncid,
  char *varname,
  char *varname2,
  char *varname3,
  char *varname4,
  ctl_t * ctl,
  met_t * met,
  float dest[EX][EY][EP],
  float scl);

/**
 * @brief Calculates the planetary boundary layer (PBL) height for each grid point.
 *
 * This function estimates the height of the planetary boundary layer
 * (PBL) based on various meteorological parameters.  The method used
 * is based on empirical relationships, such as those proposed by
 * Vogelezang and Holtslag (1996) or Seidel et al. (2012).  It
 * computes the PBL height by analyzing the vertical profiles of
 * temperature, wind speed, humidity, and pressure.
 *
 * @param met A pointer to a structure containing meteorological data.
 *
 * The function performs the following steps:
 * - Sets timer for performance monitoring.
 * - Iterates over each grid point to calculate the PBL height.
 * - Determines the bottom level of the PBL based on pressure and a specified thickness.
 * - Finds the lowest model level near the bottom of the PBL.
 * - Interpolates meteorological variables to the near-surface level.
 * - Computes virtual potential temperature (theta_v) at the surface.
 * - Initializes variables for Richardson number calculation.
 * - Loops over vertical levels to calculate Richardson number and identify the PBL height.
 * - Determines the PBL height based on the critical Richardson number criterion.
 * - Stores the calculated PBL height in the meteorological data structure.
 *
 * @note This function plays a crucial role in atmospheric modeling applications for characterizing the height of the boundary layer, which influences atmospheric dynamics and pollutant dispersion.
 *
 * @author Lars Hoffmann
 */
void read_met_pbl(
  met_t * met);

/**
 * @brief Applies periodic boundary conditions to meteorological data along longitudinal axis.
 *
 * This function applies periodic boundary conditions to
 * meteorological data along the longitudinal axis.  It checks if the
 * difference between the last and first longitudes and the difference
 * between the second and first longitudes are approximately equal to
 * 360 degrees, indicating periodicity.  If the condition is met, the
 * function increases the longitude counter, sets the longitude value
 * for the new grid point, and copies meteorological data from the
 * first grid point to the last grid point to ensure periodicity.
 *
 * @param met A pointer to a structure containing meteorological data.
 *
 * The function performs the following steps:
 * - Sets timer for performance monitoring.
 * - Checks if the difference between the last and first longitudes and the difference between the second and first longitudes
 *   are approximately equal to 360 degrees, indicating periodicity.
 * - If periodicity is confirmed:
 *   - Increases the longitude counter.
 *   - Sets the longitude value for the new grid point by adding the difference between the second and first longitudes
 *     to the longitude of the penultimate grid point.
 *   - Copies meteorological data from the first grid point to the last grid point to ensure periodicity:
 *     - Surface variables (e.g., pressure, temperature, wind speed, land-sea mask, sea surface temperature) are copied.
 *     - Meteorological variables at each pressure level are copied.
 *     - Meteorological variables at each hybrid pressure level are copied.
 *
 * @note This function is useful for generating continuous meteorological fields over a periodic domain, which is common in atmospheric modeling, especially for global simulations.
 *
 * @author Lars Hoffmann
 */
void read_met_periodic(
  met_t * met);

/**
 * @brief Applies a fix for polar winds in meteorological data.
 *
 * This function applies a fix for polar winds in meteorological data,
 * particularly focusing on the u and v wind components.  It checks if
 * the latitudes at the top and bottom of the grid are close to the
 * poles.  If so, it transforms the winds at 89-degree latitude into
 * Cartesian coordinates, takes their mean, and replaces the winds at
 * 90-degree latitude with this mean, effectively fixing the
 * unrealistic behavior of winds at the poles.
 *
 * @param met A pointer to a structure containing meteorological data.
 *
 * The function performs the following steps:
 * - Sets a timer for performance monitoring.
 * - Checks if the latitudes at the top and bottom of the grid are close to the poles (within 0.001 degree latitude of the poles).
 * - For each hemisphere (north and south):
 *   - Sets latitude indices for 89 degrees and 90 degrees.
 *   - Determines the sign of longitude adjustments based on the hemisphere.
 *   - Constructs lookup tables for cosine and sine of longitudes.
 *   - Loops over pressure levels and performs the following operations:
 *     - Transforms u and v wind components at 89 degrees latitude into Cartesian coordinates and calculates their mean.
 *     - Replaces u and v wind components at 90 degrees latitude with the calculated mean, effectively fixing the polar winds.
 *
 * @note This function is useful for correcting unrealistic behavior of winds near the poles in meteorological data, which can affect various atmospheric simulations.
 * @note Based on a Python code provided by Jens-Uwe Grooß.
 *
 * @author Lars Hoffmann
 */
void read_met_polar_winds(
  met_t * met);

/**
 * @brief Calculates potential vorticity (PV) from meteorological data.
 *
 * This function calculates the potential vorticity (PV) from the
 * provided meteorological data.  It employs finite difference methods
 * to estimate gradients of temperature, wind components, and pressure
 * in longitude, latitude, and pressure dimensions. These gradients
 * are then used to compute PV at each grid point.  Additionally, a
 * fix for polar regions is applied to ensure smoothness of PV values
 * in these regions.
 *
 * @param met A pointer to a structure containing meteorological data.
 *
 * The function performs the following steps:
 * - Sets a timer for performance monitoring.
 * - Computes powers for pressure calculation.
 * - Loops over grid points in longitude:
 *   - Sets latitude indices.
 *   - Loops over grid points in latitude:
 *     - Sets indices and auxiliary variables.
 *     - Loops over pressure levels:
 *       - Computes gradients in longitude, latitude, and pressure.
 *       - Calculates PV using computed gradients.
 * - Applies a fix for polar regions to ensure smoothness of PV values.
 *
 * @note Potential vorticity is a fundamental quantity in atmospheric dynamics, representing the potential of a fluid parcel to rotate due to changes in pressure, temperature, and wind fields.
 * @note Based on a Python code by Mathew Barlow (https://github.com/mathewbarlow/potential-vorticity).
 *
 * @author Lars Hoffmann
 */
void read_met_pv(
  met_t * met);

/**
 * @brief Calculates the total column ozone from meteorological ozone data.
 *
 * This function calculates the total column ozone from the provided
 * meteorological ozone data.  It integrates ozone concentrations over
 * altitude to obtain the column ozone density.  The result is then
 * converted to Dobson units, which represent the thickness of the
 * ozone layer if compressed into one layer at standard temperature
 * and pressure.
 *
 * @param met A pointer to a structure containing meteorological ozone data.
 *
 * The function performs the following steps:
 * - Sets a timer for performance monitoring.
 * - Loops over columns in longitude and latitude:
 *   - Integrates ozone concentrations over altitude.
 *   - Converts the integrated ozone density to Dobson units.
 *
 * @note Total column ozone is a critical metric for understanding ozone distribution in the atmosphere, with implications for climate, air quality, and UV radiation.
 *
 * @author Lars Hoffmann
 */
void read_met_ozone(
  met_t * met);

/**
 * @brief Downsamples meteorological data based on specified parameters.
 *
 * This function downsamples meteorological data based on the provided
 * control parameters.  It reduces the resolution of meteorological
 * data by averaging over specified intervals in longitude, latitude,
 * and altitude.
 *
 * @param ctl A pointer to a structure containing control parameters for downsampling.
 * @param met A pointer to a structure containing meteorological data to be downsampled.
 *
 * The function performs the following steps:
 * - Checks if downsampling parameters are set to a value less than or equal to 1, if so, returns without downsampling.
 * - Sets a timer for performance monitoring.
 * - Allocates memory for a temporary meteorological data structure.
 * - Copies metadata from the original structure to the temporary structure.
 * - Performs downsampling by smoothing over specified intervals:
 *   - Computes weighted averages over the specified intervals.
 *   - Updates the temporary structure with the smoothed values.
 * - Downsamples the smoothed data:
 *   - Updates longitude and latitude arrays with downsampled values.
 *   - Stores downsampled meteorological variables in the original structure.
 * - Frees memory allocated for the temporary structure.
 *
 * @note Downsampling meteorological data can be useful for reducing computational cost while preserving essential features for modeling and analysis.
 *
 * @author Lars Hoffmann
 */
void read_met_sample(
  ctl_t * ctl,
  met_t * met);

/**
 * @brief Reads surface meteorological data from a netCDF file and stores it in the meteorological data structure.
 *
 * This function reads various surface meteorological variables from a
 * netCDF file and stores them in the provided meteorological data
 * structure.  Depending on the configuration, it may read data for
 * surface pressure, geopotential height, temperature, zonal and
 * meridional wind components, land-sea mask, and sea surface
 * temperature.
 *
 * @param ncid NetCDF file identifier.
 * @param met A pointer to the meteorological data structure to store the read data.
 * @param ctl A pointer to a structure containing control parameters.
 *
 * The function performs the following steps:
 * - Sets a timer for performance monitoring.
 * - Reads surface meteorological data based on the configuration:
 *   - For MPTRAC meteorological data:
 *     - Reads surface pressure from "lnsp", "ps", or "sp" variables.
 *     - Converts surface pressure to Pa if necessary.
 *     - Reads geopotential height at the surface from "z" or "zm" variables.
 *     - Reads surface temperature from "t2m" or "2t" variables.
 *     - Reads zonal wind at the surface from "u10m" or "10u" variables.
 *     - Reads meridional wind at the surface from "v10m" or "10v" variables.
 *     - Reads land-sea mask from "lsm" variable.
 *     - Reads sea surface temperature from "sstk" or "sst" variables.
 *   - For CLaMS meteorological data:
 *     - Reads surface pressure from "ps" variable.
 *     - Reads geopotential height at the surface using the lowest level of the 3-D data field.
 *     - Reads surface temperature from "t2" variable.
 *     - Reads zonal wind at the surface from "u10" variable.
 *     - Reads meridional wind at the surface from "v10" variable.
 *     - Reads land-sea mask from "lsm" variable.
 *     - Reads sea surface temperature from "sstk" variable.
 *
 * @note The function handles different variable names and units according to the specified meteorological data source (MPTRAC or CLaMS).
 *
 * @authors Lars Hoffmann
 * @authors Jan Clemens
 */
void read_met_surface(
  int ncid,
  met_t * met,
  ctl_t * ctl);

/**
 * @brief Calculates the tropopause and related meteorological variables based on various methods and stores the results in the meteorological data structure.
 *
 * This function calculates the tropopause and related meteorological
 * variables using different methods specified by the control
 * parameters.  The calculated tropopause pressure is stored in the
 * provided meteorological data structure.
 *
 * @param ctl A pointer to a structure containing control parameters.
 * @param clim A pointer to the climatological data structure.
 * @param met A pointer to the meteorological data structure to store the calculated tropopause pressure and related variables.
 *
 * The function performs the following steps:
 * - Sets a timer for performance monitoring.
 * - Retrieves altitude and pressure profiles from the meteorological data structure.
 * - Depending on the control parameters (ctl->met_tropo), it calculates the tropopause using one of the following methods:
 *   - If ctl->met_tropo == 0, it does not calculate the tropopause and assigns NaN values to the tropopause pressure.
 *   - If ctl->met_tropo == 1, it uses tropopause climatology to estimate the tropopause pressure based on latitude and time.
 *   - If ctl->met_tropo == 2, it calculates the tropopause based on the cold point method, finding the altitude where the temperature is at a minimum.
 *   - If ctl->met_tropo == 3 or ctl->met_tropo == 4, it calculates the tropopause using the WMO definition, which involves identifying a sharp temperature lapse rate between two pressure levels.
 *   - If ctl->met_tropo == 5, it calculates the dynamical tropopause based on potential vorticity and potential temperature profiles.
 * - Interpolates temperature, geopotential height, and water vapor content to the tropopause pressure level using spatial interpolation.
 * - Stores the interpolated values in the meteorological data structure.
 *
 * @note The function supports parallelization using OpenMP directives to improve performance.
 *
 * @author Lars Hoffmann
 */
void read_met_tropo(
  ctl_t * ctl,
  clim_t * clim,
  met_t * met);

/**
 * @brief Reads observation data from a file and stores it in arrays.
 *
 * This function reads observation data from a specified file in
 * either ASCII or NetCDF format, depending on the value of the
 * OBS_TYPE control parameter. It stores the time, altitude,
 * longitude, latitude, and observation values in the provided arrays.
 *
 * @param filename The path to the observation data file.
 * @param ctl A pointer to a structure containing control parameters.
 * @param rt An array to store the time values of the observations.
 * @param rz An array to store the altitude values of the observations.
 * @param rlon An array to store the longitude values of the observations.
 * @param rlat An array to store the latitude values of the observations.
 * @param robs An array to store the observation values.
 * @param nobs A pointer to an integer variable to store the number of observations read.
 *
 * The function performs the following steps:
 * - Logs an informational message indicating the observation data file being read.
 * - Reads the observation data from the file based on the OBS_TYPE control parameter:
 *   - If ctl->obs_type == 0, it reads the data from an ASCII file using the read_obs_asc function.
 *   - If ctl->obs_type == 1, it reads the data from a NetCDF file using the read_obs_nc function.
 *   - If ctl->obs_type is neither 0 nor 1, it generates an error message indicating that the OBS_TYPE must be set to 0 or 1.
 * - Checks if the time values are in ascending order and generates an error message if not.
 * - Logs statistical information about the observation data, including the number of observations, time range, altitude range, longitude range, latitude range, and observation value range.
 *
 * @note The function assumes that the observation data file is formatted correctly and that the arrays provided have sufficient memory allocated to store the data.
 *
 * @author Lars Hoffmann
 * @author Mingzhao Liu
 */
void read_obs(
  const char *filename,
  ctl_t * ctl,
  double *rt,
  double *rz,
  double *rlon,
  double *rlat,
  double *robs,
  int *nobs);

/**
 * @brief Reads observation data from an ASCII file.
 *
 * This function reads observation data from a specified ASCII
 * file. It extracts time, altitude, longitude, latitude, and
 * observation values from each line of the file and stores them in
 * the provided arrays.
 *
 * @param filename The path to the ASCII file containing the observation data.
 * @param rt An array to store the time values of the observations.
 * @param rz An array to store the altitude values of the observations.
 * @param rlon An array to store the longitude values of the observations.
 * @param rlat An array to store the latitude values of the observations.
 * @param robs An array to store the observation values.
 * @param nobs A pointer to an integer variable to store the number of observations read.
 *
 * The function performs the following steps:
 * - Attempts to open the specified observation data file in read mode.
 * - Reads each line of the file and parses it to extract the time, altitude, longitude, latitude, and observation values using the sscanf function.
 * - Stores the extracted values in the respective arrays.
 * - Checks if the number of observations exceeds the maximum allowed limit (NOBS) and generates an error message if so.
 * - Closes the observation data file after reading all data.
 *
 * @note The function assumes that the observation data file is properly formatted and that the arrays provided have sufficient memory allocated to store the data.
 *
 * @author Lars Hoffmann
 */
void read_obs_asc(
  const char *filename,
  double *rt,
  double *rz,
  double *rlon,
  double *rlat,
  double *robs,
  int *nobs);

/**
 * @brief Reads observation data from a NetCDF file.
 *
 * This function reads observation data from a specified NetCDF
 * file. It extracts time, altitude, longitude, latitude, and
 * observation values from the variables in the NetCDF file and stores
 * them in the provided arrays.
 *
 * @param filename The path to the NetCDF file containing the observation data.
 * @param rt An array to store the time values of the observations.
 * @param rz An array to store the altitude values of the observations.
 * @param rlon An array to store the longitude values of the observations.
 * @param rlat An array to store the latitude values of the observations.
 * @param robs An array to store the observation values.
 * @param nobs A pointer to an integer variable to store the number of observations read.
 *
 * The function performs the following steps:
 * - Attempts to open the specified NetCDF file in read-only mode using the nc_open function.
 * - Queries the dimensions of the 'nobs' variable in the NetCDF file to determine the number of observations using the NC_INQ_DIM macro.
 * - Reads the 'time', 'alt', 'lon', 'lat', and 'obs' variables from the NetCDF file using the NC_GET_DOUBLE macro and stores them in the respective arrays.
 * - Closes the NetCDF file after reading all data using the nc_close function.
 *
 * @note The function assumes that the NetCDF file contains the required variables ('time', 'alt', 'lon', 'lat', 'obs') and that the arrays provided have sufficient memory allocated to store the data.
 *
 * @author Lars Hoffmann
 */
void read_obs_nc(
  const char *filename,
  double *rt,
  double *rz,
  double *rlon,
  double *rlat,
  double *robs,
  int *nobs);

/**
 * @brief Scans a control file or command-line arguments for a specified variable.
 *
 * This function scans either a control file or command-line arguments
 * for a specified variable name and retrieves its value.  It searches
 * for the variable name in the control file or command-line arguments
 * and returns its corresponding value.  If the variable is not found,
 * it returns a default value specified by the user.
 *
 * @param filename The name of the control file to be scanned. If NULL, only command-line arguments will be scanned.
 * @param argc The number of command-line arguments.
 * @param argv An array of command-line arguments.
 * @param varname The name of the variable to be searched.
 * @param arridx The index of the variable array, if applicable. Set to -1 if not an array.
 * @param defvalue The default value to be returned if the variable is not found.
 * @param value A pointer to a character array to store the retrieved value.
 * @return The retrieved value of the variable as a double.
 *
 * The function performs the following steps:
 * - Attempts to open the specified control file in read mode using the fopen function. If the filename ends with a '-', the file is not opened.
 * - Constructs the full variable name based on the variable name and array index provided.
 * - Reads data from the control file, searching for the full variable name. If found, it sets the contain flag to 1 and breaks the loop.
 * - Searches through the command-line arguments for the full variable name. If found, it sets the value and contain flag and breaks the loop.
 * - Closes the control file if opened.
 * - If the variable is not found, it sets the value to the default value provided or throws an error if no default value is provided.
 * - Writes the variable name and its value to the log.
 * - Copies the retrieved value to the value parameter if it is not NULL.
 * - Returns the retrieved value as a double after converting it from a string using the atof function.
 *
 * @note This function assumes that the variable names and their values in the control file or command-line arguments are separated by whitespace.
 *
 * @author Lars Hoffmann
 */
double scan_ctl(
  const char *filename,
  int argc,
  char *argv[],
  const char *varname,
  int arridx,
  const char *defvalue,
  char *value);

/**
 * @brief Calculates the sedimentation velocity of a particle in air.
 *
 * This function calculates the sedimentation velocity of a particle
 * in air using the given parameters.
 *
 * @param p The atmospheric pressure [hPa].
 * @param T The temperature [K].
 * @param rp The radius of the particle [microns].
 * @param rhop The density of the particle [kg/m^3].
 * @return The sedimentation velocity of the particle [m/s].
 *
 * The function performs the following steps:
 * - Converts the radius of the particle from microns to meters.
 * - Calculates the density of dry air using the given atmospheric pressure and temperature.
 * - Calculates the dynamic viscosity of air using Sutherland's formula.
 * - Calculates the thermal velocity of an air molecule using the given temperature.
 * - Calculates the mean free path of an air molecule.
 * - Computes the Knudsen number for air based on the ratio of mean free path to particle radius.
 * - Applies the Cunningham slip-flow correction factor to account for particle size.
 * - Computes the sedimentation velocity of the particle based on the difference in densities between the particle and air, incorporating the slip-flow correction.
 *
 * @note This function assumes that the ideal gas law and Stokes' law are applicable for calculating the sedimentation velocity of the particle.
 *
 * @author Lars Hoffmann
 */
double sedi(
  const double p,
  const double T,
  const double rp,
  const double rhop);

/**
 * @brief Performs spline interpolation or linear interpolation.
 *
 * This function interpolates a set of data points using either cubic
 * spline interpolation or linear interpolation, depending on the
 * specified method.
 *
 * @param x The array of x-coordinates of the data points.
 * @param y The array of y-coordinates of the data points.
 * @param n The number of data points.
 * @param x2 The array of x-coordinates where interpolation is required.
 * @param y2 The array to store the interpolated y-values.
 * @param n2 The number of points to interpolate.
 * @param method The interpolation method: 1 for cubic spline, 0 for linear interpolation.
 *
 * If the method is set to 1 (cubic spline interpolation):
 * - The function initializes a cubic spline interpolator using GSL.
 * - It interpolates the y-values at the specified x-coordinates using the spline.
 * - The interpolated y-values are stored in the provided y2 array.
 *
 * If the method is set to 0 (linear interpolation):
 * - The function performs linear interpolation between adjacent data points.
 * - It locates the interval where each interpolation point falls and calculates the interpolated y-value using linear interpolation.
 * - The interpolated y-values are stored in the provided y2 array.
 *
 * @note The x-coordinates in both arrays (x and x2) must be sorted in ascending order.
 *
 * @author Lars Hoffmann
 */
void spline(
  const double *x,
  const double *y,
  const int n,
  const double *x2,
  double *y2,
  const int n2,
  const int method);

/**
 * @brief Calculates the standard deviation of a set of data.
 *
 * This function calculates the standard deviation of a set of
 * floating-point data values.
 *
 * @param data Pointer to the array of data values.
 * @param n Number of data values in the array.
 * @return The standard deviation of the data values. If the number of data values is less than or equal to 0, returns 0.
 *
 * The standard deviation is calculated using the formula:
 *
 * \f[ \sigma = \sqrt{\frac{\sum_{i=1}^{n} (x_i - \bar{x})^2}{n}} \f]
 *
 * where:
 * - \f$ \sigma \f$ is the standard deviation,
 * - \f$ x_i \f$ is each data value,
 * - \f$ \bar{x} \f$ is the mean of the data values, and 
 * - \f$ n \f$ is the total number of data values.
 *
 * @author Lars Hoffmann
 */
float stddev(
  const float *data,
  const int n);

/**
 * @brief Calculates the solar zenith angle.
 *
 * This function calculates the solar zenith angle, which is the angle
 * between the zenith (straight up) and the line connecting the
 * observer to the center of the sun.
 *
 * @param sec Seconds elapsed since 2000-01-01T12:00Z.
 * @param lon Observer's longitude in degrees.
 * @param lat Observer's latitude in degrees.
 * @return The solar zenith angle in radians.
 *
 * The solar zenith angle is calculated based on the observer's position (longitude and latitude) and the time specified
 * in seconds elapsed since 2000-01-01T12:00Z.
 *
 * @note This function assumes that the input longitude and latitude are given in degrees.
 *
 * @author Lars Hoffmann
 */
double sza_calc(
  const double sec,
  const double lon,
  const double lat);

/**
 * @brief Converts time components to seconds since January 1, 2000, 12:00:00 UTC.
 *
 * This function calculates the number of seconds elapsed since
 * January 1, 2000, 12:00:00 UTC, based on the provided year, month,
 * day, hour, minute, and second. It also includes a fractional part
 * to represent the remaining seconds.
 *
 * @param year The year.
 * @param mon The month (1-12).
 * @param day The day of the month (1-31).
 * @param hour The hour of the day (0-23).
 * @param min The minute (0-59).
 * @param sec The second (0-59).
 * @param remain The fractional part of seconds.
 * @param jsec Pointer to store the calculated number of seconds since January 1, 2000, 12:00:00 UTC.
 *
 * The function calculates the time elapsed since January 1, 2000, 12:00:00 UTC, up to the specified time and includes
 * any fractional seconds indicated by the "remain" parameter.
 *
 * @note The function uses the timegm function, which is similar to mktime but operates in UTC.
 *
 * @author Lars Hoffmann
 */
void time2jsec(
  const int year,
  const int mon,
  const int day,
  const int hour,
  const int min,
  const int sec,
  const double remain,
  double *jsec);

/**
 * @brief Measures and reports elapsed time for named and grouped timers.
 *
 * The `timer` function measures elapsed time for a specified named
 * timer and an optional group of timers, accumulating time statistics
 * such as minimum, maximum, and mean elapsed times. It also provides
 * an option to log the timing statistics to an output.
 *
 * @param name A string representing the name of the timer.
 * @param group A string representing the group to which the timer belongs.
 * @param output An integer flag indicating whether to report the timing statistics (non-zero to report).
 *
 * The function keeps track of multiple timers and groups. When called, it:
 * - Gets the current time and calculates the elapsed time since the last call.
 * - Adds the elapsed time to the current timers' statistics.
 * - Reports the statistics if the `output` parameter is non-zero.
 * - Identifies the IDs of the next timer and group based on the provided `name` and `group`.
 * - Checks if the `name` and `group` are new, and if so, initializes them.
 * - Saves the starting time for the next measurement.
 *
 * @note The function uses OpenMP's `omp_get_wtime()` to get the current wall time.
 * @note The function maintains static arrays and variables to store timer names, groups, and statistics.
 * @note The maximum number of timers and groups is defined by the `NTIMER` macro.
 * 
 * @warning If the number of timers or groups exceeds `NTIMER`, the function will trigger an error message.
 * 
 * @author Lars Hoffmann
 */
void timer(
  const char *name,
  const char *group,
  int output);

/**
 * @brief Extracts and converts a timestamp from a filename to Julian seconds.
 *
 * The `time_from_filename` function parses a given filename to
 * extract a timestamp and converts it to Julian seconds. The
 * timestamp is expected to follow a specific format and position
 * within the filename, defined by the `offset` parameter.
 *
 * @param filename A string representing the filename containing the timestamp.
 * @param offset An integer indicating the position from the end of the filename where the timestamp starts.
 * 
 * @return The time in Julian seconds as a double.
 *
 * The function performs the following steps:
 * - Extracts the year, month, day, hour, and minute components of the timestamp from the filename using the given offset.
 * - Validates the extracted components to ensure they represent a valid date and time.
 * - Converts the validated date and time components to Julian seconds using the `time2jsec` function.
 * - Returns the computed time in Julian seconds.
 *
 * @note The expected format of the timestamp in the filename is `YYYY-MM-DD_HH-MM` (e.g., "2023-05-27_14-45").
 *
 * @warning If the extracted components do not represent a valid date and time, the function will trigger an error message.
 * 
 * @author Lars Hoffmann
 */
double time_from_filename(
  const char *filename,
  int offset);

/**
 * @brief Computes the weighting factor for a given pressure with respect to the tropopause.
 *
 * The `tropo_weight` function calculates a weighting factor that
 * indicates how much a given pressure `p` is influenced by the
 * tropopause based on climatological data.
 *
 * @param clim A pointer to a `clim_t` structure containing climatological tropopause data.
 * @param t The time parameter, typically representing the time of year or specific temporal context.
 * @param lat The latitude for which the weighting factor is being calculated.
 * @param p The pressure for which the weighting factor is to be computed.
 * 
 * @return A double representing the weighting factor.
 * 
 * The function performs the following steps:
 * - Calculates the tropopause pressure `pt` using the `clim_tropo` function.
 * - Determines the pressure range around the tropopause: `p1` (lower bound) and `p0` (upper bound).
 * - Computes the weighting factor based on the given pressure `p`:
 *   - Returns `1` if `p` is greater than `p0` (troposphere).
 *   - Returns `0` if `p` is less than `p1` (stratosphere).
 *   - Linearly interpolates between `1` and `0` for pressures between `p0` and `p1`.
 *
 * @note The `clim_tropo` function is assumed to provide the tropopause pressure based on climatological data.
 * @note The constants used in the pressure range calculation (`0.866877899` and its reciprocal) are specific to the function's logic.
 *
 * @author Lars Hoffmann
 */
double tropo_weight(
  const clim_t * clim,
  const double t,
  const double lat,
  const double p);

/**
 * @brief Writes air parcel data to a file in various formats.
 *
 * The `write_atm` function writes the air parcel data stored in the
 * `atm` structure to a file specified by `filename`. The format of
 * the output file is determined by the `atm_type_out` field in the
 * `ctl` control structure.
 *
 * @param filename A string representing the name of the file to write the data to.
 * @param ctl A pointer to a `ctl_t` structure containing control parameters.
 * @param atm A pointer to an `atm_t` structure containing atmospheric data.
 * @param t The current time, used for certain output formats.
 * 
 * The function performs the following steps:
 * - Sets a timer for the write operation using the `SELECT_TIMER` macro.
 * - Logs the beginning of the write operation with the specified filename.
 * - Depending on the `atm_type_out` value in the `ctl` structure, writes the data in one of the following formats:
 *   - ASCII (`atm_type_out == 0`): Calls `write_atm_asc`.
 *   - Binary (`atm_type_out == 1`): Calls `write_atm_bin`.
 *   - netCDF (`atm_type_out == 2`): Calls `write_atm_nc`.
 *   - CLaMS trajectory (`atm_type_out == 3`): Calls `write_atm_clams_traj`.
 *   - CLaMS position (`atm_type_out == 4`): Calls `write_atm_clams`.
 * - If the `atm_type_out` value is not supported, triggers an error message.
 * - Logs various statistics about the atmospheric data, including the number of particles,
 *   time range, altitude range, pressure range, longitude range, and latitude range.
 * - Logs the range for each quantity specified in the `ctl` structure.
 *
 * @author Lars Hoffmann
 */
void write_atm(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm,
  double t);

/**
 * @brief Writes air parcel data to an ASCII file or gnuplot.
 *
 * The `write_atm_asc` function writes the atmospheric data stored in
 * the `atm` structure to an ASCII file specified by `filename` or to
 * pipe to gnuplot if requested.
 *
 * @param filename A string representing the name of the file to write the data to.
 * @param ctl A pointer to a `ctl_t` structure containing control parameters.
 * @param atm A pointer to an `atm_t` structure containing atmospheric data.
 * @param t The current time used for filtering and timestamping.
 * 
 * The function performs the following steps:
 * - Sets the time interval for the output data based on the control parameters.
 * - Checks if gnuplot output is requested and, if so, creates a pipe to gnuplot and sets up the plot.
 * - If gnuplot output is not requested, creates an ASCII file for writing.
 * - Writes the header information to the output file, including the description of each column.
 * - Iterates over the particles in the `atm` structure, filtering by time if specified, and writes the data to the output file.
 * - Closes the output file or gnuplot pipe.
 *
 * @author Lars Hoffmann
 */
void write_atm_asc(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm,
  double t);

/**
 * @brief Writes air parcel data to a binary file.
 *
 * The `write_atm_bin` function writes the air parcel data stored in
 * the `atm` structure to a binary file specified by `filename`. The
 * function includes versioning information and ensures that all
 * relevant data arrays are written in a consistent binary format.
 *
 * @param filename A string representing the name of the file to write the data to.
 * @param ctl A pointer to a `ctl_t` structure containing control parameters.
 * @param atm A pointer to an `atm_t` structure containing atmospheric data.
 * 
 * The function performs the following steps:
 * - Creates the binary file for writing. If the file cannot be created, it triggers an error message.
 * - Writes a version number for the binary data format.
 * - Writes the number of particles to the file.
 * - Writes the time, pressure, longitude, and latitude arrays to the file.
 * - Iterates over the quantities specified in the `ctl` structure and writes each quantity array to the file.
 * - Writes a final flag to indicate the end of the binary data.
 * - Closes the file.
 *
 * @author Lars Hoffmann
 */
void write_atm_bin(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

/**
 * @brief Writes air parcel data to a NetCDF file in the CLaMS format.
 *
 * The `write_atm_clams` function creates a NetCDF file and writes air
 * parcel data into it. The data includes time, latitude, longitude,
 * pressure, and other specified quantities. The function defines the
 * dimensions and variables, sets global attributes, and writes the
 * data to the file.
 *
 * @param filename A string representing the name of the output file.
 * @param ctl A pointer to a `ctl_t` structure containing control parameters.
 * @param atm A pointer to an `atm_t` structure containing atmospheric data.
 *
 * The function performs the following steps:
 * - Creates the NetCDF file with the specified filename.
 * - Defines the dimensions for time and the number of particles (NPARTS).
 * - Defines variables for time, latitude, longitude, pressure, zeta, and other quantities.
 * - Sets global attributes for the vertical coordinate name and model.
 * - Writes the data into the NetCDF file.
 * - Closes the NetCDF file after writing.
 *
 * @author Jan Clemens
 */
void write_atm_clams(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

/**
 * @brief Writes CLaMS trajectory data to a NetCDF file.
 *
 * The `write_atm_clams_traj` function writes trajectory data for the
 * CLaMS model to a NetCDF file. The file is created and populated
 * with data including time, latitude, longitude, pressure, and other
 * quantities. The function also handles the creation of a final
 * initialization file at the last time step.
 *
 * @param dirname A string representing the directory name where the file will be created.
 * @param ctl A pointer to a `ctl_t` structure containing control parameters.
 * @param atm A pointer to an `atm_t` structure containing atmospheric data.
 * @param t The current time in seconds since a reference epoch.
 *
 * The function performs the following steps:
 * - Determines the start and stop times of the calculation.
 * - Constructs the output filename based on the start and stop times.
 * - Defines the hyperslab for the trajectory file.
 * - Creates the NetCDF file if it's the first time step and defines dimensions and variables.
 * - Writes the trajectory data to the NetCDF file.
 * - At the last time step, creates an initialization file with the final data.
 *
 * @author Jan Clemens
 */
void write_atm_clams_traj(
  const char *dirname,
  ctl_t * ctl,
  atm_t * atm,
  double t);

/**
 * @brief Writes air parcel data to a NetCDF file.
 *
 * The `write_atm_nc` function creates a NetCDF file and writes air
 * parcel data into it.  The data includes time, pressure, longitude,
 * latitude, and other specified quantities.  The function defines the
 * dimensions and variables, sets global attributes, and writes the
 * data to the file.
 *
 * @param filename A string representing the name of the output file.
 * @param ctl A pointer to a `ctl_t` structure containing control parameters.
 * @param atm A pointer to an `atm_t` structure containing atmospheric data.
 *
 * The function performs the following steps:
 * - Creates the NetCDF file with the specified filename.
 * - Defines the dimension for the number of observations (obs).
 * - Defines variables for time, pressure, longitude, latitude, and other quantities.
 * - Sets global attributes for the feature type.
 * - Writes the data into the NetCDF file.
 * - Closes the NetCDF file after writing.
 *
 * @author Lars Hoffmann
 */
void write_atm_nc(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

/**
 * @brief Writes Critical Success Index (CSI) data to a file.
 *
 * The `write_csi` function processes air parcel and observation data
 * to calculate and write various verification statistics, including
 * the Critical Success Index (CSI), to a specified output file at
 * regular intervals. The statistics include measures such as the
 * number of hits, misses, and false alarms, bias, probability of
 * detection, false alarm rate, equitable threat score, and
 * correlation coefficients.
 *
 * @param filename A string representing the name of the output file.
 * @param ctl A pointer to a `ctl_t` structure containing control parameters.
 * @param atm A pointer to an `atm_t` structure containing atmospheric data.
 * @param t A double representing the current time.
 *
 * The function performs the following steps:
 * - Initializes resources and sets up the output file if the current time is the start time.
 * - Reads observation data and kernel data if provided.
 * - Sets grid box sizes and horizontal coordinates.
 * - Allocates memory for mean and count arrays.
 * - Loops over observations and model data to accumulate mean values and counts.
 * - Analyzes the grid cells to calculate CSI and other statistics.
 * - Writes the calculated statistics to the output file at specified intervals.
 * - Frees allocated resources and closes the file when the processing is complete.
 *
 * @author Lars Hoffmann
 */
void write_csi(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm,
  double t);

/**
 * @brief Writes ensemble data to a file.
 *
 * The `write_ens` function processes air parcel data to calculate
 * ensemble means and standard deviations for various quantities and
 * writes them to a specified output file. It handles ensemble members
 * and calculates statistics such as means and standard deviations for
 * each ensemble, along with latitude, longitude, altitude, and time
 * information.
 *
 * @param filename A string representing the name of the output file.
 * @param ctl A pointer to a `ctl_t` structure containing control parameters.
 * @param atm A pointer to an `atm_t` structure containing atmospheric data.
 * @param t A double representing the current time.
 *
 * The function performs the following steps:
 * - Initializes resources and sets up necessary variables.
 * - Sets a time interval for processing data.
 * - Loops over air parcels to accumulate means and standard deviations 
 *   for each ensemble member.
 * - Creates an output file and writes header information.
 * - Writes ensemble data, including time, altitude, latitude, longitude, 
 *   means, standard deviations, and the number of members.
 * - Closes the output file.
 *
 * @author Lars Hoffmann
 */
void write_ens(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm,
  double t);

/**
 * @brief Writes grid data to a file in ASCII or netCDF format.
 *
 * The `write_grid` function processes air parcel data to calculate
 * various grid-based statistics such as column density, mean, and
 * standard deviation for specified quantities. It then writes this
 * data to a specified output file either in ASCII or netCDF format
 * based on the configuration parameters provided in the `ctl`
 * structure.
 *
 * @param filename A string representing the name of the output file.
 * @param ctl A pointer to a `ctl_t` structure containing control parameters.
 * @param met0 A pointer to a `met_t` structure containing meteorological data
 *             for the initial time step.
 * @param met1 A pointer to a `met_t` structure containing meteorological data
 *             for the final time step.
 * @param atm A pointer to an `atm_t` structure containing atmospheric data.
 * @param t A double representing the current time.
 *
 * The function performs the following steps:
 * - Initializes resources and sets up necessary variables.
 * - Reads kernel data if it is specified in the control parameters.
 * - Allocates memory for various arrays to store grid data.
 * - Determines the grid box size and sets up vertical and horizontal
 *   coordinates.
 * - Sets a time interval for output data processing.
 * - Calculates grid box indices for atmospheric model data.
 * - Averages data within each grid box.
 * - Calculates column density and volume mixing ratio.
 * - Writes data to the output file either in ASCII or netCDF format based
 *   on the specified `grid_type` in the control parameters.
 * - Frees allocated memory.
 *
 * @note The function supports parallel processing using OpenMP for certain
 *       computational tasks to improve performance.
 *
 * @author Lars Hoffmann
 */
void write_grid(
  const char *filename,
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double t);

/**
 * @brief Writes grid data to an ASCII file.
 *
 * The `write_grid_asc` function writes gridded air parcel data,
 * including column density, mean and standard deviation for specified
 * quantities, and volume mixing ratio (if available), to an ASCII
 * file. The function also supports writing gnuplot commands to
 * generate plots if requested in the control parameters.
 *
 * @param filename A string representing the name of the output file.
 * @param ctl A pointer to a `ctl_t` structure containing control parameters.
 * @param cd An array of doubles representing column density values.
 * @param mean An array of arrays of doubles representing the mean values for
 *             specified quantities.
 * @param sigma An array of arrays of doubles representing the standard
 *              deviation values for specified quantities.
 * @param vmr_impl An array of doubles representing the volume mixing ratio
 *                 (implicit) values.
 * @param t A double representing the current time.
 * @param z An array of doubles representing vertical coordinates (altitude).
 * @param lon An array of doubles representing longitudinal coordinates.
 * @param lat An array of doubles representing latitudinal coordinates.
 * @param area An array of doubles representing surface area values.
 * @param dz A double representing the layer depth.
 * @param np An array of integers representing the number of particles.
 *
 * The function performs the following steps:
 * - Checks if gnuplot output is requested in the control parameters and sets
 *   up a gnuplot pipe if needed.
 * - If gnuplot output is requested, sets the plot filename and time string,
 *   and dumps gnuplot file contents to the pipe.
 * - Otherwise, creates the output file for writing in ASCII format.
 * - Writes the header information to the output file, including column labels.
 * - Writes the grid data to the output file, including time, altitude,
 *   coordinates, surface area, layer depth, column density, volume mixing
 *   ratio, number of particles, mean values for specified quantities, and
 *   standard deviation values if requested.
 * - Closes the output file.
 *
 * @note The function supports writing gnuplot commands to generate plots if
 *       requested in the control parameters. It also supports writing mean and
 *       standard deviation values for specified quantities if requested.
 *
 * @author Lars Hoffmann
 */
void write_grid_asc(
  const char *filename,
  ctl_t * ctl,
  double *cd,
  double *mean[NQ],
  double *sigma[NQ],
  double *vmr_impl,
  double t,
  double *z,
  double *lon,
  double *lat,
  double *area,
  double dz,
  int *np);

/**
 * @brief Writes grid data to a NetCDF file.
 *
 * The `write_grid_nc` function writes gridded air parcel data,
 * including column density, mean and standard deviation for specified
 * quantities, and volume mixing ratio (if available), to a NetCDF
 * file. NetCDF is a self-describing, machine-independent data format
 * for storing scientific data.
 *
 * @param filename A string representing the name of the output file.
 * @param ctl A pointer to a `ctl_t` structure containing control parameters.
 * @param cd An array of doubles representing column density values.
 * @param mean An array of arrays of doubles representing the mean values for
 *             specified quantities.
 * @param sigma An array of arrays of doubles representing the standard
 *              deviation values for specified quantities.
 * @param vmr_impl An array of doubles representing the volume mixing ratio
 *                 (implicit) values.
 * @param t A double representing the current time.
 * @param z An array of doubles representing vertical coordinates (altitude).
 * @param lon An array of doubles representing longitudinal coordinates.
 * @param lat An array of doubles representing latitudinal coordinates.
 * @param area An array of doubles representing surface area values.
 * @param dz A double representing the layer depth.
 * @param np An array of integers representing the number of particles.
 *
 * The function performs the following steps:
 * - Allocates memory for temporary arrays required for writing data.
 * - Creates a NetCDF file with the specified filename.
 * - Defines dimensions and variables in the NetCDF file, along with their
 *   attributes.
 * - Writes the data arrays to the NetCDF file.
 * - Closes the NetCDF file.
 * - Frees allocated memory.
 *
 * @note NetCDF files are commonly used in scientific computing and can be
 *       accessed by various programming languages and software packages.
 *       Additionally, the function supports writing mean and standard
 *       deviation values for specified quantities if requested.
 *
 * @author Lars Hoffmann
 */
void write_grid_nc(
  const char *filename,
  ctl_t * ctl,
  double *cd,
  double *mean[NQ],
  double *sigma[NQ],
  double *vmr_impl,
  double t,
  double *z,
  double *lon,
  double *lat,
  double *area,
  double dz,
  int *np);

/**
 * @brief Writes meteorological data to a binary file.
 *
 * The `write_met` function writes meteorological data, including
 * surface data and level data, to a binary file specified by the
 * filename parameter. The data is written in binary format for
 * efficient storage and transfer.
 *
 * @param filename A string representing the name of the output file.
 * @param ctl A pointer to a `ctl_t` structure containing control parameters.
 * @param met A pointer to a `met_t` structure containing meteorological data.
 *
 * The function performs the following steps:
 * - Sets a timer to measure the execution time.
 * - Writes information about the meteorological data being written to the log.
 * - Checks if compression flags are enabled and displays error messages if
 *   required compression methods are not available.
 * - Writes binary data to the file, including the type and version of the data,
 *   grid data (time, dimensions, coordinates), surface data, and level data.
 * - Closes the file after writing all data.
 *
 * @return Returns an integer value indicating the success or failure of the
 *         function. A value of 0 indicates successful execution.
 *
 * @note This function supports writing meteorological data in different binary
 *       formats, including uncompressed and compressed formats using various
 *       compression algorithms such as packing, ZFP, ZSTD, and cmultiscale. The specific
 *       compression method used depends on the configuration specified in the
 *       `ctl_t` structure.
 *
 * @author Lars Hoffmann
 */
int write_met(
  const char *filename,
  ctl_t * ctl,
  met_t * met);

/**
 * @brief Writes a 2-dimensional meteorological variable to a binary file.
 *
 * The `write_met_bin_2d` function writes a 2-dimensional
 * meteorological variable to a binary file specified by the `out`
 * parameter. The variable data is provided in a 2-dimensional array
 * `var` with maximum dimensions `EX` by `EY`.  The variable name is
 * provided as a string in the `varname` parameter.
 *
 * @param out A pointer to a FILE structure representing the output file.
 * @param met A pointer to a `met_t` structure containing meteorological data.
 * @param var An array of floats representing the 2-dimensional variable data.
 * @param varname A string containing the name of the variable being written.
 *
 * The function performs the following steps:
 * - Allocates memory for a temporary buffer to hold the variable data.
 * - Copies the variable data from the 2-dimensional array `var` to the temporary
 *   buffer `help`.
 * - Writes the uncompressed variable data to the binary file specified by `out`.
 * - Logs a message indicating the successful writing of the variable data.
 * - Frees the allocated memory.
 *
 * @note This function is typically used to write surface data or other
 *       2-dimensional meteorological variables to a binary file.
 *
 * @author Lars Hoffmann
 */
void write_met_bin_2d(
  FILE * out,
  met_t * met,
  float var[EX][EY],
  char *varname);

/**
 * @brief Writes a 3-dimensional meteorological variable to a binary file.
 *
 * The `write_met_bin_3d` function writes a 3-dimensional
 * meteorological variable to a binary file specified by the `out`
 * parameter. The variable data is provided in a 3-dimensional array
 * `var` with maximum dimensions `EX` by `EY` by `EP`. The variable
 * name is provided as a string in the `varname` parameter.
 * Additionally, the function takes parameters for specifying the
 * compression precision and tolerance.
 *
 * @param out A pointer to a FILE structure representing the output file.
 * @param ctl A pointer to a `ctl_t` structure containing control parameters.
 * @param met A pointer to a `met_t` structure containing meteorological data.
 * @param var An array of floats representing the 3-dimensional variable data.
 * @param varname A string containing the name of the variable being written.
 * @param precision An integer specifying the precision of compression (for certain compression methods).
 * @param tolerance A double specifying the tolerance for compression (for certain compression methods).
 *
 * The function performs the following steps:
 * - Allocates memory for a temporary buffer to hold the variable data.
 * - Copies the variable data from the 3-dimensional array `var` to the temporary
 *   buffer `help`.
 * - Writes the variable data to the binary file specified by `out` using the specified compression method
 *   (uncompressed, packed, zfp, zstd, cmultiscale).
 * - Logs a message indicating the successful writing of the variable data.
 * - Frees the allocated memory.
 *
 * @note This function is typically used to write level data or other
 *       3-dimensional meteorological variables to a binary file.
 *
 * @note Depending on the value of `ctl->met_type`, the function writes the
 *       variable data using different compression methods. If `ctl->met_type`
 *       is not supported, an error message is logged.
 *
 * @author Lars Hoffmann
 */
void write_met_bin_3d(
  FILE * out,
  ctl_t * ctl,
  met_t * met,
  float var[EX][EY][EP],
  char *varname,
  int precision,
  double tolerance);

/**
 * @brief Writes various types of output data to files in a specified directory.
 *
 * The `write_output` function writes various types of output data to
 * files in the directory specified by the `dirname` parameter. The
 * function takes control parameters (`ctl`), two meteorological data
 * structures (`met0` and `met1`), an atmospheric data structure
 * (`atm`), and a time value (`t`) as input.
 *
 * @param dirname A string representing the directory path where output files will be written.
 * @param ctl A pointer to a `ctl_t` structure containing control parameters.
 * @param met0 A pointer to a `met_t` structure representing the first set of meteorological data.
 * @param met1 A pointer to a `met_t` structure representing the second set of meteorological data.
 * @param atm A pointer to an `atm_t` structure representing atmospheric data.
 * @param t A double value representing the time at which the output is being written.
 *
 * The function performs the following steps:
 * - Parses the input time (`t`) to extract year, month, day, hour, minute, and second.
 * - Updates host memory if necessary based on control parameters.
 * - Writes atmospheric data to files if specified by control parameters.
 * - Writes gridded data to files if specified by control parameters.
 * - Writes CSI (Critical Success Index) data to files if specified by control parameters.
 * - Writes ensemble data to files if specified by control parameters.
 * - Writes profile data to files if specified by control parameters.
 * - Writes sample data to files if specified by control parameters.
 * - Writes station data to files if specified by control parameters.
 * - Writes VTK (Visualization Toolkit) data to files if specified by control parameters.
 *
 * @note This function orchestrates the writing of various types of output data to files
 *       based on control parameters and the current simulation time.
 *
 * @author Lars Hoffmann
 */
void write_output(
  const char *dirname,
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double t);

/**
 * @brief Writes profile data to a specified file.
 *
 * The `write_prof` function writes profile data to a file specified
 * by the `filename` parameter. It takes control parameters (`ctl`),
 * two meteorological data structures (`met0` and `met1`), an
 * atmospheric data structure (`atm`), and a time value (`t`) as
 * input.
 *
 * @param filename A string representing the filename where the profile data will be written.
 * @param ctl A pointer to a `ctl_t` structure containing control parameters.
 * @param met0 A pointer to a `met_t` structure representing the first set of meteorological data.
 * @param met1 A pointer to a `met_t` structure representing the second set of meteorological data.
 * @param atm A pointer to an `atm_t` structure representing atmospheric data.
 * @param t A double value representing the time at which the profile data is being written.
 *
 * The function performs the following steps:
 * - Initializes variables and allocates memory if it's the start of the simulation.
 * - Reads observation data and creates a new output file if necessary.
 * - Writes header information to the output file.
 * - Sets grid box size and vertical coordinates.
 * - Processes observations and model data within the specified time interval.
 * - Calculates and writes output data for each grid cell.
 * - Finalizes by closing the output file and freeing allocated memory if it's the end of the simulation.
 *
 * @note This function writes profile data to a file, including time, altitude, coordinates,
 *       atmospheric properties, observed data, and the number of observations.
 *
 * @author Lars Hoffmann
 */
void write_prof(
  const char *filename,
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double t);

/**
 * @brief Writes sample data to a specified file.
 *
 * The `write_sample` function writes sample data to a file specified
 * by the `filename` parameter. It takes control parameters (`ctl`),
 * two meteorological data structures (`met0` and `met1`), an
 * atmospheric data structure (`atm`), and a time value (`t`) as
 * input.
 *
 * @param filename A string representing the filename where the sample data will be written.
 * @param ctl A pointer to a `ctl_t` structure containing control parameters.
 * @param met0 A pointer to a `met_t` structure representing the first set of meteorological data.
 * @param met1 A pointer to a `met_t` structure representing the second set of meteorological data.
 * @param atm A pointer to an `atm_t` structure representing atmospheric data.
 * @param t A double value representing the time at which the sample data is being written.
 *
 * The function performs the following steps:
 * - Initializes variables and allocates memory if it's the start of the simulation.
 * - Reads observation data and kernel data if necessary.
 * - Creates a new output file and writes header information to it.
 * - Sets latitude range, squared radius, and area.
 * - Processes observations and calculates sample data within the specified time interval.
 * - Writes output data for each observation.
 * - Finalizes by closing the output file and freeing allocated memory if it's the end of the simulation.
 *
 * @note This function writes sample data to a file, including time, altitude, coordinates,
 *       surface area, layer depth, number of particles, column density, volume mixing ratio,
 *       and observed data.
 *
 * @author Lars Hoffmann
 */
void write_sample(
  const char *filename,
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double t);

/**
 * @brief Writes station data to a specified file.
 *
 * The `write_station` function writes station data to a file
 * specified by the `filename` parameter. It takes control parameters
 * (`ctl`), an atmospheric data structure (`atm`), and a time value
 * (`t`) as input.
 *
 * @param filename A string representing the filename where the station data will be written.
 * @param ctl A pointer to a `ctl_t` structure containing control parameters.
 * @param atm A pointer to an `atm_t` structure representing atmospheric data.
 * @param t A double value representing the time at which the station data is being written.
 *
 * The function performs the following steps:
 * - Initializes variables and opens a new file if it's the start of the simulation.
 * - Writes header information to the output file.
 * - Sets geolocation and search radius for station data.
 * - Processes air parcels and writes station data within the specified time interval and search radius.
 * - Writes station data for each air parcel satisfying the criteria.
 * - Closes the output file if it's the end of the simulation.
 *
 * @note This function writes station data to a file, including time, altitude, longitude, latitude,
 *       and additional quantities specified in the control parameters.
 *
 * @author Lars Hoffmann
 */
void write_station(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm,
  double t);

/**
 * @brief Writes VTK (Visualization Toolkit) data to a specified file.
 *
 * The `write_vtk` function writes VTK data to a file specified by the
 * `filename` parameter. It takes control parameters (`ctl`), an
 * atmospheric data structure (`atm`), and a time value (`t`) as
 * input.
 *
 * @param filename A string representing the filename where the VTK data will be written.
 * @param ctl A pointer to a `ctl_t` structure containing control parameters.
 * @param atm A pointer to an `atm_t` structure representing atmospheric data.
 * @param t A double value representing the time at which the VTK data is being written.
 *
 * The function performs the following steps:
 * - Sets a timer and logs information about writing VTK data.
 * - Sets a time interval for output based on the specified time and control parameters.
 * - Creates a new file and checks if the file creation was successful.
 * - Counts the number of data points to be written.
 * - Writes the VTK header, including metadata.
 * - Writes point coordinates based on the sphere or Cartesian coordinate system.
 * - Writes point data for each quantity specified in the control parameters.
 * - Closes the output file.
 *
 * @note This function writes VTK data in ASCII format, including point coordinates
 *       and associated scalar data for visualization purposes.
 *
 * @author Lars Hoffmann
 */
void write_vtk(
  const char *filename,
  ctl_t * ctl,
  atm_t * atm,
  double t);

/* ------------------------------------------------------------
   OpenACC routines...
   ------------------------------------------------------------ */

#ifdef _OPENACC
#pragma acc routine (clim_oh)
#pragma acc routine (clim_photo)
#pragma acc routine (clim_tropo)
#pragma acc routine (clim_ts)
#pragma acc routine (clim_zm)
#pragma acc routine (intpol_met_4d_coord)
#pragma acc routine (intpol_met_space_3d)
#pragma acc routine (intpol_met_space_3d_ml)
#pragma acc routine (intpol_met_space_2d)
#pragma acc routine (intpol_met_time_3d)
#pragma acc routine (intpol_met_time_3d_ml)
#pragma acc routine (intpol_met_time_2d)
#pragma acc routine (kernel_weight)
#pragma acc routine (lapse_rate)
#pragma acc routine (locate_irr)
#pragma acc routine (locate_irr_float)
#pragma acc routine (locate_reg)
#pragma acc routine (locate_vert)
#pragma acc routine (locate_irr_3d)
#pragma acc routine (nat_temperature)
#pragma acc routine (sedi)
#pragma acc routine (stddev)
#pragma acc routine (sza_calc)
#pragma acc routine (tropo_weight)
#endif

#endif /* LIBTRAC_H */
