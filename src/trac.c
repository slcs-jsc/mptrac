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

  Copyright (C) 2013-2018 Forschungszentrum Juelich GmbH
*/

/*!
  \file
  Lagrangian particle dispersion model.
*/

#include "libtrac.h"

#ifdef MPI
#include "mpi.h"
#endif

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Calculate advection of air parcels. */
void module_advection(
  met_t const * met0,
  met_t const * met1,
  atm_t * atm,
  int const ip,
  double const dt);

/*! Calculate exponential decay of particle mass. */
void module_decay(
  ctl_t const * ctl,
  met_t const * met0,
  met_t const * met1,
  atm_t * atm,
  int const ip,
  double const dt);

/*! Calculate mesoscale diffusion. */
void module_diffusion_meso(
  ctl_t const * ctl,
  met_t const * met0,
  met_t const * met1,
  atm_t * atm,
  int const ip,
  double const dt,
  gsl_rng * rng,
  double cache[EX][EY][EP][4]);

/*! Calculate turbulent diffusion. */
void module_diffusion_turb(
  ctl_t const * ctl,
  atm_t * atm,
  int const ip,
  double const dt,
  gsl_rng * rng);

/*! Force air parcels to stay on isosurface. */
void module_isosurf(
  ctl_t const * ctl,
  met_t const * met0,
  met_t const * met1,
  atm_t * atm,
  int const ip);

/*! Interpolate meteorological data for air parcel positions. */
void module_meteo(
  ctl_t const * ctl,
  met_t const * met0,
  met_t const * met1,
  atm_t * atm,
  int const ip);

/*! Check position of air parcels. */
void module_position(
  met_t const * met0,
  met_t const * met1,
  atm_t * atm,
  int const ip);

/*! Calculate sedimentation of air parcels. */
void module_sedi(
  ctl_t const * ctl,
  met_t const * met0,
  met_t const * met1,
  atm_t * atm,
  int const ip,
  double const dt);

/*! Write simulation output. */
void write_output(
  char const * dirname,
  ctl_t const * ctl,
  met_t const * met0,
  met_t const * met1,
  atm_t * atm,
  double const t,
  update_atm_meteo_data_t const key);

/*! Check if module_meteo needs to be executed at all */
update_atm_meteo_data_t need_meteo_update_at_all(
  ctl_t const * ctl);

/*! Check if module_meteo needs to be executed in this time iteration */
update_atm_meteo_data_t need_meteo_update_now(
  ctl_t const * ctl,
  update_atm_meteo_data_t const static_key,
  double const time_now);

  double accumulate_8winds(
    float const wind[EX][EY][EP],
    int const ix,
    int const iy,
    int const iz,
    double * sum_of_squares);

double get_standard_deviation(
  int const n,
  double const sum,
  double const sum_squares);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int const argc,
  char const *argv[]) {

  ctl_t ctl;

  atm_t *atm;

  met_t *met0, *met1;

  gsl_rng *rng[NTHREADS];

  FILE *dirlist;

  char dirname[LEN], filename[2 * LEN];

  double *dt, t;

  int i, ip, ntask = -1, rank = 0, size = 1;



  update_atm_meteo_data_t update_meteo_at_all = 0;
  update_atm_meteo_data_t update_meteo_now = 0;
  int do_turbulent_diffusion = 0;
  int do_mesoscale_diffusion = 0;
  int do_sedimentation = 0;
  int do_decay = 0, do_isosurface = 0;

  double (*meso_cache)[EY][EP][4]; /* allocate EX of these */
  int ix, iy, iz, i4;

  ALLOC(meso_cache, double, EX*EY*EP*4);
  /* mark all cache entries as not computed yet */
  for (ix = 0; ix < EX; ++ix) {
    for (iy = 0; iy < EY; ++iy) {
      for (iz = 0; iz < EP; ++iz) {
        for (i4 = 0; i4 < 4; ++i4) {
          meso_cache[ix][iy][iz][i4] = GSL_NAN;
        }
      }
    }
  }


#ifdef MPI
  /* Initialize MPI... */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <dirlist> <ctl> <atm_in> <metbase>");

  /* Open directory list... */
  if (!(dirlist = fopen(argv[1], "r")))
    ERRMSG("Cannot open directory list!");

  /* Loop over directories... */
  while (fscanf(dirlist, "%s", dirname) != EOF) {

    /* MPI parallelization... */
    if ((++ntask) % size != rank)
      continue;

    /* ------------------------------------------------------------
       Initialize model run...
       ------------------------------------------------------------ */

    /* Set timers... */
    START_TIMER(TIMER_TOTAL);
    START_TIMER(TIMER_INIT);

    /* Allocate... */
    ALLOC(atm, atm_t, 1);
    ALLOC(met0, met_t, 1);
    ALLOC(met1, met_t, 1);
    ALLOC(dt, double, NP);

    /* Initialize random number generators... */
    gsl_rng_env_setup();
    if (omp_get_max_threads() > NTHREADS)
      ERRMSG("Too many threads!");
    for (i = 0; i < NTHREADS; i++) {
      rng[i] = gsl_rng_alloc(gsl_rng_default);
      gsl_rng_set(rng[i], gsl_rng_default_seed + (long unsigned) i);
    }

    /* Read control parameters... */
    sprintf(filename, "%s/%s", dirname, argv[2]);
    read_ctl(filename, argc, argv, &ctl);

    /* Read atmospheric data... */
    sprintf(filename, "%s/%s", dirname, argv[3]);
    read_atm(filename, &ctl, atm);

    /* Set start time... */
    if (ctl.direction == 1) {
      ctl.t_start = gsl_stats_min(atm->time, 1, (size_t) atm->np);
      if (ctl.t_stop > 1e99)
	ctl.t_stop = gsl_stats_max(atm->time, 1, (size_t) atm->np);
    } else {
      ctl.t_start = gsl_stats_max(atm->time, 1, (size_t) atm->np);
      if (ctl.t_stop > 1e99)
	ctl.t_stop = gsl_stats_min(atm->time, 1, (size_t) atm->np);
    }

    /* Check time interval... */
    if (ctl.direction * (ctl.t_stop - ctl.t_start) <= 0)
      ERRMSG("Nothing to do!");

    /* Round start time... */
    if (ctl.direction == 1)
      ctl.t_start = floor(ctl.t_start / ctl.dt_mod) * ctl.dt_mod;
    else
      ctl.t_start = ceil(ctl.t_start / ctl.dt_mod) * ctl.dt_mod;

    /* Check if output in general is wanted */
    update_meteo_at_all = need_meteo_update_at_all(&ctl);

    do_turbulent_diffusion = (ctl.turb_dx_trop > 0 ||
      ctl.turb_dz_trop > 0 || ctl.turb_dx_strat > 0 || ctl.turb_dz_strat > 0);
    do_mesoscale_diffusion = (ctl.turb_mesox > 0 || ctl.turb_mesoz > 0);
    do_sedimentation = (ctl.qnt_r >= 0 && ctl.qnt_rho >= 0);
    do_decay = ((ctl.tdec_trop > 0 || ctl.tdec_strat > 0) && ctl.qnt_m >= 0);
    do_isosurface = (ctl.isosurf >= 1 && ctl.isosurf <= 4);

    /* Set timers... */
    STOP_TIMER(TIMER_INIT);

    /* ------------------------------------------------------------
       Loop over timesteps...
       ------------------------------------------------------------ */

    /* Loop over timesteps... */
    for (t = ctl.t_start; ctl.direction * (t - ctl.t_stop) < ctl.dt_mod;
         t += ctl.direction * ctl.dt_mod) {

      /* Adjust length of final time step... */
      if (ctl.direction * (t - ctl.t_stop) > 0) {
        t = ctl.t_stop;
      }
      
      update_meteo_now = need_meteo_update_now(&ctl, update_meteo_at_all, t);

      /* Set time steps for air parcels... */
#pragma omp parallel for default(shared) private(ip)
      for (ip = 0; ip < atm->np; ip++) {
        dt[ip] = ((ctl.direction * (atm->time[ip] - ctl.t_start) >= 0
                && ctl.direction * (atm->time[ip] - ctl.t_stop) <= 0
                && ctl.direction * (atm->time[ip] - t) < 0))
                ? (t - atm->time[ip]) : (GSL_NAN);
      }

      /* Get meteorological data... */
      START_TIMER(TIMER_INPUT);
      get_met(&ctl, argv[4], t, met0, met1);
      if (ctl.dt_mod > fabs(met0->lon[1] - met0->lon[0]) * 111132. / 150.)
        printf("Warning: Violation of CFL criterion! Set DT_MOD <= %g s!\n",
                fabs(met0->lon[1] - met0->lon[0]) * 111132. / 150.);
      STOP_TIMER(TIMER_INPUT);

      /* Initialize isosurface... */
      START_TIMER(TIMER_ISOSURF);
      if (ctl.isosurf >= 1 && ctl.isosurf <= 4 && t == ctl.t_start) {
        module_isosurf(&ctl, met0, met1, atm, -1);
      }
      STOP_TIMER(TIMER_ISOSURF);

      /* Advection... */
      START_TIMER(TIMER_ADVECT);
#pragma omp parallel for default(shared) private(ip)
      for (ip = 0; ip < atm->np; ip++) {
        if (gsl_finite(dt[ip])) {
          module_advection(met0, met1, atm, ip, dt[ip]);
        }
      }
      STOP_TIMER(TIMER_ADVECT);

      /* Turbulent diffusion... */
      if (do_turbulent_diffusion) {
        START_TIMER(TIMER_DIFFTURB);
#pragma omp parallel for default(shared) private(ip)
        for (ip = 0; ip < atm->np; ip++) {
          if (gsl_finite(dt[ip])) {
            module_diffusion_turb(&ctl, atm, ip, dt[ip], rng[omp_get_thread_num()]);
          }
        }
        STOP_TIMER(TIMER_DIFFTURB);
      }

      /* Mesoscale diffusion... */
      if (do_mesoscale_diffusion) {
        START_TIMER(TIMER_DIFFMESO);
#pragma omp parallel for default(shared) private(ip)
        for (ip = 0; ip < atm->np; ip++) {
          if (gsl_finite(dt[ip])) {
            module_diffusion_meso(&ctl, met0, met1, atm, ip, dt[ip],
              rng[omp_get_thread_num()], meso_cache);
          }
        }
        STOP_TIMER(TIMER_DIFFMESO);
      }

      /* Sedimentation... */
      if (do_sedimentation) {
        START_TIMER(TIMER_SEDI);
#pragma omp parallel for default(shared) private(ip)
        for (ip = 0; ip < atm->np; ip++) {
          if (gsl_finite(dt[ip])) {
            module_sedi(&ctl, met0, met1, atm, ip, dt[ip]);
          }
        }
        STOP_TIMER(TIMER_SEDI);
      }

      /* Isosurface... */
      if (do_isosurface) {
        START_TIMER(TIMER_ISOSURF);
#pragma omp parallel for default(shared) private(ip)
        for (ip = 0; ip < atm->np; ip++) {
          module_isosurf(&ctl, met0, met1, atm, ip);
        }
        STOP_TIMER(TIMER_ISOSURF);
      }

      /* Position... */
      START_TIMER(TIMER_POSITION);
#pragma omp parallel for default(shared) private(ip)
      for (ip = 0; ip < atm->np; ip++) {
        module_position(met0, met1, atm, ip);
      }
      STOP_TIMER(TIMER_POSITION);

      /* Meteorological data... */
      if (update_meteo_now) {
        START_TIMER(TIMER_METEO);
#pragma omp parallel for default(shared) private(ip)
        for (ip = 0; ip < atm->np; ip++) {
          module_meteo(&ctl, met0, met1, atm, ip);
        }
        STOP_TIMER(TIMER_METEO);
      }

      /* Decay... */
      if (do_decay) {
        START_TIMER(TIMER_DECAY);
#pragma omp parallel for default(shared) private(ip)
        for (ip = 0; ip < atm->np; ip++) {
          if (gsl_finite(dt[ip])) {
            module_decay(&ctl, met0, met1, atm, ip, dt[ip]);
          }
        }
        STOP_TIMER(TIMER_DECAY);
      }

      if (update_meteo_now) {
        /* Write output... */
        START_TIMER(TIMER_OUTPUT);
        write_output(dirname, &ctl, met0, met1, atm, t, update_meteo_now);
        STOP_TIMER(TIMER_OUTPUT);
      }

    } /* end of Loop over timesteps */

    /* ------------------------------------------------------------
       Finalize model run...
       ------------------------------------------------------------ */

    /* Report memory usage... */
    printf("MEMORY_ATM = %g MByte\n", sizeof(atm_t) / 1024. / 1024.);
    printf("MEMORY_METEO = %g MByte\n", 2. * sizeof(met_t) / 1024. / 1024.);
    printf("MEMORY_DYNAMIC = %g MByte\n",
	   4 * NP * sizeof(double) / 1024. / 1024.);
    printf("MEMORY_STATIC = %g MByte\n",
	   ((3 * GX * GY + 4 * GX * GY * GZ) * sizeof(double)
	    + (EX * EY + EX * EY * EP) * sizeof(float)
	    + (GX * GY + GX * GY * GZ) * sizeof(int)) / 1024. / 1024.);

    /* Report problem size... */
    printf("SIZE_NP = %d\n", atm->np);
    printf("SIZE_TASKS = %d\n", size);
    printf("SIZE_THREADS = %d\n", omp_get_max_threads());

    /* Report timers... */
    STOP_TIMER(TIMER_TOTAL);
    PRINT_TIMER(TIMER_TOTAL);
    PRINT_TIMER(TIMER_INIT);
    PRINT_TIMER(TIMER_INPUT);
    PRINT_TIMER(TIMER_OUTPUT);
    PRINT_TIMER(TIMER_ADVECT);
    PRINT_TIMER(TIMER_DECAY);
    PRINT_TIMER(TIMER_DIFFMESO);
    PRINT_TIMER(TIMER_DIFFTURB);
    PRINT_TIMER(TIMER_ISOSURF);
    PRINT_TIMER(TIMER_METEO);
    PRINT_TIMER(TIMER_POSITION);
    PRINT_TIMER(TIMER_SEDI);

    /* Free random number generators... */
    for (i = 0; i < NTHREADS; i++)
      gsl_rng_free(rng[i]);

    /* Free... */
    free(atm);
    free(met0);
    free(met1);
    free(dt);
  }

#ifdef MPI
  /* Finalize MPI... */
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void module_advection(
  met_t const * met0,
  met_t const * met1,
  atm_t * atm,
  int const ip,
  double const dt) {

  double v[3], xm[4], xo[4];

  /* load old positions */
  xo[0] = atm->lon[ip];
  xo[1] = atm->lat[ip];
  xo[2] = atm->p[ip];
  xo[3] = atm->time[ip];
  
  /* Interpolate meteorological data... */
  intpol_met_time(met0, met1, xo[3], xo[2], xo[0], xo[1], 
    NULL, NULL, NULL, NULL, &v[0], &v[1], &v[2], NULL, NULL, NULL);

  /* Get position of the mid point... */
  xm[0] = xo[0] + DX2DEG(0.5 * dt * v[0] / 1000., xo[1]);
  xm[1] = xo[1] + DY2DEG(0.5 * dt * v[1] / 1000.);
  xm[2] = xo[2] + 0.5 * dt * v[2];
  xm[3] = xo[3] + 0.5 * dt;

  /* Interpolate meteorological data for mid point... */
  intpol_met_time(met0, met1, xm[3], xm[2], xm[0], xm[1], 
    NULL, NULL, NULL, NULL, &v[0], &v[1], &v[2], NULL, NULL, NULL);

  /* Save new position... */
  atm->lon[ip]  = xo[0] + DX2DEG(dt * v[0] / 1000., xm[1]);
  atm->lat[ip]  = xo[1] + DY2DEG(dt * v[1] / 1000.);
  atm->p[ip]    = xo[2] + dt * v[2];
  atm->time[ip] = xo[3] + dt;
}

/*****************************************************************************/

void module_decay(
  ctl_t const * ctl,
  met_t const * met0,
  met_t const * met1,
  atm_t * atm,
  int const ip,
  double const dt) {

  double ps, pt, tdec;

  /* Set constant lifetime... */
  if (ctl->tdec_trop == ctl->tdec_strat)
    tdec = ctl->tdec_trop;

  /* Set altitude-dependent lifetime... */
  else {

    /* Get surface pressure... */
    intpol_met_time(met0, met1, atm->time[ip], atm->p[ip],
		    atm->lon[ip], atm->lat[ip], &ps, NULL, NULL, NULL,
		    NULL, NULL, NULL, NULL, NULL, NULL);

    /* Get tropopause pressure... */
    pt = clim_tropo(atm->time[ip], atm->lat[ip]);

    /* Set lifetime... */
    if (atm->p[ip] <= pt)
      tdec = ctl->tdec_strat;
    else
      tdec = LIN(ps, ctl->tdec_trop, pt, ctl->tdec_strat, atm->p[ip]);
  }

  /* Calculate exponential decay... */
  atm->q[ctl->qnt_m][ip] *= exp(-dt / tdec);
}

/*****************************************************************************/


double get_standard_deviation(int const n, double const sum, double const sum_squares) {
  if (n < 2) return 0;
  return sqrt(GSL_MAX(0.0, n*sum_squares - sum*sum)/((double)n*((double)n - 1.0)));
}

double accumulate_8winds(float const wind[EX][EY][EP],
  int const ix, int const iy, int const iz, double * sum_of_squares) {
  double s, s2;
  s = 0; s2 = 0;
#define accumulate_stats(sum, sum_squares, new) \
                                    sum += (new); \
                            sum_squares += (new)*(new)
  accumulate_stats(s, s2, wind[ix + 0][iy + 0][iz + 0]);
  accumulate_stats(s, s2, wind[ix + 0][iy + 0][iz + 1]);
  accumulate_stats(s, s2, wind[ix + 0][iy + 1][iz + 0]);
  accumulate_stats(s, s2, wind[ix + 0][iy + 1][iz + 1]);
  accumulate_stats(s, s2, wind[ix + 1][iy + 0][iz + 0]);
  accumulate_stats(s, s2, wind[ix + 1][iy + 0][iz + 1]);
  accumulate_stats(s, s2, wind[ix + 1][iy + 1][iz + 0]);
  accumulate_stats(s, s2, wind[ix + 1][iy + 1][iz + 1]);
#undef accumulate_stats
  *sum_of_squares += s2;
  return s;
}

void module_diffusion_meso(
  ctl_t const * ctl,
  met_t const * met0,
  met_t const * met1,
  atm_t * atm,
  int const ip,
  double const dt,
  gsl_rng * rng,
  double cache[EX][EY][EP][4]) { /* right-most dim 4 contains {tref,usig,vsig,wsig} */

#ifndef EIGHT_WINDS
  double u[16], v[16], w[16];
#else
  double s, ss; /* new method, results differ silghtly */
#endif
  double reference_time, r, rs, usig, vsig, wsig;

  int ix, iy, iz, iz0, iz1;

  /* Get indices... */
  ix = locate_reg(met0->lon, met0->nx, atm->lon[ip]);
  iy = locate_reg(met0->lat, met0->ny, atm->lat[ip]);
  iz0 = locate_irr(met0->p, met0->np, atm->p[ip]);
  iz1 = locate_irr(met1->p, met1->np, atm->p[ip]); /* we can omit this extra line if we can be sure that the pressure grid is also constant over time */

  /* assume that met1 has the same grids! */
/*
  assert(locate(met1->lon, met1->nx, atm->lon[ip]) == ix);
  assert(locate(met1->lat, met1->ny, atm->lat[ip]) == iy);
  assert(locate(met1->p,   met1->np, atm->p[ip])   == iz);
*/

  /* take the reference time of any of the two,
  but always take either met0 or met1 */
reference_time = met0->time;

iz = iz0; /* check it this is a good choice */

if (reference_time == cache[ix][iy][iz][0]) { /* yes, we do test doubles for exact equality here */

  /* cache hit, load data from cache cell */
  usig = cache[ix][iy][iz][1];
  vsig = cache[ix][iy][iz][2];
  wsig = cache[ix][iy][iz][3];

} else {
  /* recompute sigmas and cache them */

  /* Collect local wind data... */
#ifndef EIGHT_WINDS
  iz = iz0;
  u[0] = met0->u[ix][iy][iz];
  u[1] = met0->u[ix + 1][iy][iz];
  u[2] = met0->u[ix][iy + 1][iz];
  u[3] = met0->u[ix + 1][iy + 1][iz];
  u[4] = met0->u[ix][iy][iz + 1];
  u[5] = met0->u[ix + 1][iy][iz + 1];
  u[6] = met0->u[ix][iy + 1][iz + 1];
  u[7] = met0->u[ix + 1][iy + 1][iz + 1];

  iz = iz1;
  u[8] = met1->u[ix][iy][iz];
  u[9] = met1->u[ix + 1][iy][iz];
  u[10] = met1->u[ix][iy + 1][iz];
  u[11] = met1->u[ix + 1][iy + 1][iz];
  u[12] = met1->u[ix][iy][iz + 1];
  u[13] = met1->u[ix + 1][iy][iz + 1];
  u[14] = met1->u[ix][iy + 1][iz + 1];
  u[15] = met1->u[ix + 1][iy + 1][iz + 1];

  /* Get standard deviations of local wind data... */
  usig = gsl_stats_sd(u, 1, 16);
#else
    ss = 0;
    s = accumulate_8winds(met0->u, ix, iy, iz, &ss)
      + accumulate_8winds(met1->u, ix, iy, iz, &ss);
    usig = get_standard_deviation(16, s, ss);
#endif

#ifndef EIGHT_WINDS
  iz = iz0;
  v[0] = met0->v[ix][iy][iz];
  v[1] = met0->v[ix + 1][iy][iz];
  v[2] = met0->v[ix][iy + 1][iz];
  v[3] = met0->v[ix + 1][iy + 1][iz];
  v[4] = met0->v[ix][iy][iz + 1];
  v[5] = met0->v[ix + 1][iy][iz + 1];
  v[6] = met0->v[ix][iy + 1][iz + 1];
  v[7] = met0->v[ix + 1][iy + 1][iz + 1];

  iz = iz1;
  v[8] = met1->v[ix][iy][iz];
  v[9] = met1->v[ix + 1][iy][iz];
  v[10] = met1->v[ix][iy + 1][iz];
  v[11] = met1->v[ix + 1][iy + 1][iz];
  v[12] = met1->v[ix][iy][iz + 1];
  v[13] = met1->v[ix + 1][iy][iz + 1];
  v[14] = met1->v[ix][iy + 1][iz + 1];
  v[15] = met1->v[ix + 1][iy + 1][iz + 1];

  /* Get standard deviations of local wind data... */
  vsig = gsl_stats_sd(v, 1, 16);
#else
    ss = 0;
    s = accumulate_8winds(met0->v, ix, iy, iz, &ss)
      + accumulate_8winds(met1->v, ix, iy, iz, &ss);
    vsig = get_standard_deviation(16, s, ss);
#endif

#ifndef EIGHT_WINDS
  iz = iz0;
  w[0] = met0->w[ix][iy][iz];
  w[1] = met0->w[ix + 1][iy][iz];
  w[2] = met0->w[ix][iy + 1][iz];
  w[3] = met0->w[ix + 1][iy + 1][iz];
  w[4] = met0->w[ix][iy][iz + 1];
  w[5] = met0->w[ix + 1][iy][iz + 1];
  w[6] = met0->w[ix][iy + 1][iz + 1];
  w[7] = met0->w[ix + 1][iy + 1][iz + 1];

  iz = iz1;
  w[8] = met1->w[ix][iy][iz];
  w[9] = met1->w[ix + 1][iy][iz];
  w[10] = met1->w[ix][iy + 1][iz];
  w[11] = met1->w[ix + 1][iy + 1][iz];
  w[12] = met1->w[ix][iy][iz + 1];
  w[13] = met1->w[ix + 1][iy][iz + 1];
  w[14] = met1->w[ix][iy + 1][iz + 1];
  w[15] = met1->w[ix + 1][iy + 1][iz + 1];

  /* Get standard deviations of local wind data... */
  wsig = gsl_stats_sd(w, 1, 16);
#else
    ss = 0;
    s = accumulate_8winds(met0->w, ix, iy, iz, &ss)
      + accumulate_8winds(met1->w, ix, iy, iz, &ss);
    wsig = get_standard_deviation(16, s, ss);
#endif

  /* store the results in the cache, update the cache time */
  iz = iz0;
  cache[ix][iy][iz][0] = reference_time;
  cache[ix][iy][iz][1] = usig;
  cache[ix][iy][iz][2] = vsig;
  cache[ix][iy][iz][3] = wsig;

}

  /* Set temporal correlations for mesoscale fluctuations... */
  r = 1 - 2 * fabs(dt) / ctl->dt_met;
  rs = sqrt(1 - r * r);

  /* Calculate horizontal mesoscale wind fluctuations... */
  if (ctl->turb_mesox > 0) {
    atm->up[ip] = (float)(r * atm->up[ip]
       + rs * gsl_ran_gaussian_ziggurat(rng, ctl->turb_mesox * usig));
    atm->lon[ip] += DX2DEG(atm->up[ip] * dt / 1000., atm->lat[ip]);

    atm->vp[ip] = (float)(r * atm->vp[ip]
       + rs * gsl_ran_gaussian_ziggurat(rng, ctl->turb_mesox * vsig));
    atm->lat[ip] += DY2DEG(atm->vp[ip] * dt / 1000.);
  }

  /* Calculate vertical mesoscale wind fluctuations... */
  if (ctl->turb_mesoz > 0) {
    atm->wp[ip] = (float)(r * atm->wp[ip]
       + rs * gsl_ran_gaussian_ziggurat(rng, ctl->turb_mesoz * wsig));
    atm->p[ip] += atm->wp[ip] * dt;
  }
}

/*****************************************************************************/

void module_diffusion_turb(
  ctl_t const * ctl,
  atm_t * atm,
  int const ip,
  double const dt,
  gsl_rng * rng) {

  double dx, dz, pt, p0, p1, w;

  /* Get tropopause pressure... */
  pt = clim_tropo(atm->time[ip], atm->lat[ip]);

  /* Get weighting factor... */
  p1 = pt * 0.866877899;
  p0 = pt / 0.866877899;
  if (atm->p[ip] > p0)
    w = 1;
  else if (atm->p[ip] < p1)
    w = 0;
  else
    w = LIN(p0, 1.0, p1, 0.0, atm->p[ip]);

  /* Set diffusivity... */
  dx = w * ctl->turb_dx_trop + (1 - w) * ctl->turb_dx_strat;
  dz = w * ctl->turb_dz_trop + (1 - w) * ctl->turb_dz_strat;

  /* Horizontal turbulent diffusion... */
  if (dx > 0) {
    atm->lon[ip]
      += DX2DEG(gsl_ran_gaussian_ziggurat(rng, sqrt(2.0 * dx * fabs(dt)))
		/ 1000., atm->lat[ip]);
    atm->lat[ip]
      += DY2DEG(gsl_ran_gaussian_ziggurat(rng, sqrt(2.0 * dx * fabs(dt)))
		/ 1000.);
  }

  /* Vertical turbulent diffusion... */
  if (dz > 0)
    atm->p[ip]
      += DZ2DP(gsl_ran_gaussian_ziggurat(rng, sqrt(2.0 * dz * fabs(dt)))
	       / 1000., atm->p[ip]);

}

/*****************************************************************************/

void module_isosurf(
  ctl_t const * ctl,
  met_t const * met0,
  met_t const * met1,
  atm_t * atm,
  int const ip) {

  static double *iso, *ps, t, *ts;

  static int idx, ip2, n;

  FILE *in;

  char line[LEN];

  /* Initialize... */
  if (ip < 0) {

    /* Allocate... */
    ALLOC(iso, double,
	  NP);
    ALLOC(ps, double,
	  NP);
    ALLOC(ts, double,
	  NP);

    /* Save pressure... */
    if (ctl->isosurf == 1)
      for (ip2 = 0; ip2 < atm->np; ip2++)
	iso[ip2] = atm->p[ip2];

    /* Save density... */
    else if (ctl->isosurf == 2)
      for (ip2 = 0; ip2 < atm->np; ip2++) {
	intpol_met_time(met0, met1, atm->time[ip2], atm->p[ip2],
			atm->lon[ip2], atm->lat[ip2], NULL, NULL, NULL,
			&t, NULL, NULL, NULL, NULL, NULL, NULL);
	iso[ip2] = atm->p[ip2] / t;
      }

    /* Save potential temperature... */
    else if (ctl->isosurf == 3)
      for (ip2 = 0; ip2 < atm->np; ip2++) {
	intpol_met_time(met0, met1, atm->time[ip2], atm->p[ip2],
			atm->lon[ip2], atm->lat[ip2], NULL, NULL, NULL,
			&t, NULL, NULL, NULL, NULL, NULL, NULL);
	iso[ip2] = THETA(atm->p[ip2], t);
      }

    /* Read balloon pressure data... */
    else if (ctl->isosurf == 4) {

      /* Write info... */
      printf("Read balloon pressure data: %s\n", ctl->balloon);

      /* Open file... */
      if (!(in = fopen(ctl->balloon, "r")))
	ERRMSG("Cannot open file!");

      /* Read pressure time series... */
      while (fgets(line, LEN, in))
	if (sscanf(line, "%lg %lg", &ts[n], &ps[n]) == 2)
	  if ((++n) > NP)
	    ERRMSG("Too many data points!");

      /* Check number of points... */
      if (n < 1)
	ERRMSG("Could not read any data!");

      /* Close file... */
      fclose(in);
    }

    /* Leave initialization... */
    return;
  }

  /* Restore pressure... */
  if (ctl->isosurf == 1)
    atm->p[ip] = iso[ip];

  /* Restore density... */
  else if (ctl->isosurf == 2) {
    intpol_met_time(met0, met1, atm->time[ip], atm->p[ip], atm->lon[ip],
		    atm->lat[ip], NULL, NULL, NULL, &t,
		    NULL, NULL, NULL, NULL, NULL, NULL);
    atm->p[ip] = iso[ip] * t;
  }

  /* Restore potential temperature... */
  else if (ctl->isosurf == 3) {
    intpol_met_time(met0, met1, atm->time[ip], atm->p[ip], atm->lon[ip],
		    atm->lat[ip], NULL, NULL, NULL, &t,
		    NULL, NULL, NULL, NULL, NULL, NULL);
    atm->p[ip] = 1000. * pow(iso[ip] / t, -1. / 0.286);
  }

  /* Interpolate pressure... */
  else if (ctl->isosurf == 4) {
    if (atm->time[ip] <= ts[0])
      atm->p[ip] = ps[0];
    else if (atm->time[ip] >= ts[n - 1])
      atm->p[ip] = ps[n - 1];
    else {
      idx = locate_irr(ts, n, atm->time[ip]);
      atm->p[ip] = LIN(ts[idx], ps[idx],
		       ts[idx + 1], ps[idx + 1], atm->time[ip]);
    }
  }
}

/*****************************************************************************/

void module_meteo(
  ctl_t const * ctl,
  met_t const * met0,
  met_t const * met1,
  atm_t * atm,
  int const ip) {

  double a, b, c, ps, pt, pv, p_hno3, p_h2o, t, u, v, w, x1, x2, h2o, o3, z;

  /* Interpolate meteorological data... */
  intpol_met_time(met0, met1, atm->time[ip], atm->p[ip], atm->lon[ip],
		  atm->lat[ip], &ps, &pt, &z, &t, &u, &v, &w, &pv, &h2o, &o3);

  /* Set surface pressure... */
  if (ctl->qnt_ps >= 0)
    atm->q[ctl->qnt_ps][ip] = ps;

  /* Set tropopause pressure... */
  if (ctl->qnt_pt >= 0)
    atm->q[ctl->qnt_pt][ip] = pt;

  /* Set pressure... */
  if (ctl->qnt_p >= 0)
    atm->q[ctl->qnt_p][ip] = atm->p[ip];

  /* Set geopotential height... */
  if (ctl->qnt_z >= 0)
    atm->q[ctl->qnt_z][ip] = z;

  /* Set temperature... */
  if (ctl->qnt_t >= 0)
    atm->q[ctl->qnt_t][ip] = t;

  /* Set zonal wind... */
  if (ctl->qnt_u >= 0)
    atm->q[ctl->qnt_u][ip] = u;

  /* Set meridional wind... */
  if (ctl->qnt_v >= 0)
    atm->q[ctl->qnt_v][ip] = v;

  /* Set vertical velocity... */
  if (ctl->qnt_w >= 0)
    atm->q[ctl->qnt_w][ip] = w;

  /* Set water vapor vmr... */
  if (ctl->qnt_h2o >= 0)
    atm->q[ctl->qnt_h2o][ip] = h2o;

  /* Set ozone vmr... */
  if (ctl->qnt_o3 >= 0)
    atm->q[ctl->qnt_o3][ip] = o3;

  /* Calculate horizontal wind... */
  if (ctl->qnt_vh >= 0)
    atm->q[ctl->qnt_vh][ip] = sqrt(u * u + v * v);

  /* Calculate vertical velocity... */
  if (ctl->qnt_vz >= 0)
    atm->q[ctl->qnt_vz][ip] = -1e3 * H0 / atm->p[ip] * w;

  /* Calculate potential temperature... */
  if (ctl->qnt_theta >= 0)
    atm->q[ctl->qnt_theta][ip] = THETA(atm->p[ip], t);

  /* Set potential vorticity... */
  if (ctl->qnt_pv >= 0)
    atm->q[ctl->qnt_pv][ip] = pv;

  /* Calculate T_ice (Marti and Mauersberger, 1993)... */
  if (ctl->qnt_tice >= 0)
    atm->q[ctl->qnt_tice][ip] =
      -2663.5 /
      (log10((ctl->psc_h2o > 0 ? ctl->psc_h2o : h2o) * atm->p[ip] * 100.) -
       12.537);

  /* Calculate T_NAT (Hanson and Mauersberger, 1988)... */
  if (ctl->qnt_tnat >= 0) {
    if (ctl->psc_hno3 > 0)
      p_hno3 = ctl->psc_hno3 * atm->p[ip] / 1.333224;
    else
      p_hno3 = clim_hno3(atm->time[ip], atm->lat[ip], atm->p[ip])
	* 1e-9 * atm->p[ip] / 1.333224;
    p_h2o = (ctl->psc_h2o > 0 ? ctl->psc_h2o : h2o) * atm->p[ip] / 1.333224;
    a = 0.009179 - 0.00088 * log10(p_h2o);
    b = (38.9855 - log10(p_hno3) - 2.7836 * log10(p_h2o)) / a;
    c = -11397.0 / a;
    x1 = (-b + sqrt(b * b - 4. * c)) / 2.;
    x2 = (-b - sqrt(b * b - 4. * c)) / 2.;
    if (x1 > 0)
      atm->q[ctl->qnt_tnat][ip] = x1;
    if (x2 > 0)
      atm->q[ctl->qnt_tnat][ip] = x2;
  }

  /* Calculate T_STS (mean of T_ice and T_NAT)... */
  if (ctl->qnt_tsts >= 0) {
    if (ctl->qnt_tice < 0 || ctl->qnt_tnat < 0)
      ERRMSG("Need T_ice and T_NAT to calculate T_STS!");
    atm->q[ctl->qnt_tsts][ip] = 0.5 * (atm->q[ctl->qnt_tice][ip]
				       + atm->q[ctl->qnt_tnat][ip]);
  }
}

/*****************************************************************************/

void module_position(
  met_t const * met0,
  met_t const * met1,
  atm_t * atm,
  int const ip) {

  double ps, lon, lon0, lat, lat0, pre, now;

  /* Load... */
  now = atm->time[ip];
  pre = atm->p[ip];
  lon = lon0 = atm->lon[ip];
  lat = lat0 = atm->lat[ip];
  
  /* Calculate modulo... */
  lon = fmod(lon, 360);
  lat = fmod(lat, 360);

  /* Check latitude... */
  while (lat < -90 || lat > 90) {
    if (lat > 90) {
      lat = 180 - lat;
      lon += 180;
    }
    if (lat < -90) {
      lat = -180 - lat;
      lon += 180;
    }
  }

  /* Store corrected latitude */
  if (lat != lat0)
    atm->lat[ip] = lat;

  /* Check longitude... */
  while (lon < -180)
    lon += 360;
  while (lon >= 180)
    lon -= 360;
  
  /* Get surface pressure... */
  intpol_met_time(met0, met1, now, pre, lon, lat, 
    &ps, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

  /* Store corrected longitude */
  if (lon != lon0)
    atm->lon[ip] = lon;
  
  /* Check pressure and store corrected pressure (maybe)... */
  if (pre > ps)
    atm->p[ip] = ps;
  else if (pre < met0->p[met0->np - 1])
    atm->p[ip] = met0->p[met0->np - 1];

}

/*****************************************************************************/

void module_sedi(
  ctl_t const * ctl,
  met_t const * met0,
  met_t const * met1,
  atm_t * atm,
  int const ip,
  double const dt) {

  /* Coefficients for Cunningham slip-flow correction (Kasten, 1968): */
  const double A = 1.249, B = 0.42, C = 0.87;

  /* Average mass of an air molecule [kg/molec]: */
  const double m = 4.8096e-26;

  double G, K, eta, lambda, p, r_p, rho, rho_p, T, v, v_p;

  /* Convert units... */
  p = 100 * atm->p[ip];
  r_p = 1e-6 * atm->q[ctl->qnt_r][ip];
  rho_p = atm->q[ctl->qnt_rho][ip];

  /* Get temperature... */
  intpol_met_time(met0, met1, atm->time[ip], atm->p[ip], atm->lon[ip],
		  atm->lat[ip], NULL, NULL, NULL, &T,
		  NULL, NULL, NULL, NULL, NULL, NULL);

  /* Density of dry air... */
  rho = p / (RA * T);

  /* Dynamic viscosity of air... */
  eta = 1.8325e-5 * (416.16 / (T + 120.)) * pow(T / 296.16, 1.5);

  /* Thermal velocity of an air molecule... */
  v = sqrt(8 * KB * T / (M_PI * m));

  /* Mean free path of an air molecule... */
  lambda = 2 * eta / (rho * v);

  /* Knudsen number for air... */
  K = lambda / r_p;

  /* Cunningham slip-flow correction... */
  G = 1 + K * (A + B * exp(-C / K));

  /* Sedimentation (fall) velocity... */
  v_p = 2. * gsl_pow_2(r_p) * (rho_p - rho) * G0 / (9. * eta) * G;

  /* Calculate pressure change... */
  atm->p[ip] += DZ2DP(v_p * dt / 1000., atm->p[ip]);
}

/*****************************************************************************/

void write_output(
  char const * dirname,
  ctl_t const * ctl,
  met_t const * met0,
  met_t const * met1,
  atm_t * atm, /* station flags are modified in write_station */
  double const t,
  update_atm_meteo_data_t const key) {

  char filename[2 * LEN];

  double r;

  int year, mon, day, hour, min, sec;

  /* Get time... */
  jsec2time(t, &year, &mon, &day, &hour, &min, &sec, &r);

  /* Write atmospheric data... */
  if (key & UPDATE_ATM_METEO_DATA_TIME_ATM) {
    sprintf(filename, "%s/%s_%04d_%02d_%02d_%02d_%02d.tab",
	    dirname, ctl->atm_basename, year, mon, day, hour, min);
    write_atm(filename, ctl, atm, t);
  }

  /* Write CSI data... */
  if (key & UPDATE_ATM_METEO_DATA_CSI) {
    sprintf(filename, "%s/%s.tab", dirname, ctl->csi_basename);
    write_csi(filename, ctl, atm, t);
  }

  /* Write ensemble data... */
  if (key & UPDATE_ATM_METEO_DATA_ENSEMBLE) {
    sprintf(filename, "%s/%s.tab", dirname, ctl->ens_basename);
    write_ens(filename, ctl, atm, t);
  }

  /* Write gridded data... */
  if (key & UPDATE_ATM_METEO_DATA_TIME_GRID) {
    sprintf(filename, "%s/%s_%04d_%02d_%02d_%02d_%02d.tab",
	    dirname, ctl->grid_basename, year, mon, day, hour, min);
    write_grid(filename, ctl, met0, met1, atm, t);
  }

  /* Write profile data... */
  if (key & UPDATE_ATM_METEO_DATA_PROFILE) {
    sprintf(filename, "%s/%s.tab", dirname, ctl->prof_basename);
    write_prof(filename, ctl, met0, met1, atm, t);
  }

  /* Write station data... */
  if (key & UPDATE_ATM_METEO_DATA_STATIONS) {
    sprintf(filename, "%s/%s.tab", dirname, ctl->stat_basename);
    write_station(filename, ctl, atm, t);
  }
}

/*****************************************************************************/

update_atm_meteo_data_t need_meteo_update_now(
  ctl_t const * ctl,
  update_atm_meteo_data_t const static_key,
  double const time_now) {
  update_atm_meteo_data_t key = static_key; /* get a non-const copy */

  if (key &    UPDATE_ATM_METEO_DATA_TIME_ATM)  /* bit flag for ATM is on */
      key &= ~(UPDATE_ATM_METEO_DATA_TIME_ATM   /* switch ATM bit flag off */
               *(0 != fmod(time_now, ctl->atm_dt_out))); /* if the modulo is non-zero */
  if (key &    UPDATE_ATM_METEO_DATA_TIME_GRID) /* bit flag for GRID is on */
      key &= ~(UPDATE_ATM_METEO_DATA_TIME_GRID  /* switch GRID bit flag off */
               *(0 != fmod(time_now, ctl->grid_dt_out))); /* if the modulo is non-zero */
  return key;
}

/*****************************************************************************/

update_atm_meteo_data_t need_meteo_update_at_all(
  ctl_t const * ctl) {
  return
    UPDATE_ATM_METEO_DATA_TIME_ATM  * ('-' != ctl->atm_basename[0])
  + UPDATE_ATM_METEO_DATA_TIME_GRID * ('-' != ctl->grid_basename[0])
  + UPDATE_ATM_METEO_DATA_CSI       * ('-' != ctl->csi_basename[0])
  + UPDATE_ATM_METEO_DATA_ENSEMBLE  * ('-' != ctl->ens_basename[0])
  + UPDATE_ATM_METEO_DATA_PROFILE   * ('-' != ctl->prof_basename[0])
  + UPDATE_ATM_METEO_DATA_STATIONS  * ('-' != ctl->stat_basename[0]);
}
