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
  Lagrangian particle dispersion model.
*/

#include "libtrac.h"

#ifdef MPI
#include "mpi.h"
#endif

/* ------------------------------------------------------------
   Defines...
   ------------------------------------------------------------ */

/*! Timer for total runtime. */
#define TIMER_TOTAL 0

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

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Set simulation time interval. */
void init_simtime(
  ctl_t * ctl,
  atm_t * atm);

/*! Calculate advection of air parcels. */
void module_advection(
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt);

/*! Calculate exponential decay of particle mass. */
void module_decay(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt);

/*! Calculate mesoscale diffusion. */
void module_diffusion_meso(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt,
  gsl_rng * rng);

/*! Calculate turbulent diffusion. */
void module_diffusion_turb(
  ctl_t * ctl,
  atm_t * atm,
  int ip,
  double dt,
  gsl_rng * rng);

/*! Force air parcels to stay on isosurface. */
void module_isosurf(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double *iso);

/*! Interpolate meteorological data for air parcel positions. */
void module_meteo(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip);

/*! Check position of air parcels. */
void module_position(
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip);

/*! Calculate sedimentation of air parcels. */
void module_sedi(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt);

/*! Write simulation output. */
void write_output(
  const char *dirname,
  ctl_t * ctl,
  atm_t * atm,
  double t);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  atm_t *atm;

  met_t *met0, *met1;

  gsl_rng *rng[NTHREADS];

  FILE *dirlist;

  char dirname[LEN], filename[LEN];

  double *dt, *iso, t, t0;

  int i, ip, ntask = 0, rank = 0, size = 1;

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
    ALLOC(dt, double,
	  NP);
    ALLOC(iso, double,
	  NP);

    /* Read control parameters... */
    sprintf(filename, "%s/%s", dirname, argv[2]);
    read_ctl(filename, argc, argv, &ctl);

    /* Initialize random number generators... */
    gsl_rng_env_setup();
    for (i = 0; i < NTHREADS; i++)
      rng[i] = gsl_rng_alloc(gsl_rng_default);

    /* Read atmospheric data... */
    sprintf(filename, "%s/%s", dirname, argv[3]);
    read_atm(filename, atm, &ctl);

    /* Get simulation time interval... */
    init_simtime(&ctl, atm);

    /* Get rounded start time... */
    if (ctl.direction == 1)
      t0 = floor(ctl.t_start / ctl.dt_mod) * ctl.dt_mod;
    else
      t0 = ceil(ctl.t_start / ctl.dt_mod) * ctl.dt_mod;

    /* Set timers... */
    STOP_TIMER(TIMER_INIT);

    /* ------------------------------------------------------------
       Loop over timesteps...
       ------------------------------------------------------------ */

    /* Loop over timesteps... */
    for (t = t0; ctl.direction * (t - ctl.t_stop) < ctl.dt_mod;
	 t += ctl.direction * ctl.dt_mod) {

      /* Adjust length of final time step... */
      if (ctl.direction * (t - ctl.t_stop) > 0)
	t = ctl.t_stop;

      /* Set time steps for air parcels... */
      for (ip = 0; ip < atm->np; ip++)
	if ((ctl.direction * (atm->time[ip] - ctl.t_start) >= 0
	     && ctl.direction * (atm->time[ip] - ctl.t_stop) <= 0
	     && ctl.direction * (atm->time[ip] - t) < 0))
	  dt[ip] = t - atm->time[ip];
	else
	  dt[ip] = GSL_NAN;

      /* Get meteorological data... */
      START_TIMER(TIMER_INPUT);
      get_met(&ctl, argv[4], t, met0, met1);
      STOP_TIMER(TIMER_INPUT);

      /* Isosurface... */
      START_TIMER(TIMER_ISOSURF);
#pragma omp parallel for default(shared) private(ip)
      for (ip = 0; ip < atm->np; ip++)
	module_isosurf(&ctl, met0, met1, atm, ip, iso);
      STOP_TIMER(TIMER_ISOSURF);

      /* Advection... */
      START_TIMER(TIMER_ADVECT);
#pragma omp parallel for default(shared) private(ip)
      for (ip = 0; ip < atm->np; ip++)
	if (gsl_finite(dt[ip]))
	  module_advection(met0, met1, atm, ip, dt[ip]);
      STOP_TIMER(TIMER_ADVECT);

      /* Turbulent diffusion... */
      START_TIMER(TIMER_DIFFTURB);
#pragma omp parallel for default(shared) private(ip)
      for (ip = 0; ip < atm->np; ip++)
	if (gsl_finite(dt[ip]))
	  module_diffusion_turb(&ctl, atm, ip, dt[ip],
				rng[omp_get_thread_num()]);
      STOP_TIMER(TIMER_DIFFTURB);

      /* Mesoscale diffusion... */
      START_TIMER(TIMER_DIFFMESO);
#pragma omp parallel for default(shared) private(ip)
      for (ip = 0; ip < atm->np; ip++)
	if (gsl_finite(dt[ip]))
	  module_diffusion_meso(&ctl, met0, met1, atm, ip, dt[ip],
				rng[omp_get_thread_num()]);
      STOP_TIMER(TIMER_DIFFMESO);

      /* Sedimentation... */
      START_TIMER(TIMER_SEDI);
#pragma omp parallel for default(shared) private(ip)
      for (ip = 0; ip < atm->np; ip++)
	if (gsl_finite(dt[ip]))
	  module_sedi(&ctl, met0, met1, atm, ip, dt[ip]);
      STOP_TIMER(TIMER_SEDI);

      /* Isosurface... */
      START_TIMER(TIMER_ISOSURF);
#pragma omp parallel for default(shared) private(ip)
      for (ip = 0; ip < atm->np; ip++)
	module_isosurf(&ctl, met0, met1, atm, ip, iso);
      STOP_TIMER(TIMER_ISOSURF);

      /* Position... */
      START_TIMER(TIMER_POSITION);
#pragma omp parallel for default(shared) private(ip)
      for (ip = 0; ip < atm->np; ip++)
	module_position(met0, met1, atm, ip);
      STOP_TIMER(TIMER_POSITION);

      /* Meteorological data... */
      START_TIMER(TIMER_METEO);
#pragma omp parallel for default(shared) private(ip)
      for (ip = 0; ip < atm->np; ip++)
	module_meteo(&ctl, met0, met1, atm, ip);
      STOP_TIMER(TIMER_METEO);

      /* Decay... */
      START_TIMER(TIMER_DECAY);
#pragma omp parallel for default(shared) private(ip)
      for (ip = 0; ip < atm->np; ip++)
	if (gsl_finite(dt[ip]))
	  module_decay(&ctl, met0, met1, atm, ip, dt[ip]);
      STOP_TIMER(TIMER_DECAY);

      /* Write output... */
      START_TIMER(TIMER_OUTPUT);
      write_output(dirname, &ctl, atm, t);
      STOP_TIMER(TIMER_OUTPUT);
    }

    /* ------------------------------------------------------------
       Finalize model run...
       ------------------------------------------------------------ */

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

    /* Report memory usage... */
    printf("MEMORY_ATM = %g MByte\n", 2. * sizeof(atm_t) / 1024. / 1024.);
    printf("MEMORY_METEO = %g MByte\n", 2. * sizeof(met_t) / 1024. / 1024.);
    printf("MEMORY_DYNAMIC = %g MByte\n",
	   NP * sizeof(double) / 1024. / 1024.);
    printf("MEMORY_STATIC = %g MByte\n",
	   (((EX + EY) + (2 + NQ) * GX * GY * GZ) * sizeof(double)
	    + (EX * EY + EX * EY * EP) * sizeof(float)
	    + (2 * GX * GY * GZ) * sizeof(int)) / 1024. / 1024.);

    /* Report problem size... */
    printf("SIZE_NP = %d\n", atm->np);
    printf("SIZE_TASKS = %d\n", size);
    printf("SIZE_THREADS = %d\n", omp_get_max_threads());

    /* Free random number generators... */
    for (i = 0; i < NTHREADS; i++)
      gsl_rng_free(rng[i]);

    /* Free... */
    free(atm);
    free(met0);
    free(met1);
    free(dt);
    free(iso);
  }

#ifdef MPI
  /* Finalize MPI... */
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void init_simtime(
  ctl_t * ctl,
  atm_t * atm) {

  /* Set inital and final time... */
  if (ctl->direction == 1) {
    if (ctl->t_start < -1e99)
      ctl->t_start = gsl_stats_min(atm->time, 1, (size_t) atm->np);
    if (ctl->t_stop < -1e99)
      ctl->t_stop = gsl_stats_max(atm->time, 1, (size_t) atm->np);
  } else if (ctl->direction == -1) {
    if (ctl->t_stop < -1e99)
      ctl->t_stop = gsl_stats_min(atm->time, 1, (size_t) atm->np);
    if (ctl->t_start < -1e99)
      ctl->t_start = gsl_stats_max(atm->time, 1, (size_t) atm->np);
  }

  /* Check time... */
  if (ctl->direction * (ctl->t_stop - ctl->t_start) <= 0)
    ERRMSG("Nothing to do!");
}

/*****************************************************************************/

void module_advection(
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt) {

  double v[3], xm[3];

  /* Interpolate meteorological data... */
  intpol_met_time(met0, met1, atm->time[ip], atm->p[ip],
		  atm->lon[ip], atm->lat[ip], NULL, NULL,
		  &v[0], &v[1], &v[2], NULL, NULL);

  /* Get position of the mid point... */
  xm[0] = atm->lon[ip] + dx2deg(0.5 * dt * v[0] / 1000., atm->lat[ip]);
  xm[1] = atm->lat[ip] + dy2deg(0.5 * dt * v[1] / 1000.);
  xm[2] = atm->p[ip] + 0.5 * dt * v[2];

  /* Interpolate meteorological data for mid point... */
  intpol_met_time(met0, met1, atm->time[ip] + 0.5 * dt,
		  xm[2], xm[0], xm[1], NULL, NULL,
		  &v[0], &v[1], &v[2], NULL, NULL);

  /* Save new position... */
  atm->time[ip] += dt;
  atm->lon[ip] += dx2deg(dt * v[0] / 1000., xm[1]);
  atm->lat[ip] += dy2deg(dt * v[1] / 1000.);
  atm->p[ip] += dt * v[2];
}

/*****************************************************************************/

void module_decay(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt) {

  double ps, pt, tdec;

  /* Check lifetime values... */
  if ((ctl->tdec_trop <= 0 && ctl->tdec_strat <= 0) || ctl->qnt_m < 0)
    return;

  /* Set constant lifetime... */
  if (ctl->tdec_trop == ctl->tdec_strat)
    tdec = ctl->tdec_trop;

  /* Set altitude-dependent lifetime... */
  else {

    /* Get surface pressure... */
    intpol_met_time(met0, met1, atm->time[ip], atm->p[ip],
		    atm->lon[ip], atm->lat[ip], &ps, NULL,
		    NULL, NULL, NULL, NULL, NULL);

    /* Get tropopause pressure... */
    pt = tropopause(atm->time[ip], atm->lat[ip]);

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

void module_diffusion_meso(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt,
  gsl_rng * rng) {

  double r, rs, u[16], v[16], w[16], usig, vsig, wsig;

  int ix, iy, iz;

  /* Calculate mesoscale velocity fluctuations... */
  if (ctl->turb_meso > 0) {

    /* Get indices... */
    ix = locate(met0->lon, met0->nx, atm->lon[ip]);
    iy = locate(met0->lat, met0->ny, atm->lat[ip]);
    iz = locate(met0->p, met0->np, atm->p[ip]);

    /* Collect local wind data... */
    u[0] = met0->u[ix][iy][iz];
    u[1] = met0->u[ix + 1][iy][iz];
    u[2] = met0->u[ix][iy + 1][iz];
    u[3] = met0->u[ix + 1][iy + 1][iz];
    u[4] = met0->u[ix][iy][iz + 1];
    u[5] = met0->u[ix + 1][iy][iz + 1];
    u[6] = met0->u[ix][iy + 1][iz + 1];
    u[7] = met0->u[ix + 1][iy + 1][iz + 1];

    v[0] = met0->v[ix][iy][iz];
    v[1] = met0->v[ix + 1][iy][iz];
    v[2] = met0->v[ix][iy + 1][iz];
    v[3] = met0->v[ix + 1][iy + 1][iz];
    v[4] = met0->v[ix][iy][iz + 1];
    v[5] = met0->v[ix + 1][iy][iz + 1];
    v[6] = met0->v[ix][iy + 1][iz + 1];
    v[7] = met0->v[ix + 1][iy + 1][iz + 1];

    w[0] = met0->w[ix][iy][iz];
    w[1] = met0->w[ix + 1][iy][iz];
    w[2] = met0->w[ix][iy + 1][iz];
    w[3] = met0->w[ix + 1][iy + 1][iz];
    w[4] = met0->w[ix][iy][iz + 1];
    w[5] = met0->w[ix + 1][iy][iz + 1];
    w[6] = met0->w[ix][iy + 1][iz + 1];
    w[7] = met0->w[ix + 1][iy + 1][iz + 1];

    /* Get indices... */
    ix = locate(met1->lon, met1->nx, atm->lon[ip]);
    iy = locate(met1->lat, met1->ny, atm->lat[ip]);
    iz = locate(met1->p, met1->np, atm->p[ip]);

    /* Collect local wind data... */
    u[8] = met1->u[ix][iy][iz];
    u[9] = met1->u[ix + 1][iy][iz];
    u[10] = met1->u[ix][iy + 1][iz];
    u[11] = met1->u[ix + 1][iy + 1][iz];
    u[12] = met1->u[ix][iy][iz + 1];
    u[13] = met1->u[ix + 1][iy][iz + 1];
    u[14] = met1->u[ix][iy + 1][iz + 1];
    u[15] = met1->u[ix + 1][iy + 1][iz + 1];

    v[8] = met1->v[ix][iy][iz];
    v[9] = met1->v[ix + 1][iy][iz];
    v[10] = met1->v[ix][iy + 1][iz];
    v[11] = met1->v[ix + 1][iy + 1][iz];
    v[12] = met1->v[ix][iy][iz + 1];
    v[13] = met1->v[ix + 1][iy][iz + 1];
    v[14] = met1->v[ix][iy + 1][iz + 1];
    v[15] = met1->v[ix + 1][iy + 1][iz + 1];

    w[8] = met1->w[ix][iy][iz];
    w[9] = met1->w[ix + 1][iy][iz];
    w[10] = met1->w[ix][iy + 1][iz];
    w[11] = met1->w[ix + 1][iy + 1][iz];
    w[12] = met1->w[ix][iy][iz + 1];
    w[13] = met1->w[ix + 1][iy][iz + 1];
    w[14] = met1->w[ix][iy + 1][iz + 1];
    w[15] = met1->w[ix + 1][iy + 1][iz + 1];

    /* Get standard deviations of local wind data... */
    usig = gsl_stats_sd(u, 1, 16);
    vsig = gsl_stats_sd(v, 1, 16);
    wsig = gsl_stats_sd(w, 1, 16);

    /* Set temporal correlations for mesoscale fluctuations... */
    r = 1 - 2 * fabs(dt) / ctl->dt_met;
    rs = sqrt(1 - r * r);

    /* Calculate mesoscale wind fluctuations... */
    atm->up[ip] =
      r * atm->up[ip] + rs * gsl_ran_gaussian_ziggurat(rng,
						       ctl->turb_meso * usig);
    atm->vp[ip] =
      r * atm->vp[ip] + rs * gsl_ran_gaussian_ziggurat(rng,
						       ctl->turb_meso * vsig);
    atm->wp[ip] =
      r * atm->wp[ip] + rs * gsl_ran_gaussian_ziggurat(rng,
						       ctl->turb_meso * wsig);

    /* Calculate air parcel displacement... */
    atm->lon[ip] += dx2deg(atm->up[ip] * dt / 1000., atm->lat[ip]);
    atm->lat[ip] += dy2deg(atm->vp[ip] * dt / 1000.);
    atm->p[ip] += atm->wp[ip] * dt;
  }
}

/*****************************************************************************/

void module_diffusion_turb(
  ctl_t * ctl,
  atm_t * atm,
  int ip,
  double dt,
  gsl_rng * rng) {

  double dx, dz, pt, p0, p1, w;

  /* Get tropopause pressure... */
  pt = tropopause(atm->time[ip], atm->lat[ip]);

  /* Get weighting factor... */
  p1 = pt * 0.866877899;
  p0 = pt / 0.866877899;
  if (atm->p[ip] > p0)
    w = 1;
  else if (atm->p[ip] < p1)
    w = 0;
  else
    w = LIN(p0, 1.0, p1, 0.0, atm->p[ip]);

  /* Set diffusivitiy... */
  dx = w * ctl->turb_dx_trop + (1 - w) * ctl->turb_dx_strat;
  dz = w * ctl->turb_dz_trop + (1 - w) * ctl->turb_dz_strat;

  /* Horizontal turbulent diffusion... */
  if (dx > 0) {
    atm->lon[ip]
      += dx2deg(gsl_ran_gaussian_ziggurat(rng, sqrt(2.0 * dx * fabs(dt)))
		/ 1000., atm->lat[ip]);
    atm->lat[ip]
      += dy2deg(gsl_ran_gaussian_ziggurat(rng, sqrt(2.0 * dx * fabs(dt)))
		/ 1000.);
  }

  /* Vertical turbulent diffusion... */
  if (dz > 0)
    atm->p[ip]
      += dz2dp(gsl_ran_gaussian_ziggurat(rng, sqrt(2.0 * dz * fabs(dt)))
	       / 1000., atm->p[ip]);
}

/*****************************************************************************/

void module_isosurf(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double *iso) {

  double t;

  /* Check control parameter... */
  if (ctl->isosurf < 1 || ctl->isosurf > 3)
    return;

  /* Set initial values... */
  if (iso[ip] == 0) {

    /* Save pressure... */
    if (ctl->isosurf == 1)
      iso[ip] = atm->p[ip];

    /* Save density... */
    else if (ctl->isosurf == 2) {
      intpol_met_time(met0, met1, atm->time[ip], atm->p[ip], atm->lon[ip],
		      atm->lat[ip], NULL, &t, NULL, NULL, NULL, NULL, NULL);
      iso[ip] = atm->p[ip] / t;
    }

    /* Save potential temperature... */
    else if (ctl->isosurf == 3) {
      intpol_met_time(met0, met1, atm->time[ip], atm->p[ip], atm->lon[ip],
		      atm->lat[ip], NULL, &t, NULL, NULL, NULL, NULL, NULL);
      iso[ip] = t * pow(P0 / atm->p[ip], 0.286);
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
		    atm->lat[ip], NULL, &t, NULL, NULL, NULL, NULL, NULL);
    atm->p[ip] = iso[ip] * t;
  }

  /* Restore potential temperature... */
  else if (ctl->isosurf == 3) {
    intpol_met_time(met0, met1, atm->time[ip], atm->p[ip], atm->lon[ip],
		    atm->lat[ip], NULL, &t, NULL, NULL, NULL, NULL, NULL);
    atm->p[ip] = P0 * pow(iso[ip] / t, -1. / 0.286);
  }
}

/*****************************************************************************/

void module_meteo(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip) {

  /* Interpolate surface pressure... */
  if (ctl->qnt_ps >= 0)
    intpol_met_time(met0, met1, atm->time[ip], atm->p[ip], atm->lon[ip],
		    atm->lat[ip], &atm->q[ctl->qnt_ps][ip], NULL,
		    NULL, NULL, NULL, NULL, NULL);

  /* Interpolate temperature... */
  if (ctl->qnt_t >= 0)
    intpol_met_time(met0, met1, atm->time[ip], atm->p[ip], atm->lon[ip],
		    atm->lat[ip], NULL, &atm->q[ctl->qnt_t][ip],
		    NULL, NULL, NULL, NULL, NULL);

  /* Interpolate zonal wind... */
  if (ctl->qnt_u >= 0)
    intpol_met_time(met0, met1, atm->time[ip], atm->p[ip], atm->lon[ip],
		    atm->lat[ip], NULL, NULL, &atm->q[ctl->qnt_u][ip],
		    NULL, NULL, NULL, NULL);

  /* Interpolate meridional wind... */
  if (ctl->qnt_v >= 0)
    intpol_met_time(met0, met1, atm->time[ip], atm->p[ip], atm->lon[ip],
		    atm->lat[ip], NULL, NULL, NULL,
		    &atm->q[ctl->qnt_v][ip], NULL, NULL, NULL);

  /* Interpolate vertical velocity... */
  if (ctl->qnt_w >= 0)
    intpol_met_time(met0, met1, atm->time[ip], atm->p[ip], atm->lon[ip],
		    atm->lat[ip], NULL, NULL, NULL, NULL,
		    &atm->q[ctl->qnt_w][ip], NULL, NULL);

  /* Interpolate water vapor vmr... */
  if (ctl->qnt_h2o >= 0)
    intpol_met_time(met0, met1, atm->time[ip], atm->p[ip], atm->lon[ip],
		    atm->lat[ip], NULL, NULL, NULL, NULL, NULL,
		    &atm->q[ctl->qnt_h2o][ip], NULL);

  /* Interpolate ozone... */
  if (ctl->qnt_o3 >= 0)
    intpol_met_time(met0, met1, atm->time[ip], atm->p[ip], atm->lon[ip],
		    atm->lat[ip], NULL, NULL, NULL, NULL, NULL, NULL,
		    &atm->q[ctl->qnt_o3][ip]);

  /* Calculate potential temperature... */
  if (ctl->qnt_theta >= 0) {
    intpol_met_time(met0, met1, atm->time[ip], atm->p[ip], atm->lon[ip],
		    atm->lat[ip], NULL, &atm->q[ctl->qnt_theta][ip],
		    NULL, NULL, NULL, NULL, NULL);
    atm->q[ctl->qnt_theta][ip] *= pow(P0 / atm->p[ip], 0.286);
  }
}

/*****************************************************************************/

void module_position(
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip) {

  double ps;

  /* Calculate modulo... */
  atm->lon[ip] = fmod(atm->lon[ip], 360);
  atm->lat[ip] = fmod(atm->lat[ip], 360);

  /* Check latitude... */
  while (atm->lat[ip] < -90 || atm->lat[ip] > 90) {
    if (atm->lat[ip] > 90) {
      atm->lat[ip] = 180 - atm->lat[ip];
      atm->lon[ip] += 180;
    }
    if (atm->lat[ip] < -90) {
      atm->lat[ip] = -180 - atm->lat[ip];
      atm->lon[ip] += 180;
    }
  }

  /* Check longitude... */
  while (atm->lon[ip] < -180)
    atm->lon[ip] += 360;
  while (atm->lon[ip] >= 180)
    atm->lon[ip] -= 360;

  /* Get surface pressure... */
  intpol_met_time(met0, met1, atm->time[ip], atm->p[ip],
		  atm->lon[ip], atm->lat[ip], &ps, NULL,
		  NULL, NULL, NULL, NULL, NULL);

  /* Check pressure... */
  if (atm->p[ip] > ps)
    atm->p[ip] = ps;
  else if (atm->p[ip] < met0->p[met0->np - 1])
    atm->p[ip] = met0->p[met0->np - 1];
}

/*****************************************************************************/

void module_sedi(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt) {

  /* Coefficients for Cunningham slip-flow correction (Kasten, 1968): */
  const double A = 1.249, B = 0.42, C = 0.87;

  /* Specific gas constant for dry air [J/(kg K)]: */
  const double R = 287.058;

  /* Average mass of an air molecule [kg/molec]: */
  const double m = 4.8096e-26;

  double G, K, eta, lambda, p, r_p, rho, rho_p, T, v, v_p;

  /* Check if parameters are available... */
  if (ctl->qnt_r < 0 || ctl->qnt_rho < 0)
    return;

  /* Convert units... */
  p = 100 * atm->p[ip];
  r_p = 1e-6 * atm->q[ctl->qnt_r][ip];
  rho_p = atm->q[ctl->qnt_rho][ip];

  /* Get temperature... */
  intpol_met_time(met0, met1, atm->time[ip], atm->p[ip], atm->lon[ip],
		  atm->lat[ip], NULL, &T, NULL, NULL, NULL, NULL, NULL);

  /* Density of dry air... */
  rho = p / (R * T);

  /* Dynamic viscosity of air... */
  eta = 1.8325e-5 * (416.16 / (T + 120.)) * pow(T / 296.16, 1.5);

  /* Thermal velocity of an air molecule... */
  v = sqrt(8 * GSL_CONST_MKSA_BOLTZMANN * T / (M_PI * m));

  /* Mean free path of an air molecule... */
  lambda = 2 * eta / (rho * v);

  /* Knudsen number for air... */
  K = lambda / r_p;

  /* Cunningham slip-flow correction... */
  G = 1 + K * (A + B * exp(-C / K));

  /* Sedimentation (fall) velocity... */
  v_p =
    2. * gsl_pow_2(r_p) * (rho_p -
			   rho) * GSL_CONST_MKSA_GRAV_ACCEL / (9. * eta) * G;

  /* Calculate pressure change... */
  atm->p[ip] += dz2dp(v_p * dt / 1000., atm->p[ip]);
}

/*****************************************************************************/

void write_output(
  const char *dirname,
  ctl_t * ctl,
  atm_t * atm,
  double t) {

  char filename[LEN];

  double r;

  int year, mon, day, hour, min, sec;

  /* Get time... */
  jsec2time(t, &year, &mon, &day, &hour, &min, &sec, &r);

  /* Write atmospheric data... */
  if (ctl->atm_basename[0] != '-' && fmod(t, ctl->atm_dt_out) == 0) {
    sprintf(filename, "%s/%s_%04d_%02d_%02d_%02d_%02d.tab",
	    dirname, ctl->atm_basename, year, mon, day, hour, min);
    write_atm(filename, atm, ctl, t);
  }

  /* Write gridded data... */
  if (ctl->grid_basename[0] != '-' && fmod(t, ctl->grid_dt_out) == 0) {
    sprintf(filename, "%s/%s_%04d_%02d_%02d_%02d_%02d.tab",
	    dirname, ctl->grid_basename, year, mon, day, hour, min);
    write_grid(filename, atm, ctl, t);
  }

  /* Write CSI data... */
  if (ctl->csi_basename[0] != '-') {
    sprintf(filename, "%s/%s.tab", dirname, ctl->csi_basename);
    write_csi(filename, atm, ctl, t);
  }

  /* Write sample data... */
  if (ctl->sample_basename[0] != '-') {
    sprintf(filename, "%s/%s.tab", dirname, ctl->sample_basename);
    write_sample(filename, atm, ctl, t);
  }

  /* Write station data... */
  if (ctl->stat_basename[0] != '-') {
    sprintf(filename, "%s/%s.tab", dirname, ctl->stat_basename);
    write_station(filename, atm, ctl, t);
  }
}
