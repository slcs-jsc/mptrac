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
   Global variables...
   ------------------------------------------------------------ */

/*! Timer for total runtime. */
int timer_total;

/*! Timer for physics calculations. */
int timer_phys;

/*! Timer for file input. */
int timer_input;

/*! Timer for file output. */
int timer_output;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Set simulation time interval. */
void init_simtime(
  ctl_t * ctl,
  atm_t * atm);

/*! Calculate advection of air parcels. */
void module_advection(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt);

/*! Calculate exponential decay of particle mass. */
void module_decay(
  ctl_t * ctl,
  atm_t * atm,
  int ip,
  double dt);

/*! Calculate turbulent and mesoscale diffusion. */
void module_diffusion(
  ctl_t * ctl,
  met_t * met_t0,
  met_t * met_t1,
  atm_t * atm,
  int ip,
  double dt,
  gsl_rng * rng);

/*! Check position of air parcels. */
void module_position(
  met_t * met,
  atm_t * atm,
  int ip);

/*! Calculate sedimentation of air parcels. */
void module_sedi(
  ctl_t * ctl,
  atm_t * atm,
  int ip,
  double dt);

/*! Write simulation output. */
void write_output(
  const char *dirname,
  ctl_t * ctl,
  atm_t * atm,
  double t,
  int force);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  atm_t *atm;

  met_t *met_t0, *met_t1;

  gsl_rng *rng[NTHREADS];

  FILE *dirlist;

  char dirname[LEN];

  double dt, t, t0;

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

    /* Create timers... */
    timer_total = CREATE_TIMER("total");
    timer_phys = CREATE_TIMER("physics");
    timer_input = CREATE_TIMER("input");
    timer_output = CREATE_TIMER("output");

    /* Start timer for total runtime... */
    START_TIMER(timer_total);

    /* Allocate... */
    ALLOC(atm, atm_t, 1);
    ALLOC(met_t0, met_t, 1);
    ALLOC(met_t1, met_t, 1);

    /* Read control parameters... */
    read_ctl(dirname, argv[2], argc, argv, &ctl);

    /* Initialize random number generators... */
    gsl_rng_env_setup();
    for (i = 0; i < NTHREADS; i++)
      rng[i] = gsl_rng_alloc(gsl_rng_default);

    /* Read atmospheric data... */
    START_TIMER(timer_input);
    read_atm(dirname, argv[3], atm, &ctl);
    STOP_TIMER(timer_input);

    /* Get simulation time interval... */
    init_simtime(&ctl, atm);

    /* Write initial output... */
    START_TIMER(timer_output);
    write_output(dirname, &ctl, atm, ctl.t_start, 1);
    STOP_TIMER(timer_output);

    /* ------------------------------------------------------------
       Loop over timesteps...
       ------------------------------------------------------------ */

    /* Get rounded start time... */
    if (ctl.direction == 1)
      t0 = floor(ctl.t_start / ctl.dt_mod) * ctl.dt_mod;
    else
      t0 = ceil(ctl.t_start / ctl.dt_mod) * ctl.dt_mod;

    /* Loop over timesteps... */
    for (t = t0; ctl.direction * (t - ctl.t_stop) < ctl.dt_mod;
	 t += ctl.direction * ctl.dt_mod) {

      /* Adjust length of final time step... */
      if (ctl.direction * (t - ctl.t_stop) > 0)
	t = ctl.t_stop;

      /* Get meteorological data... */
      START_TIMER(timer_input);
      get_met(t, ctl.direction, argv[4], ctl.dt_met, ctl.red_met, met_t0,
	      met_t1);
      STOP_TIMER(timer_input);

      /* Loop over air parcels... */
      START_TIMER(timer_phys);
#pragma omp parallel for default(shared) private(dt,ip)
      for (ip = 0; ip < atm->np; ip++)
	if ((ctl.direction * (atm->time[ip] - ctl.t_start) >= 0
	     && ctl.direction * (atm->time[ip] - ctl.t_stop) <= 0
	     && ctl.direction * (atm->time[ip] - t) < 0)) {

	  /* Set time step... */
	  dt = t - atm->time[ip];

	  /* Calculate advection... */
	  module_advection(&ctl, met_t0, met_t1, atm, ip, dt);

	  /* Calculate diffusion... */
	  module_diffusion(&ctl, met_t0, met_t1, atm, ip, dt,
			   rng[omp_get_thread_num()]);

	  /* Calculate sedimentation... */
	  module_sedi(&ctl, atm, ip, dt);

	  /* Check position... */
	  module_position(met_t0, atm, ip);

	  /* Calculate decay of mass... */
	  module_decay(&ctl, atm, ip, dt);
	}
      STOP_TIMER(timer_phys);

      /* Write output... */
      START_TIMER(timer_output);
      write_output(dirname, &ctl, atm, t, t == ctl.t_stop);
      STOP_TIMER(timer_output);
    }

    /* ------------------------------------------------------------
       Finalize model run...
       ------------------------------------------------------------ */

    /* Free random number generators... */
    for (i = 0; i < NTHREADS; i++)
      gsl_rng_free(rng[i]);

    /* Free... */
    free(atm);
    free(met_t0);
    free(met_t1);

    /* Stop timer for total runtime... */
    STOP_TIMER(timer_total);

    /* Report timers... */
    PRINT_TIMER(timer_total);
    PRINT_TIMER(timer_phys);
    PRINT_TIMER(timer_input);
    PRINT_TIMER(timer_output);
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
  ctl_t * ctl,
  met_t * met_t0,
  met_t * met_t1,
  atm_t * atm,
  int ip,
  double dt) {

  double t0, t1, v[3], x[3], xout[3];

  /* Copy air parcel data... */
  x[0] = atm->lon[ip];
  x[1] = atm->lat[ip];
  x[2] = atm->p[ip];

  /* Interpolate meteorological data... */
  intpol_met_time(met_t0, met_t1, atm->time[ip], x[2], x[0], x[1], &t0,
		  &v[0], &v[1], &v[2]);

  /* Get position of the mid point... */
  xout[0] = x[0] + dx2deg(0.5 * dt * v[0] / 1000., x[1]);
  xout[1] = x[1] + dy2deg(0.5 * dt * v[1] / 1000.);
  xout[2] = x[2] + 0.5 * dt * v[2];

  /* Interpolate meteorological data for mid point... */
  intpol_met_time(met_t0, met_t1, atm->time[ip] + 0.5 * dt,
		  xout[2], xout[0], xout[1], &t1, &v[0], &v[1], &v[2]);

  /* Get new position... */
  xout[0] = x[0] + dx2deg(dt * v[0] / 1000., x[1]);
  xout[1] = x[1] + dy2deg(dt * v[1] / 1000.);
  xout[2] = x[2] + dt * v[2];

  /* Save new position... */
  atm->time[ip] += dt;
  atm->lon[ip] = xout[0];
  atm->lat[ip] = xout[1];
  atm->p[ip] = xout[2];

  /* Extrapolate temperature... */
  if (ctl->qnt_temp >= 0)
    atm->q[ctl->qnt_temp][ip] = 2. * t1 - t0;
}

/*****************************************************************************/

void module_decay(
  ctl_t * ctl,
  atm_t * atm,
  int ip,
  double dt) {

  /* Calculate exponential decay... */
  if (ctl->t12 > 0 && ctl->qnt_mass >= 0)
    atm->q[ctl->qnt_mass][ip] *= exp(-dt / (ctl->t12 * 86400. / log(2.)));
}

/*****************************************************************************/

void module_diffusion(
  ctl_t * ctl,
  met_t * met_t0,
  met_t * met_t1,
  atm_t * atm,
  int ip,
  double dt,
  gsl_rng * rng) {

  double r, rs, u[16], v[16], w[16], usig, vsig, wsig;

  int ix, iy, iz;

  /* Calculate mesoscale velocity fluctuations... */
  if (ctl->turb_meso > 0) {

    /* Get indices... */
    ix = locate(met_t0->lon, met_t0->nx, atm->lon[ip]);
    iy = locate(met_t0->lat, met_t0->ny, atm->lat[ip]);
    iz = locate(met_t0->p, met_t0->np, atm->p[ip]);

    /* Collect local wind data... */
    u[0] = met_t0->u[ix][iy][iz];
    u[1] = met_t0->u[ix + 1][iy][iz];
    u[2] = met_t0->u[ix][iy + 1][iz];
    u[3] = met_t0->u[ix + 1][iy + 1][iz];
    u[4] = met_t0->u[ix][iy][iz + 1];
    u[5] = met_t0->u[ix + 1][iy][iz + 1];
    u[6] = met_t0->u[ix][iy + 1][iz + 1];
    u[7] = met_t0->u[ix + 1][iy + 1][iz + 1];

    v[0] = met_t0->v[ix][iy][iz];
    v[1] = met_t0->v[ix + 1][iy][iz];
    v[2] = met_t0->v[ix][iy + 1][iz];
    v[3] = met_t0->v[ix + 1][iy + 1][iz];
    v[4] = met_t0->v[ix][iy][iz + 1];
    v[5] = met_t0->v[ix + 1][iy][iz + 1];
    v[6] = met_t0->v[ix][iy + 1][iz + 1];
    v[7] = met_t0->v[ix + 1][iy + 1][iz + 1];

    w[0] = met_t0->w[ix][iy][iz];
    w[1] = met_t0->w[ix + 1][iy][iz];
    w[2] = met_t0->w[ix][iy + 1][iz];
    w[3] = met_t0->w[ix + 1][iy + 1][iz];
    w[4] = met_t0->w[ix][iy][iz + 1];
    w[5] = met_t0->w[ix + 1][iy][iz + 1];
    w[6] = met_t0->w[ix][iy + 1][iz + 1];
    w[7] = met_t0->w[ix + 1][iy + 1][iz + 1];

    /* Get indices... */
    ix = locate(met_t1->lon, met_t1->nx, atm->lon[ip]);
    iy = locate(met_t1->lat, met_t1->ny, atm->lat[ip]);
    iz = locate(met_t1->p, met_t1->np, atm->p[ip]);

    /* Collect local wind data... */
    u[8] = met_t1->u[ix][iy][iz];
    u[9] = met_t1->u[ix + 1][iy][iz];
    u[10] = met_t1->u[ix][iy + 1][iz];
    u[11] = met_t1->u[ix + 1][iy + 1][iz];
    u[12] = met_t1->u[ix][iy][iz + 1];
    u[13] = met_t1->u[ix + 1][iy][iz + 1];
    u[14] = met_t1->u[ix][iy + 1][iz + 1];
    u[15] = met_t1->u[ix + 1][iy + 1][iz + 1];

    v[8] = met_t1->v[ix][iy][iz];
    v[9] = met_t1->v[ix + 1][iy][iz];
    v[10] = met_t1->v[ix][iy + 1][iz];
    v[11] = met_t1->v[ix + 1][iy + 1][iz];
    v[12] = met_t1->v[ix][iy][iz + 1];
    v[13] = met_t1->v[ix + 1][iy][iz + 1];
    v[14] = met_t1->v[ix][iy + 1][iz + 1];
    v[15] = met_t1->v[ix + 1][iy + 1][iz + 1];

    w[8] = met_t1->w[ix][iy][iz];
    w[9] = met_t1->w[ix + 1][iy][iz];
    w[10] = met_t1->w[ix][iy + 1][iz];
    w[11] = met_t1->w[ix + 1][iy + 1][iz];
    w[12] = met_t1->w[ix][iy][iz + 1];
    w[13] = met_t1->w[ix + 1][iy][iz + 1];
    w[14] = met_t1->w[ix][iy + 1][iz + 1];
    w[15] = met_t1->w[ix + 1][iy + 1][iz + 1];

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

  /* Horizontal turbulent diffusion... */
  if (ctl->turb_dx > 0) {
    atm->lon[ip]
      +=
      dx2deg(gsl_ran_gaussian_ziggurat
	     (rng, sqrt(2.0 * ctl->turb_dx * fabs(dt)))
	     / 1000., atm->lat[ip]);
    atm->lat[ip]
      +=
      dy2deg(gsl_ran_gaussian_ziggurat
	     (rng, sqrt(2.0 * ctl->turb_dx * fabs(dt)))
	     / 1000.);
  }

  /* Vertical turbulent diffusion... */
  if (ctl->turb_dz > 0)
    atm->p[ip]
      +=
      dz2dp(gsl_ran_gaussian_ziggurat
	    (rng, sqrt(2.0 * ctl->turb_dz * fabs(dt)))
	    / 1000., atm->p[ip]);
}

/*****************************************************************************/

void module_position(
  met_t * met,
  atm_t * atm,
  int ip) {

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

  /* Check pressure... */
  if (atm->p[ip] > met->p[0])
    atm->p[ip] = met->p[0];
  else if (atm->p[ip] < met->p[met->np - 1])
    atm->p[ip] = met->p[met->np - 1];
}

/*****************************************************************************/

void module_sedi(
  ctl_t * ctl,
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
  if (ctl->qnt_r_p < 0 || ctl->qnt_rho_p < 0 || ctl->qnt_temp < 0)
    return;

  /* Convert units... */
  p = 100 * atm->p[ip];
  T = atm->q[ctl->qnt_temp][ip];
  r_p = 1e-6 * atm->q[ctl->qnt_r_p][ip];
  rho_p = atm->q[ctl->qnt_rho_p][ip];

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
  double t,
  int force) {

  char filename[LEN];

  double r;

  int year, mon, day, hour, min, sec;

  /* Get time... */
  jsec2time(t, &year, &mon, &day, &hour, &min, &sec, &r);

  /* Write atmospheric data... */
  if (ctl->atm_dt_out > 0 && (fmod(t, ctl->atm_dt_out) == 0 || force)) {
    sprintf(filename, "%s_%04d_%02d_%02d_%02d_%02d.%s",
	    ctl->atm_basename, year, mon, day, hour, min,
	    (ctl->atm_oformat == 0 ? "tab" : "nc"));
    write_atm(dirname, filename, atm, ctl);
  }

  /* Write CSI data... */
  if (ctl->csi_dt_update > 0 && (fmod(t, ctl->csi_dt_update) == 0 || force)) {
    sprintf(filename, "%s_%04d_%02d_%02d_%02d_%02d.tab",
	    ctl->csi_basename, year, mon, day, hour, min);
    write_csi(dirname, filename, atm, ctl, t, ctl->csi_dt_update,
	      fmod(t, ctl->csi_dt_out) == 0 || force);
  }

  /* Write gridded data... */
  if (ctl->grid_dt_out > 0 && (fmod(t, ctl->grid_dt_out) == 0 || force)) {
    sprintf(filename, "%s_%04d_%02d_%02d_%02d_%02d.tab",
	    ctl->grid_basename, year, mon, day, hour, min);
    write_grid(dirname, filename, atm, ctl, t, ctl->grid_dt_out);
  }

  /* Write station data... */
  if (ctl->stat_dt_out > 0 && (fmod(t, ctl->stat_dt_out) == 0 || force)) {
    sprintf(filename, "%s_%04d_%02d_%02d_%02d_%02d.tab",
	    ctl->stat_basename, year, mon, day, hour, min);
    write_station(dirname, filename, atm, ctl, t, ctl->stat_dt_out);
  }
}
