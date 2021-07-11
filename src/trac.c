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
  Lagrangian particle dispersion model.
*/

#include "libtrac.h"

/* ------------------------------------------------------------
   Global variables...
   ------------------------------------------------------------ */

#ifdef _OPENACC
curandGenerator_t rng;
#else
static gsl_rng *rng[NTHREADS];
#endif

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Calculate advection of air parcels. */
void module_advection(
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/*! Apply boundary conditions. */
void module_bound_cond(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/*! Calculate convection of air parcels. */
void module_convection(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt,
  double *rs);

/*! Calculate exponential decay of particle mass. */
void module_decay(
  ctl_t * ctl,
  atm_t * atm,
  double *dt);

/*! Calculate mesoscale diffusion. */
void module_diffusion_meso(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  cache_t * cache,
  double *dt,
  double *rs);

/*! Calculate turbulent diffusion. */
void module_diffusion_turb(
  ctl_t * ctl,
  atm_t * atm,
  double *dt,
  double *rs);

/*! Calculate dry deposition. */
void module_dry_deposition(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/*! Initialize isosurface module. */
void module_isosurf_init(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  cache_t * cache);

/*! Force air parcels to stay on isosurface. */
void module_isosurf(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  cache_t * cache);

/*! Interpolate meteorological data for air parcel positions. */
void module_meteo(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm);

/*! Check position of air parcels. */
void module_position(
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/*! Initialize random number generator... */
void module_rng_init(
  int ntask);

/*! Generate random numbers. */
void module_rng(
  double *rs,
  size_t n,
  int method);

/*! Calculate sedimentation of air parcels. */
void module_sedi(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/*! Calculate OH chemistry. */
void module_oh_chem(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/*! Calculate wet deposition. */
void module_wet_deposition(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/*! Write simulation output. */
void write_output(
  const char *dirname,
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
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

  cache_t *cache;

  met_t *met0, *met1;

  FILE *dirlist;

  char dirname[LEN], filename[2 * LEN];

  double *dt, *rs, t;

  int num_devices = 0, ntask = -1, rank = 0, size = 1;

  /* Start timers... */
  START_TIMERS;

  /* Initialize MPI... */
#ifdef MPI
  SELECT_TIMER("MPI_INIT", NVTX_CPU);
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  /* Initialize GPUs... */
#ifdef _OPENACC
  SELECT_TIMER("ACC_INIT", NVTX_GPU);
  num_devices = acc_get_num_devices(acc_device_nvidia);
  if (num_devices <= 0)
    ERRMSG("Not running on a GPU device!");
  int device_num = rank % num_devices;
  acc_set_device_num(device_num, acc_device_nvidia);
  acc_device_t device_type = acc_get_device_type();
  acc_init(device_type);
#endif

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <dirlist> <ctl> <atm_in>");

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

    /* Allocate... */
    SELECT_TIMER("ALLOC", NVTX_CPU);
    ALLOC(atm, atm_t, 1);
    ALLOC(cache, cache_t, 1);
    ALLOC(met0, met_t, 1);
    ALLOC(met1, met_t, 1);
    ALLOC(dt, double,
	  NP);
    ALLOC(rs, double,
	  3 * NP + 1);

    /* Create data region on GPUs... */
#ifdef _OPENACC
    SELECT_TIMER("CREATE_DATA_REGION", NVTX_GPU);
#pragma acc enter data create(atm[:1],cache[:1],ctl,met0[:1],met1[:1],dt[:NP],rs[:3*NP])
#endif

    /* Read control parameters... */
    sprintf(filename, "%s/%s", dirname, argv[2]);
    read_ctl(filename, argc, argv, &ctl);

    /* Read atmospheric data... */
    sprintf(filename, "%s/%s", dirname, argv[3]);
    if (!read_atm(filename, &ctl, atm))
      ERRMSG("Cannot open file!");

    /* Set start time... */
    SELECT_TIMER("TIMESTEPS", NVTX_CPU);
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
      ERRMSG("Nothing to do! Check T_STOP and DIRECTION!");

    /* Round start time... */
    if (ctl.direction == 1)
      ctl.t_start = floor(ctl.t_start / ctl.dt_mod) * ctl.dt_mod;
    else
      ctl.t_start = ceil(ctl.t_start / ctl.dt_mod) * ctl.dt_mod;

    /* Update GPU... */
#ifdef _OPENACC
    SELECT_TIMER("UPDATE_DEVICE", NVTX_H2D);
#pragma acc update device(atm[:1],ctl)
#endif

    /* Initialize random number generator... */
    module_rng_init(ntask);

    /* Initialize meteorological data... */
    get_met(&ctl, ctl.t_start, &met0, &met1);
    if (ctl.dt_mod > fabs(met0->lon[1] - met0->lon[0]) * 111132. / 150.)
      WARN("Violation of CFL criterion! Check DT_MOD!");

    /* Initialize isosurface... */
    if (ctl.isosurf >= 1 && ctl.isosurf <= 4)
      module_isosurf_init(&ctl, met0, met1, atm, cache);

    /* Update GPU... */
#ifdef _OPENACC
    SELECT_TIMER("UPDATE_DEVICE", NVTX_H2D);
#pragma acc update device(cache[:1])
#endif

    /* ------------------------------------------------------------
       Loop over timesteps...
       ------------------------------------------------------------ */

    /* Loop over timesteps... */
    for (t = ctl.t_start; ctl.direction * (t - ctl.t_stop) < ctl.dt_mod;
	 t += ctl.direction * ctl.dt_mod) {

      /* Adjust length of final time step... */
      if (ctl.direction * (t - ctl.t_stop) > 0)
	t = ctl.t_stop;

      /* Set time steps for air parcels... */
      SELECT_TIMER("TIMESTEPS", NVTX_GPU);
#ifdef _OPENACC
#pragma acc parallel loop independent gang vector present(ctl,atm,atm->time,dt)
#endif
      for (int ip = 0; ip < atm->np; ip++) {
	double atmtime = atm->time[ip];
	double tstart = ctl.t_start;
	double tstop = ctl.t_stop;
	int dir = ctl.direction;
	if ((dir * (atmtime - tstart) >= 0 && dir * (atmtime - tstop) <= 0
	     && dir * (atmtime - t) < 0))
	  dt[ip] = t - atmtime;
	else
	  dt[ip] = 0;
      }

      /* Get meteorological data... */
      if (t != ctl.t_start)
	get_met(&ctl, t, &met0, &met1);

      /* Check initial positions... */
      module_position(met0, met1, atm, dt);

      /* Advection... */
      module_advection(met0, met1, atm, dt);

      /* Turbulent diffusion... */
      if (ctl.turb_dx_trop > 0 || ctl.turb_dz_trop > 0
	  || ctl.turb_dx_strat > 0 || ctl.turb_dz_strat > 0)
	module_diffusion_turb(&ctl, atm, dt, rs);

      /* Mesoscale diffusion... */
      if (ctl.turb_mesox > 0 || ctl.turb_mesoz > 0)
	module_diffusion_meso(&ctl, met0, met1, atm, cache, dt, rs);

      /* Convection... */
      if (ctl.conv_cape >= 0
	  && (ctl.conv_dt <= 0 || fmod(t, ctl.conv_dt) == 0))
	module_convection(&ctl, met0, met1, atm, dt, rs);

      /* Sedimentation... */
      if (ctl.qnt_r >= 0 && ctl.qnt_rho >= 0)
	module_sedi(&ctl, met0, met1, atm, dt);

      /* Isosurface... */
      if (ctl.isosurf >= 1 && ctl.isosurf <= 4)
	module_isosurf(&ctl, met0, met1, atm, cache);

      /* Check final positions... */
      module_position(met0, met1, atm, dt);

      /* Interpolate meteorological data... */
      if (ctl.met_dt_out > 0
	  && (ctl.met_dt_out < ctl.dt_mod || fmod(t, ctl.met_dt_out) == 0))
	module_meteo(&ctl, met0, met1, atm);

      /* Decay of particle mass... */
      if (ctl.tdec_trop > 0 && ctl.tdec_strat > 0)
	module_decay(&ctl, atm, dt);

      /* OH chemistry... */
      if (ctl.oh_chem[0] > 0 && ctl.oh_chem[2] > 0)
	module_oh_chem(&ctl, met0, met1, atm, dt);

      /* Dry deposition... */
      if (ctl.dry_depo[0] > 0)
	module_dry_deposition(&ctl, met0, met1, atm, dt);

      /* Wet deposition... */
      if ((ctl.wet_depo[0] > 0 || ctl.wet_depo[2] > 0)
	  && (ctl.wet_depo[4] > 0 || ctl.wet_depo[6] > 0))
	module_wet_deposition(&ctl, met0, met1, atm, dt);

      /* Boundary conditions... */
      if (ctl.bound_mass > 0 || ctl.bound_vmr > 0)
	module_bound_cond(&ctl, met0, met1, atm, dt);

      /* Write output... */
      write_output(dirname, &ctl, met0, met1, atm, t);
    }

    /* ------------------------------------------------------------
       Finalize model run...
       ------------------------------------------------------------ */

    /* Report problem size... */
    LOG(1, "SIZE_NP = %d", atm->np);
    LOG(1, "SIZE_MPI_TASKS = %d", size);
    LOG(1, "SIZE_OMP_THREADS = %d", omp_get_max_threads());
    LOG(1, "SIZE_ACC_DEVICES = %d", num_devices);

    /* Report memory usage... */
    LOG(1, "MEMORY_ATM = %g MByte", sizeof(atm_t) / 1024. / 1024.);
    LOG(1, "MEMORY_CACHE = %g MByte", sizeof(cache_t) / 1024. / 1024.);
    LOG(1, "MEMORY_METEO = %g MByte", 2 * sizeof(met_t) / 1024. / 1024.);
    LOG(1, "MEMORY_DYNAMIC = %g MByte", (sizeof(met_t)
					 + 4 * NP * sizeof(double)
					 + EX * EY * EP * sizeof(float)) /
	1024. / 1024.);
    LOG(1, "MEMORY_STATIC = %g MByte", (EX * EY * sizeof(double)
					+ EX * EY * EP * sizeof(float)
					+ 4 * GX * GY * GZ * sizeof(double)
					+ 2 * GX * GY * GZ * sizeof(int)
					+ 2 * GX * GY * sizeof(double)
					+ GX * GY * sizeof(int)) / 1024. /
	1024.);

    /* Delete data region on GPUs... */
#ifdef _OPENACC
    SELECT_TIMER("DELETE_DATA_REGION", NVTX_GPU);
#pragma acc exit data delete(ctl,atm,cache,met0,met1,dt,rs)
#endif

    /* Free... */
    SELECT_TIMER("FREE", NVTX_CPU);
    free(atm);
    free(cache);
    free(met0);
    free(met1);
    free(dt);
    free(rs);

    /* Report timers... */
    PRINT_TIMERS;
  }

  /* Finalize MPI... */
#ifdef MPI
  MPI_Finalize();
#endif

  /* Stop timers... */
  STOP_TIMERS;

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void module_advection(
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt) {

  /* Set timer... */
  SELECT_TIMER("MODULE_ADVECTION", NVTX_GPU);

#ifdef _OPENACC
#pragma acc data present(met0,met1,atm,dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++)
    if (dt[ip] != 0) {

      double u, v, w;

      /* Interpolate meteorological data... */
      INTPOL_INIT;
      INTPOL_3D(u, 1);
      INTPOL_3D(v, 0);
      INTPOL_3D(w, 0);

      /* Get position of the mid point... */
      double dtm = atm->time[ip] + 0.5 * dt[ip];
      double xm0 =
	atm->lon[ip] + DX2DEG(0.5 * dt[ip] * u / 1000., atm->lat[ip]);
      double xm1 = atm->lat[ip] + DY2DEG(0.5 * dt[ip] * v / 1000.);
      double xm2 = atm->p[ip] + 0.5 * dt[ip] * w;

      /* Interpolate meteorological data for mid point... */
      intpol_met_time_3d(met0, met0->u, met1, met1->u,
			 dtm, xm2, xm0, xm1, &u, ci, cw, 1);
      intpol_met_time_3d(met0, met0->v, met1, met1->v,
			 dtm, xm2, xm0, xm1, &v, ci, cw, 0);
      intpol_met_time_3d(met0, met0->w, met1, met1->w,
			 dtm, xm2, xm0, xm1, &w, ci, cw, 0);

      /* Save new position... */
      atm->time[ip] += dt[ip];
      atm->lon[ip] += DX2DEG(dt[ip] * u / 1000., xm1);
      atm->lat[ip] += DY2DEG(dt[ip] * v / 1000.);
      atm->p[ip] += dt[ip] * w;
    }
}

/*****************************************************************************/

void module_bound_cond(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt) {

  /* Set timer... */
  SELECT_TIMER("MODULE_BOUNDCOND", NVTX_GPU);

  /* Check quantity flags... */
  if (ctl->qnt_m < 0 && ctl->qnt_vmr < 0)
    ERRMSG("Module needs quantity mass or volume mixing ratio!");

#ifdef _OPENACC
#pragma acc data present(ctl,met0,met1,atm,dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++)
    if (dt[ip] != 0) {

      double ps;

      /* Check latitude and pressure range... */
      if (atm->lat[ip] < ctl->bound_lat0 || atm->lat[ip] > ctl->bound_lat1
	  || atm->p[ip] > ctl->bound_p0 || atm->p[ip] < ctl->bound_p1)
	continue;

      /* Check surface layer... */
      if (ctl->bound_dps > 0) {

	/* Get surface pressure... */
	INTPOL_INIT;
	INTPOL_2D(ps, 1);

	/* Check whether particle is above the surface layer... */
	if (atm->p[ip] < ps - ctl->bound_dps)
	  continue;
      }

      /* Set mass and volume mixing ratio... */
      if (ctl->qnt_m >= 0)
	atm->q[ctl->qnt_m][ip] = ctl->bound_mass;
      if (ctl->qnt_vmr >= 0)
	atm->q[ctl->qnt_vmr][ip] = ctl->bound_vmr;
    }
}

/*****************************************************************************/

void module_convection(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt,
  double *rs) {

  /* Set timer... */
  SELECT_TIMER("MODULE_CONVECTION", NVTX_GPU);

  /* Create random numbers... */
  module_rng(rs, (size_t) atm->np, 0);

#ifdef _OPENACC
#pragma acc data present(ctl,met0,met1,atm,dt,rs)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++)
    if (dt[ip] != 0) {

      double cape, pbot, pel, ps, ptop;

      /* Interpolate CAPE... */
      INTPOL_INIT;
      INTPOL_2D(cape, 1);

      /* Check threshold... */
      if (isfinite(cape) && cape >= ctl->conv_cape) {

	/* Interpolate equilibrium level... */
	INTPOL_2D(pel, 0);

	/* Check whether particle is above cloud top... */
	if (!isfinite(pel) || atm->p[ip] < pel)
	  continue;

	/* Set pressure range for mixing... */
	pbot = ptop = atm->p[ip];
	if (ctl->conv_mix_bot == 1) {
	  INTPOL_2D(ps, 0);
	  pbot = ps;
	}
	if (ctl->conv_mix_top == 1)
	  ptop = pel;

	/* Mix in entire in cloud column... */
	atm->p[ip] = pbot + (ptop - pbot) * rs[ip];
      }
    }
}

/*****************************************************************************/

void module_decay(
  ctl_t * ctl,
  atm_t * atm,
  double *dt) {

  /* Set timer... */
  SELECT_TIMER("MODULE_DECAY", NVTX_GPU);

  /* Check quantity flags... */
  if (ctl->qnt_m < 0 && ctl->qnt_vmr < 0)
    ERRMSG("Module needs quantity mass or volume mixing ratio!");

#ifdef _OPENACC
#pragma acc data present(ctl,atm,dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++)
    if (dt[ip] != 0) {

      /* Get tropopause pressure... */
      double pt = clim_tropo(atm->time[ip], atm->lat[ip]);

      /* Get weighting factor... */
      double w;
      double p1 = pt * 0.866877899;
      double p0 = pt / 0.866877899;
      if (atm->p[ip] > p0)
	w = 1;
      else if (atm->p[ip] < p1)
	w = 0;
      else
	w = LIN(p0, 1.0, p1, 0.0, atm->p[ip]);

      /* Set lifetime... */
      double tdec = w * ctl->tdec_trop + (1 - w) * ctl->tdec_strat;

      /* Calculate exponential decay... */
      double aux = exp(-dt[ip] / tdec);
      if (ctl->qnt_m >= 0)
	atm->q[ctl->qnt_m][ip] *= aux;
      if (ctl->qnt_vmr >= 0)
	atm->q[ctl->qnt_vmr][ip] *= aux;
    }
}

/*****************************************************************************/

void module_diffusion_meso(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  cache_t * cache,
  double *dt,
  double *rs) {

  /* Set timer... */
  SELECT_TIMER("MODULE_TURBMESO", NVTX_GPU);

  /* Create random numbers... */
  module_rng(rs, 3 * (size_t) atm->np, 1);

#ifdef _OPENACC
#pragma acc data present(ctl,met0,met1,atm,cache,dt,rs)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++)
    if (dt[ip] != 0) {

      double u[16], v[16], w[16];

      /* Get indices... */
      int ix = locate_reg(met0->lon, met0->nx, atm->lon[ip]);
      int iy = locate_reg(met0->lat, met0->ny, atm->lat[ip]);
      int iz = locate_irr(met0->p, met0->np, atm->p[ip]);

      /* Caching of wind standard deviations... */
      if (cache->tsig[ix][iy][iz] != met0->time) {

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
	cache->usig[ix][iy][iz] = (float) stddev(u, 16);
	cache->vsig[ix][iy][iz] = (float) stddev(v, 16);
	cache->wsig[ix][iy][iz] = (float) stddev(w, 16);
	cache->tsig[ix][iy][iz] = met0->time;
      }

      /* Set temporal correlations for mesoscale fluctuations... */
      double r = 1 - 2 * fabs(dt[ip]) / ctl->dt_met;
      double r2 = sqrt(1 - r * r);

      /* Calculate horizontal mesoscale wind fluctuations... */
      if (ctl->turb_mesox > 0) {
	cache->up[ip] = (float)
	  (r * cache->up[ip]
	   + r2 * rs[3 * ip] * ctl->turb_mesox * cache->usig[ix][iy][iz]);
	atm->lon[ip] += DX2DEG(cache->up[ip] * dt[ip] / 1000., atm->lat[ip]);

	cache->vp[ip] = (float)
	  (r * cache->vp[ip]
	   + r2 * rs[3 * ip + 1] * ctl->turb_mesox * cache->vsig[ix][iy][iz]);
	atm->lat[ip] += DY2DEG(cache->vp[ip] * dt[ip] / 1000.);
      }

      /* Calculate vertical mesoscale wind fluctuations... */
      if (ctl->turb_mesoz > 0) {
	cache->wp[ip] = (float)
	  (r * cache->wp[ip]
	   + r2 * rs[3 * ip + 2] * ctl->turb_mesoz * cache->wsig[ix][iy][iz]);
	atm->p[ip] += cache->wp[ip] * dt[ip];
      }
    }
}

/*****************************************************************************/

void module_diffusion_turb(
  ctl_t * ctl,
  atm_t * atm,
  double *dt,
  double *rs) {

  /* Set timer... */
  SELECT_TIMER("MODULE_TURBDIFF", NVTX_GPU);

  /* Create random numbers... */
  module_rng(rs, 3 * (size_t) atm->np, 1);

#ifdef _OPENACC
#pragma acc data present(ctl,atm,dt,rs)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++)
    if (dt[ip] != 0) {

      /* Get tropopause pressure... */
      double pt = clim_tropo(atm->time[ip], atm->lat[ip]);

      /* Get weighting factor... */
      double w;
      double p1 = pt * 0.866877899;
      double p0 = pt / 0.866877899;
      if (atm->p[ip] > p0)
	w = 1;
      else if (atm->p[ip] < p1)
	w = 0;
      else
	w = LIN(p0, 1.0, p1, 0.0, atm->p[ip]);

      /* Set diffusivity... */
      double dx = w * ctl->turb_dx_trop + (1 - w) * ctl->turb_dx_strat;
      double dz = w * ctl->turb_dz_trop + (1 - w) * ctl->turb_dz_strat;

      /* Horizontal turbulent diffusion... */
      if (dx > 0) {
	double sigma = sqrt(2.0 * dx * fabs(dt[ip]));
	atm->lon[ip] += DX2DEG(rs[3 * ip] * sigma / 1000., atm->lat[ip]);
	atm->lat[ip] += DY2DEG(rs[3 * ip + 1] * sigma / 1000.);
      }

      /* Vertical turbulent diffusion... */
      if (dz > 0) {
	double sigma = sqrt(2.0 * dz * fabs(dt[ip]));
	atm->p[ip]
	  += DZ2DP(rs[3 * ip + 2] * sigma / 1000., atm->p[ip]);
      }
    }
}

/*****************************************************************************/

void module_dry_deposition(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt) {

  /* Set timer... */
  SELECT_TIMER("MODULE_DRYDEPO", NVTX_GPU);

  /* Width of the surface layer [hPa]. */
  const double dp = 30.;

  /* Check quantity flags... */
  if (ctl->qnt_m < 0 && ctl->qnt_vmr < 0)
    ERRMSG("Module needs quantity mass or volume mixing ratio!");

#ifdef _OPENACC
#pragma acc data present(ctl,met0,met1,atm,dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++)
    if (dt[ip] != 0) {

      double ps, t, v_dep;

      /* Get surface pressure... */
      INTPOL_INIT;
      INTPOL_2D(ps, 1);

      /* Check whether particle is above the surface layer... */
      if (atm->p[ip] < ps - dp)
	continue;

      /* Set width of surface layer... */
      double dz = 1000. * (Z(ps - dp) - Z(ps));

      /* Calculate sedimentation velocity for particles... */
      if (ctl->qnt_r > 0 && ctl->qnt_rho > 0) {

	/* Get temperature... */
	INTPOL_3D(t, 1);

	/* Set deposition velocity... */
	v_dep = sedi(atm->p[ip], t, atm->q[ctl->qnt_r][ip],
		     atm->q[ctl->qnt_rho][ip]);
      }

      /* Use explicit sedimentation velocity for gases... */
      else
	v_dep = ctl->dry_depo[0];

      /* Calculate loss of mass based on deposition velocity... */
      double aux = exp(-dt[ip] * v_dep / dz);
      if (ctl->qnt_m >= 0)
	atm->q[ctl->qnt_m][ip] *= aux;
      if (ctl->qnt_vmr >= 0)
	atm->q[ctl->qnt_vmr][ip] *= aux;
    }
}

/*****************************************************************************/

void module_isosurf_init(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  cache_t * cache) {

  FILE *in;

  char line[LEN];

  double t;

  /* Set timer... */
  SELECT_TIMER("MODULE_ISOSURF", NVTX_GPU);

  /* Init... */
  INTPOL_INIT;

  /* Save pressure... */
  if (ctl->isosurf == 1)
    for (int ip = 0; ip < atm->np; ip++)
      cache->iso_var[ip] = atm->p[ip];

  /* Save density... */
  else if (ctl->isosurf == 2)
    for (int ip = 0; ip < atm->np; ip++) {
      INTPOL_3D(t, 1);
      cache->iso_var[ip] = atm->p[ip] / t;
    }

  /* Save potential temperature... */
  else if (ctl->isosurf == 3)
    for (int ip = 0; ip < atm->np; ip++) {
      INTPOL_3D(t, 1);
      cache->iso_var[ip] = THETA(atm->p[ip], t);
    }

  /* Read balloon pressure data... */
  else if (ctl->isosurf == 4) {

    /* Write info... */
    LOG(1, "Read balloon pressure data: %s", ctl->balloon);

    /* Open file... */
    if (!(in = fopen(ctl->balloon, "r")))
      ERRMSG("Cannot open file!");

    /* Read pressure time series... */
    while (fgets(line, LEN, in))
      if (sscanf(line, "%lg %lg", &(cache->iso_ts[cache->iso_n]),
		 &(cache->iso_ps[cache->iso_n])) == 2)
	if ((++cache->iso_n) > NP)
	  ERRMSG("Too many data points!");

    /* Check number of points... */
    if (cache->iso_n < 1)
      ERRMSG("Could not read any data!");

    /* Close file... */
    fclose(in);
  }
}

/*****************************************************************************/

void module_isosurf(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  cache_t * cache) {

  /* Set timer... */
  SELECT_TIMER("MODULE_ISOSURF", NVTX_GPU);

#ifdef _OPENACC
#pragma acc data present(ctl,met0,met1,atm,cache)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++) {

    double t;

    /* Init... */
    INTPOL_INIT;

    /* Restore pressure... */
    if (ctl->isosurf == 1)
      atm->p[ip] = cache->iso_var[ip];

    /* Restore density... */
    else if (ctl->isosurf == 2) {
      INTPOL_3D(t, 1);
      atm->p[ip] = cache->iso_var[ip] * t;
    }

    /* Restore potential temperature... */
    else if (ctl->isosurf == 3) {
      INTPOL_3D(t, 1);
      atm->p[ip] = 1000. * pow(cache->iso_var[ip] / t, -1. / 0.286);
    }

    /* Interpolate pressure... */
    else if (ctl->isosurf == 4) {
      if (atm->time[ip] <= cache->iso_ts[0])
	atm->p[ip] = cache->iso_ps[0];
      else if (atm->time[ip] >= cache->iso_ts[cache->iso_n - 1])
	atm->p[ip] = cache->iso_ps[cache->iso_n - 1];
      else {
	int idx = locate_irr(cache->iso_ts, cache->iso_n, atm->time[ip]);
	atm->p[ip] = LIN(cache->iso_ts[idx], cache->iso_ps[idx],
			 cache->iso_ts[idx + 1], cache->iso_ps[idx + 1],
			 atm->time[ip]);
      }
    }
  }
}

/*****************************************************************************/

void module_meteo(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm) {

  /* Set timer... */
  SELECT_TIMER("MODULE_METEO", NVTX_GPU);

  /* Check quantity flags... */
  if (ctl->qnt_tsts >= 0)
    if (ctl->qnt_tice < 0 || ctl->qnt_tnat < 0)
      ERRMSG("Need T_ice and T_NAT to calculate T_STS!");

#ifdef _OPENACC
#pragma acc data present(ctl,met0,met1,atm)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++) {

    double ps, ts, zs, us, vs, pbl, pt, pc, cl, plcl, plfc, pel, cape, pv, t,
      tt, u, v, w, h2o, h2ot, o3, lwc, iwc, z, zt;

    /* Interpolate meteorological data... */
    INTPOL_INIT;
    INTPOL_TIME_ALL(atm->time[ip], atm->p[ip], atm->lon[ip], atm->lat[ip]);

    /* Set quantities... */
    SET_ATM(qnt_ps, ps);
    SET_ATM(qnt_ts, ts);
    SET_ATM(qnt_zs, zs);
    SET_ATM(qnt_us, us);
    SET_ATM(qnt_vs, vs);
    SET_ATM(qnt_pbl, pbl);
    SET_ATM(qnt_pt, pt);
    SET_ATM(qnt_tt, tt);
    SET_ATM(qnt_zt, zt);
    SET_ATM(qnt_h2ot, h2ot);
    SET_ATM(qnt_p, atm->p[ip]);
    SET_ATM(qnt_z, z);
    SET_ATM(qnt_t, t);
    SET_ATM(qnt_u, u);
    SET_ATM(qnt_v, v);
    SET_ATM(qnt_w, w);
    SET_ATM(qnt_h2o, h2o);
    SET_ATM(qnt_o3, o3);
    SET_ATM(qnt_lwc, lwc);
    SET_ATM(qnt_iwc, iwc);
    SET_ATM(qnt_pc, pc);
    SET_ATM(qnt_cl, cl);
    SET_ATM(qnt_plcl, plcl);
    SET_ATM(qnt_plfc, plfc);
    SET_ATM(qnt_pel, pel);
    SET_ATM(qnt_cape, cape);
    SET_ATM(qnt_hno3, clim_hno3(atm->time[ip], atm->lat[ip], atm->p[ip]));
    SET_ATM(qnt_oh, clim_oh(atm->time[ip], atm->lat[ip], atm->p[ip]));
    SET_ATM(qnt_vh, sqrt(u * u + v * v));
    SET_ATM(qnt_vz, -1e3 * H0 / atm->p[ip] * w);
    SET_ATM(qnt_psat, PSAT(t));
    SET_ATM(qnt_psice, PSICE(t));
    SET_ATM(qnt_pw, PW(atm->p[ip], h2o));
    SET_ATM(qnt_sh, SH(h2o));
    SET_ATM(qnt_rh, RH(atm->p[ip], t, h2o));
    SET_ATM(qnt_rhice, RHICE(atm->p[ip], t, h2o));
    SET_ATM(qnt_theta, THETA(atm->p[ip], t));
    SET_ATM(qnt_zeta, ZETA(ps, atm->p[ip], t));
    SET_ATM(qnt_tvirt, TVIRT(t, h2o));
    SET_ATM(qnt_lapse, lapse_rate(t, h2o));
    SET_ATM(qnt_pv, pv);
    SET_ATM(qnt_tdew, TDEW(atm->p[ip], h2o));
    SET_ATM(qnt_tice, TICE(atm->p[ip], h2o));
    SET_ATM(qnt_tnat,
	    nat_temperature(atm->p[ip], h2o,
			    clim_hno3(atm->time[ip], atm->lat[ip],
				      atm->p[ip]) * 1e-9));
    SET_ATM(qnt_tsts,
	    0.5 * (atm->q[ctl->qnt_tice][ip] + atm->q[ctl->qnt_tnat][ip]));
  }
}

/*****************************************************************************/

void module_position(
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt) {

  /* Set timer... */
  SELECT_TIMER("MODULE_POSITION", NVTX_GPU);

#ifdef _OPENACC
#pragma acc data present(met0,met1,atm,dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++)
    if (dt[ip] != 0) {

      /* Init... */
      double ps;
      INTPOL_INIT;

      /* Calculate modulo... */
      atm->lon[ip] = FMOD(atm->lon[ip], 360.);
      atm->lat[ip] = FMOD(atm->lat[ip], 360.);

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
      if (atm->p[ip] < met0->p[met0->np - 1])
	atm->p[ip] = met0->p[met0->np - 1];
      else if (atm->p[ip] > 300.) {
	INTPOL_2D(ps, 1);
	if (atm->p[ip] > ps)
	  atm->p[ip] = ps;
      }
    }
}

/*****************************************************************************/

void module_rng_init(
  int ntask) {

  /* Initialize random number generator... */
#ifdef _OPENACC

  if (curandCreateGenerator(&rng, CURAND_RNG_PSEUDO_DEFAULT)
      != CURAND_STATUS_SUCCESS)
    ERRMSG("Cannot create random number generator!");
  if (curandSetPseudoRandomGeneratorSeed(rng, ntask)
      != CURAND_STATUS_SUCCESS)
    ERRMSG("Cannot set seed for random number generator!");
  if (curandSetStream(rng, (cudaStream_t) acc_get_cuda_stream(acc_async_sync))
      != CURAND_STATUS_SUCCESS)
    ERRMSG("Cannot set stream for random number generator!");

#else

  gsl_rng_env_setup();
  if (omp_get_max_threads() > NTHREADS)
    ERRMSG("Too many threads!");
  for (int i = 0; i < NTHREADS; i++) {
    rng[i] = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng[i], gsl_rng_default_seed
		+ (long unsigned) (ntask * NTHREADS + i));
  }

#endif
}

/*****************************************************************************/

void module_rng(
  double *rs,
  size_t n,
  int method) {

#ifdef _OPENACC

#pragma acc host_data use_device(rs)
  {
    /* Uniform distribution... */
    if (method == 0) {
      if (curandGenerateUniformDouble(rng, rs, (n < 4 ? 4 : n))
	  != CURAND_STATUS_SUCCESS)
	ERRMSG("Cannot create random numbers!");
    }

    /* Normal distribution... */
    else if (method == 1) {
      if (curandGenerateNormalDouble(rng, rs, (n < 4 ? 4 : n), 0.0, 1.0)
	  != CURAND_STATUS_SUCCESS)
	ERRMSG("Cannot create random numbers!");
    }
  }

#else

  /* Uniform distribution... */
  if (method == 0) {
#pragma omp parallel for default(shared)
    for (size_t i = 0; i < n; ++i)
      rs[i] = gsl_rng_uniform(rng[omp_get_thread_num()]);
  }

  /* Normal distribution... */
  else if (method == 1) {
#pragma omp parallel for default(shared)
    for (size_t i = 0; i < n; ++i)
      rs[i] = gsl_ran_gaussian_ziggurat(rng[omp_get_thread_num()], 1.0);
  }
#endif

}

/*****************************************************************************/

void module_sedi(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt) {

  /* Set timer... */
  SELECT_TIMER("MODULE_SEDI", NVTX_GPU);

#ifdef _OPENACC
#pragma acc data present(ctl,met0,met1,atm,dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++)
    if (dt[ip] != 0) {

      /* Get temperature... */
      double t;
      INTPOL_INIT;
      INTPOL_3D(t, 1);

      /* Sedimentation velocity... */
      double v_s = sedi(atm->p[ip], t, atm->q[ctl->qnt_r][ip],
			atm->q[ctl->qnt_rho][ip]);

      /* Calculate pressure change... */
      atm->p[ip] += DZ2DP(v_s * dt[ip] / 1000., atm->p[ip]);
    }
}

/*****************************************************************************/

void module_oh_chem(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt) {

  /* Set timer... */
  SELECT_TIMER("MODULE_OHCHEM", NVTX_GPU);

  /* Check quantity flags... */
  if (ctl->qnt_m < 0 && ctl->qnt_vmr < 0)
    ERRMSG("Module needs quantity mass or volume mixing ratio!");

#ifdef _OPENACC
#pragma acc data present(ctl,met0,met1,atm,dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++)
    if (dt[ip] != 0) {

      /* Get temperature... */
      double t;
      INTPOL_INIT;
      INTPOL_3D(t, 1);

      /* Calculate molecular density (IUPAC Data Sheet I.A4.86 SOx15)... */
      double M = 7.243e21 * (atm->p[ip] / 1000.) / t;

      /* Calculate rate coefficient for X + OH + M -> XOH + M
         (JPL Publication 15-10) ... */
      double k0 = ctl->oh_chem[0] *
	(ctl->oh_chem[1] > 0 ? pow(t / 300., -ctl->oh_chem[1]) : 1.);
      double ki = ctl->oh_chem[2] *
	(ctl->oh_chem[3] > 0 ? pow(t / 300., -ctl->oh_chem[3]) : 1.);
      double c = log10(k0 * M / ki);
      double k = k0 * M / (1. + k0 * M / ki) * pow(0.6, 1. / (1. + c * c));

      /* Calculate exponential decay... */
      double aux
	= exp(-dt[ip] * k * clim_oh(atm->time[ip], atm->lat[ip], atm->p[ip]));
      if (ctl->qnt_m >= 0)
	atm->q[ctl->qnt_m][ip] *= aux;
      if (ctl->qnt_vmr >= 0)
	atm->q[ctl->qnt_vmr][ip] *= aux;
    }
}

/*****************************************************************************/

void module_wet_deposition(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt) {

  /* Set timer... */
  SELECT_TIMER("MODULE_WETDEPO", NVTX_GPU);

  /* Check quantity flags... */
  if (ctl->qnt_m < 0 && ctl->qnt_vmr < 0)
    ERRMSG("Module needs quantity mass or volume mixing ratio!");

#ifdef _OPENACC
#pragma acc data present(ctl,met0,met1,atm,dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++)
    if (dt[ip] != 0) {

      double cl, dz, h, lambda = 0, t, iwc, lwc, pc;

      int inside;

      /* Check whether particle is below cloud top... */
      INTPOL_INIT;
      INTPOL_2D(pc, 1);
      if (!isfinite(pc) || atm->p[ip] <= pc)
	continue;

      /* Estimate precipitation rate (Pisso et al., 2019)... */
      INTPOL_2D(cl, 0);
      double Is = pow(2. * cl, 1. / 0.36);
      if (Is < 0.01)
	continue;

      /* Check whether particle is inside or below cloud... */
      INTPOL_3D(lwc, 1);
      INTPOL_3D(iwc, 0);
      inside = (iwc > 0 || lwc > 0);

      /* Calculate in-cloud scavenging coefficient... */
      if (inside) {

	/* Use exponential dependency for particles... */
	if (ctl->wet_depo[0] > 0)
	  lambda = ctl->wet_depo[0] * pow(Is, ctl->wet_depo[1]);

	/* Use Henry's law for gases... */
	else if (ctl->wet_depo[2] > 0) {

	  /* Get temperature... */
	  INTPOL_3D(t, 0);

	  /* Get Henry's constant (Sander, 2015)... */
	  h = ctl->wet_depo[2]
	    * exp(ctl->wet_depo[3] * (1. / t - 1. / 298.15));

	  /* Estimate depth of cloud layer... */
	  dz = 1e3 * Z(pc);

	  /* Calculate scavenging coefficient (Draxler and Hess, 1997)... */
	  lambda = h * RI * t * Is / 3.6e6 * dz;
	}
      }

      /* Calculate below-cloud scavenging coefficient... */
      else {

	/* Use exponential dependency for particles... */
	if (ctl->wet_depo[4] > 0)
	  lambda = ctl->wet_depo[4] * pow(Is, ctl->wet_depo[5]);

	/* Use Henry's law for gases... */
	else if (ctl->wet_depo[6] > 0) {

	  /* Get temperature... */
	  INTPOL_3D(t, 0);

	  /* Get Henry's constant (Sander, 2015)... */
	  h = ctl->wet_depo[6]
	    * exp(ctl->wet_depo[7] * (1. / t - 1. / 298.15));

	  /* Estimate depth of cloud layer... */
	  dz = 1e3 * Z(pc);

	  /* Calculate scavenging coefficient (Draxler and Hess, 1997)... */
	  lambda = h * RI * t * Is / 3.6e6 * dz;
	}
      }

      /* Calculate exponential decay of mass... */
      double aux = exp(-dt[ip] * lambda);
      if (ctl->qnt_m >= 0)
	atm->q[ctl->qnt_m][ip] *= aux;
      if (ctl->qnt_vmr >= 0)
	atm->q[ctl->qnt_vmr][ip] *= aux;
    }
}

/*****************************************************************************/

void write_output(
  const char *dirname,
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double t) {

  char filename[2 * LEN];

  double r;

  int year, mon, day, hour, min, sec;

  /* Get time... */
  jsec2time(t, &year, &mon, &day, &hour, &min, &sec, &r);

  /* Update host... */
#ifdef _OPENACC
  if ((ctl->atm_basename[0] != '-' && fmod(t, ctl->atm_dt_out) == 0)
      || (ctl->grid_basename[0] != '-' && fmod(t, ctl->grid_dt_out) == 0)
      || ctl->csi_basename[0] != '-' || ctl->ens_basename[0] != '-'
      || ctl->prof_basename[0] != '-' || ctl->sample_basename[0] != '-'
      || ctl->stat_basename[0] != '-') {
    SELECT_TIMER("UPDATE_HOST", NVTX_D2H);
#pragma acc update host(atm[:1])
  }
#endif

  /* Write atmospheric data... */
  if (ctl->atm_basename[0] != '-' && fmod(t, ctl->atm_dt_out) == 0) {
    sprintf(filename, "%s/%s_%04d_%02d_%02d_%02d_%02d.tab",
	    dirname, ctl->atm_basename, year, mon, day, hour, min);
    write_atm(filename, ctl, atm, t);
  }

  /* Write gridded data... */
  if (ctl->grid_basename[0] != '-' && fmod(t, ctl->grid_dt_out) == 0) {
    sprintf(filename, "%s/%s_%04d_%02d_%02d_%02d_%02d.tab",
	    dirname, ctl->grid_basename, year, mon, day, hour, min);
    write_grid(filename, ctl, met0, met1, atm, t);
  }

  /* Write CSI data... */
  if (ctl->csi_basename[0] != '-') {
    sprintf(filename, "%s/%s.tab", dirname, ctl->csi_basename);
    write_csi(filename, ctl, atm, t);
  }

  /* Write ensemble data... */
  if (ctl->ens_basename[0] != '-') {
    sprintf(filename, "%s/%s.tab", dirname, ctl->ens_basename);
    write_ens(filename, ctl, atm, t);
  }

  /* Write profile data... */
  if (ctl->prof_basename[0] != '-') {
    sprintf(filename, "%s/%s.tab", dirname, ctl->prof_basename);
    write_prof(filename, ctl, met0, met1, atm, t);
  }

  /* Write sample data... */
  if (ctl->sample_basename[0] != '-') {
    sprintf(filename, "%s/%s.tab", dirname, ctl->sample_basename);
    write_sample(filename, ctl, met0, met1, atm, t);
  }

  /* Write station data... */
  if (ctl->stat_basename[0] != '-') {
    sprintf(filename, "%s/%s.tab", dirname, ctl->stat_basename);
    write_station(filename, ctl, atm, t);
  }
}
