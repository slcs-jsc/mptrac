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

  Copyright (C) 2013-2023 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Lagrangian particle dispersion model.
*/

#include "libtrac.h"

#ifdef KPP
#include "chem_Parameters.h"
#include "chem_Global.h"
#include "chem_Sparse.h"
#include "kpp_chem.h"
#endif

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
void module_advect(
  ctl_t * ctl,
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
  clim_t * clim,
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
  clim_t * clim,
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

/*! Interpolate meteo data for air parcel positions. */
void module_meteo(
  ctl_t * ctl,
  clim_t * clim,
  met_t * met0,
  met_t * met1,
  atm_t * atm);

/*! Calculate OH chemistry. */
void module_oh_chem(
  ctl_t * ctl,
  clim_t * clim,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/*! Calculate H2O2 chemistry. */
void module_h2o2_chem(
  ctl_t * ctl,
  clim_t * clim,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/*! Calculate grid data for H2O2 chemistry module. */
void module_h2o2_chemgrid(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double t);

/*! KPP chemistry module. */
void module_kpp_chem(
  ctl_t * ctl,
  clim_t * clim,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/*! Calculate grid data for KPP chemistry module. */
void module_kpp_chemgrid(
  ctl_t * ctl,
  clim_t * clim,
  atm_t * atm,
  met_t * met0,
  met_t * met1,
  double t);

/*! Check position of air parcels. */
void module_position(
  ctl_t * ctl,
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

/*! Sort particles according to box index. */
void module_sort(
  ctl_t * ctl,
  met_t * met0,
  atm_t * atm);

/*! Helper function for sorting module. */
void module_sort_help(
  double *a,
  int *p,
  int np);

/*! Calculate time steps. */
void module_timesteps(
  ctl_t * ctl,
  atm_t * atm,
  met_t * met0,
  double *dt,
  double t);

/*! Initialize timesteps. */
void module_timesteps_init(
  ctl_t * ctl,
  atm_t * atm);

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

  clim_t *clim;

  met_t *met0, *met1;

#ifdef ASYNCIO
  met_t *met0TMP, *met1TMP, *mets;
  ctl_t ctlTMP;
#endif

  FILE *dirlist;

  char dirname[LEN], filename[2 * LEN];

  double *dt, *rs, t;

  int ntask = -1, rank = 0, size = 1;

  /* Start timers... */
  START_TIMERS;

  /* Initialize MPI... */
#ifdef MPI
  SELECT_TIMER("MPI_INIT", "INIT", NVTX_CPU);
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  /* Initialize GPUs... */
#ifdef _OPENACC
  SELECT_TIMER("ACC_INIT", "INIT", NVTX_GPU);
  if (acc_get_num_devices(acc_device_nvidia) <= 0)
    ERRMSG("Not running on a GPU device!");
  acc_set_device_num(rank % acc_get_num_devices(acc_device_nvidia),
		     acc_device_nvidia);
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
  while (fscanf(dirlist, "%4999s", dirname) != EOF) {

    /* MPI parallelization... */
    if ((++ntask) % size != rank)
      continue;

    /* ------------------------------------------------------------
       Initialize model run...
       ------------------------------------------------------------ */

    /* Allocate... */
    SELECT_TIMER("ALLOC", "MEMORY", NVTX_CPU);
    ALLOC(atm, atm_t, 1);
    ALLOC(cache, cache_t, 1);
    ALLOC(clim, clim_t, 1);
    ALLOC(met0, met_t, 1);
    ALLOC(met1, met_t, 1);
#ifdef ASYNCIO
    ALLOC(met0TMP, met_t, 1);
    ALLOC(met1TMP, met_t, 1);
#endif
    ALLOC(dt, double,
	  NP);
    ALLOC(rs, double,
	  3 * NP + 1);

    /* Create data region on GPUs... */
#ifdef _OPENACC
    SELECT_TIMER("CREATE_DATA_REGION", "MEMORY", NVTX_GPU);
#ifdef ASYNCIO
#pragma acc enter data create(atm[:1], cache[:1], clim[:1], ctl, ctlTMP, met0[:1], met1[:1], met0TMP[:1], met1TMP[:1], dt[:NP], rs[:3 * NP])
#else
#pragma acc enter data create(atm[:1], cache[:1], clim[:1], ctl, met0[:1], met1[:1], dt[:NP], rs[:3 * NP])
#endif
#endif

    /* Read control parameters... */
    sprintf(filename, "%s/%s", dirname, argv[2]);
    read_ctl(filename, argc, argv, &ctl);

    /* Read climatological data... */
    read_clim(&ctl, clim);

    /* Read atmospheric data... */
    sprintf(filename, "%s/%s", dirname, argv[3]);
    if (!read_atm(filename, &ctl, atm))
      ERRMSG("Cannot open file!");

    /* Initialize timesteps... */
    module_timesteps_init(&ctl, atm);

    /* Update GPU... */
#ifdef _OPENACC
    SELECT_TIMER("UPDATE_DEVICE", "MEMORY", NVTX_H2D);
#pragma acc update device(atm[:1], clim[:1], ctl)
#endif

    /* Initialize random number generator... */
    module_rng_init(ntask);

    /* Initialize meteo data... */
#ifdef ASYNCIO
    ctlTMP = ctl;
#endif
    get_met(&ctl, clim, ctl.t_start, &met0, &met1);
    if (ctl.dt_mod > fabs(met0->lon[1] - met0->lon[0]) * 111132. / 150.)
      WARN("Violation of CFL criterion! Check DT_MOD!");
#ifdef ASYNCIO
    get_met(&ctlTMP, clim, ctlTMP.t_start, &met0TMP, &met1TMP);
#endif

    /* Initialize isosurface... */
    if (ctl.isosurf >= 1 && ctl.isosurf <= 4)
      module_isosurf_init(&ctl, met0, met1, atm, cache);

    /* Update GPU... */
#ifdef _OPENACC
    SELECT_TIMER("UPDATE_DEVICE", "MEMORY", NVTX_H2D);
#pragma acc update device(cache[:1])
#endif

    /* ------------------------------------------------------------
       Loop over timesteps...
       ------------------------------------------------------------ */

    /* Loop over timesteps... */
#ifdef ASYNCIO
    omp_set_nested(1);
    // omp_set_dynamic(0);
    int ompTrdnum = omp_get_max_threads();
#endif
    for (t = ctl.t_start; ctl.direction * (t - ctl.t_stop) < ctl.dt_mod;
	 t += ctl.direction * ctl.dt_mod) {
#ifdef ASYNCIO
#pragma omp parallel num_threads(2)
      {
#endif

	/* Adjust length of final time step... */
	if (ctl.direction * (t - ctl.t_stop) > 0)
	  t = ctl.t_stop;

	/* Set time steps of air parcels... */
	module_timesteps(&ctl, atm, met0, dt, t);

	/* Get meteo data... */
#ifdef ASYNCIO
#pragma acc wait(5)
#pragma omp barrier
	if (omp_get_thread_num() == 0) {

	  /* Pointer swap... */
	  if (t != ctl.t_start) {
	    mets = met0;
	    met0 = met0TMP;
	    met0TMP = mets;

	    mets = met1;
	    met1 = met1TMP;
	    met1TMP = mets;
	  }
#endif
#ifndef ASYNCIO
	  if (t != ctl.t_start)
	    get_met(&ctl, clim, t, &met0, &met1);
#endif

	  /* Sort particles... */
	  if (ctl.sort_dt > 0 && fmod(t, ctl.sort_dt) == 0)
	    module_sort(&ctl, met0, atm);
	  
	  /* Check positions (initial)... */
	  module_position(&ctl, met0, met1, atm, dt);
	  
	  /* Advection... */
	  module_advect(&ctl, met0, met1, atm, dt);

	  /* Turbulent diffusion... */
	  if (ctl.turb_dx_trop > 0 || ctl.turb_dz_trop > 0
	      || ctl.turb_dx_strat > 0 || ctl.turb_dz_strat > 0)
	    module_diffusion_turb(&ctl, clim, atm, dt, rs);

	  /* Mesoscale diffusion... */
	  if (ctl.turb_mesox > 0 || ctl.turb_mesoz > 0)
	    module_diffusion_meso(&ctl, met0, met1, atm, cache, dt, rs);

	  /* Convection... */
	  if (ctl.conv_cape >= 0
	      && (ctl.conv_dt <= 0 || fmod(t, ctl.conv_dt) == 0))
	    module_convection(&ctl, met0, met1, atm, dt, rs);

	  /* Sedimentation... */
	  if (ctl.qnt_rp >= 0 && ctl.qnt_rhop >= 0)
	    module_sedi(&ctl, met0, met1, atm, dt);

	  /* Isosurface... */
	  if (ctl.isosurf >= 1 && ctl.isosurf <= 4)
	    module_isosurf(&ctl, met0, met1, atm, cache);

	  /* Check positions (final)... */
	  module_position(&ctl, met0, met1, atm, dt);

	  /* Interpolate meteo data... */
	  if (ctl.met_dt_out > 0
	      && (ctl.met_dt_out < ctl.dt_mod
		  || fmod(t, ctl.met_dt_out) == 0))
	    module_meteo(&ctl, clim, met0, met1, atm);
	  
	  /* Check boundary conditions (initial)... */
	  if (ctl.bound_mass >= 0 || ctl.bound_vmr >= 0
	      || ctl.kpp_chem_bound == 1)
	    module_bound_cond(&ctl, met0, met1, atm, dt);

	  /* Decay of particle mass... */
	  if (ctl.tdec_trop > 0 && ctl.tdec_strat > 0)
	    module_decay(&ctl, clim, atm, dt);

	  /* OH chemistry... */
	  if (ctl.clim_oh_filename[0] != '-' && ctl.oh_chem_reaction != 0)
	    module_oh_chem(&ctl, clim, met0, met1, atm, dt);

	  /* H2O2 chemistry (for SO2 aqueous phase oxidation)... */
	  if (ctl.clim_h2o2_filename[0] != '-' && ctl.h2o2_chem_reaction != 0) {
	    module_h2o2_chemgrid(&ctl, met0, met1, atm, t);
	    module_h2o2_chem(&ctl, clim, met0, met1, atm, dt);
	  }

#ifdef KPP
	  /* KPP chemistry... */
	  if (ctl.kpp_chem == 1) {
	    module_kpp_chemgrid(&ctl, clim, atm, met0, met1, t);
	    module_kpp_chem(&ctl, clim, met0, met1, atm, dt);
	  }
#endif

	  /* Wet deposition... */
	  if ((ctl.wet_depo_ic_a > 0 || ctl.wet_depo_ic_h[0] > 0)
	      && (ctl.wet_depo_bc_a > 0 || ctl.wet_depo_bc_h[0] > 0))
	    module_wet_deposition(&ctl, met0, met1, atm, dt);

	  /* Dry deposition... */
	  if (ctl.dry_depo_vdep > 0)
	    module_dry_deposition(&ctl, met0, met1, atm, dt);
	  
	  /* Check boundary conditions (final)... */
	  if (ctl.bound_mass >= 0 || ctl.bound_vmr >= 0
	      || ctl.kpp_chem_bound == 1)
	    module_bound_cond(&ctl, met0, met1, atm, dt);

	  /* Write output... */
	  write_output(dirname, &ctl, met0, met1, atm, t);
#ifdef ASYNCIO
	} else {
	  omp_set_num_threads(ompTrdnum);
	  if (ctl.direction * (t - ctl.t_stop + ctl.direction * ctl.dt_mod) <
	      ctl.dt_mod)
	    get_met(&ctl, clim, t + (ctl.direction * ctl.dt_mod), &met0TMP,
		    &met1TMP);
	}
      }
#endif
    }

#ifdef ASYNCIO
    omp_set_num_threads(ompTrdnum);
#endif

    /* ------------------------------------------------------------
       Finalize model run...
       ------------------------------------------------------------ */

    /* Report problem size... */
    LOG(1, "SIZE_NP = %d", atm->np);
    LOG(1, "SIZE_MPI_TASKS = %d", size);
    LOG(1, "SIZE_OMP_THREADS = %d", omp_get_max_threads());
#ifdef _OPENACC
    LOG(1, "SIZE_ACC_DEVICES = %d", acc_get_num_devices(acc_device_nvidia));
#endif

    /* Report memory usage... */
    LOG(1, "MEMORY_ATM = %g MByte", sizeof(atm_t) / 1024. / 1024.);
    LOG(1, "MEMORY_CACHE = %g MByte", sizeof(cache_t) / 1024. / 1024.);
    LOG(1, "MEMORY_CLIM = %g MByte", sizeof(clim_t) / 1024. / 1024.);
    LOG(1, "MEMORY_METEO = %g MByte", 2 * sizeof(met_t) / 1024. / 1024.);
    LOG(1, "MEMORY_DYNAMIC = %g MByte", (3 * NP * sizeof(int)
					 + 4 * NP * sizeof(double)
					 + EX * EY * EP * sizeof(float)) /
	1024. / 1024.);
    LOG(1, "MEMORY_STATIC = %g MByte", (EX * EY * EP * sizeof(float)) /
	1024. / 1024.);

    /* Delete data region on GPUs... */
#ifdef _OPENACC
    SELECT_TIMER("DELETE_DATA_REGION", "MEMORY", NVTX_GPU);
#ifdef ASYNCIO
#pragma acc exit data delete (ctl, atm, cache, clim, met0, met1, dt, rs, met0TMP, met1TMP)
#else
#pragma acc exit data delete (ctl, atm, cache, clim, met0, met1, dt, rs)
#endif
#endif

    /* Free... */
    SELECT_TIMER("FREE", "MEMORY", NVTX_CPU);
    free(atm);
    free(cache);
    free(clim);
    free(met0);
    free(met1);
#ifdef ASYNCIO
    free(met0TMP);
    free(met1TMP);
#endif
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

void module_advect(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt) {

  /* Set timer... */
  SELECT_TIMER("MODULE_ADVECTION", "PHYSICS", NVTX_GPU);

  const int np = atm->np;
#ifdef _OPENACC
#pragma acc data present(ctl,met0,met1,atm,dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < np; ip++)
    if (dt[ip] != 0) {

      /* Init... */
      double dts, u[4], um = 0, v[4], vm = 0, w[4], wm = 0, x[3];

      /* Loop over integration nodes... */
      for (int i = 0; i < ctl->advect; i++) {

	/* Set position... */
	if (i == 0) {
	  dts = 0.0;
	  x[0] = atm->lon[ip];
	  x[1] = atm->lat[ip];
	  x[2] = atm->p[ip];
	} else {
	  dts = (i == 3 ? 1.0 : 0.5) * dt[ip];
	  x[0] = atm->lon[ip] + DX2DEG(dts * u[i - 1] / 1000., atm->lat[ip]);
	  x[1] = atm->lat[ip] + DY2DEG(dts * v[i - 1] / 1000.);
	  x[2] = atm->p[ip] + dts * w[i - 1];
	}
	double tm = atm->time[ip] + dts;

	/* Interpolate meteo data... */
#ifdef UVW
	intpol_met_time_uvw(met0, met1, tm, x[2], x[0], x[1],
			    &u[i], &v[i], &w[i]);
#else
	INTPOL_INIT;
	intpol_met_time_3d(met0, met0->u, met1, met1->u, tm,
			   x[2], x[0], x[1], &u[i], ci, cw, 1);
	intpol_met_time_3d(met0, met0->v, met1, met1->v, tm,
			   x[2], x[0], x[1], &v[i], ci, cw, 0);
	intpol_met_time_3d(met0, met0->w, met1, met1->w, tm,
			   x[2], x[0], x[1], &w[i], ci, cw, 0);
#endif

	/* Get mean wind... */
	double k = 1.0;
	if (ctl->advect == 2)
	  k = (i == 0 ? 0.0 : 1.0);
	else if (ctl->advect == 4)
	  k = (i == 0 || i == 3 ? 1.0 / 6.0 : 2.0 / 6.0);
	um += k * u[i];
	vm += k * v[i];
	wm += k * w[i];
      }

      /* Set new position... */
      atm->time[ip] += dt[ip];
      atm->lon[ip] += DX2DEG(dt[ip] * um / 1000.,
			     (ctl->advect == 2 ? x[1] : atm->lat[ip]));
      atm->lat[ip] += DY2DEG(dt[ip] * vm / 1000.);
      atm->p[ip] += dt[ip] * wm;
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
  SELECT_TIMER("MODULE_BOUNDCOND", "PHYSICS", NVTX_GPU);

  /* Check quantity flags... */
  if (ctl->qnt_m < 0 && ctl->qnt_vmr < 0 && ctl->kpp_chem_bound == 0)
    ERRMSG("Module needs quantity mass or volume mixing ratio!");

  const int np = atm->np;
#ifdef _OPENACC
#pragma acc data present(ctl, met0, met1, atm, dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < np; ip++)
    if (dt[ip] != 0) {

      /* Check latitude and pressure range... */
      if (atm->lat[ip] < ctl->bound_lat0 || atm->lat[ip] > ctl->bound_lat1
	  || atm->p[ip] > ctl->bound_p0 || atm->p[ip] < ctl->bound_p1)
	continue;

      /* Check surface layer... */
      if (ctl->bound_dps > 0 || ctl->bound_dzs > 0
	  || ctl->bound_zetas > 0 || ctl->bound_pbl) {

	/* Get surface pressure... */
	double ps;
	INTPOL_INIT;
	INTPOL_2D(ps, 1);

	/* Check pressure... */
	if (ctl->bound_dps > 0 && atm->p[ip] < ps - ctl->bound_dps)
	  continue;

	/* Check height... */
	if (ctl->bound_dzs > 0 && Z(atm->p[ip]) > Z(ps) + ctl->bound_dzs)
	  continue;

	/* Check zeta range... */
	if (ctl->bound_zetas > 0) {
	  double t;
	  INTPOL_3D(t, 1);
	  if (ZETA(ps, atm->p[ip], t) > ctl->bound_zetas)
	    continue;
	}

	/* Check planetary boundary layer... */
	if (ctl->bound_pbl) {
	  double pbl;
	  INTPOL_2D(pbl, 0);
	  if (atm->p[ip] < pbl)
	    continue;
	}
      }

      /* Set mass and volume mixing ratio... */
      if (ctl->qnt_m >= 0 && ctl->bound_mass >= 0)
	atm->q[ctl->qnt_m][ip] =
	  ctl->bound_mass + ctl->bound_mass_trend * atm->time[ip];
      if (ctl->qnt_vmr >= 0 && ctl->bound_vmr >= 0)
	atm->q[ctl->qnt_vmr][ip] =
	  ctl->bound_vmr + ctl->bound_vmr_trend * atm->time[ip];
#ifdef KPP
      if (ctl->kpp_chem_bound == 1)
	kpp_chem_bound_cond(ctl, atm, met0, met1, ip);
#endif
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
  SELECT_TIMER("MODULE_CONVECTION", "PHYSICS", NVTX_GPU);

  /* Create random numbers... */
  module_rng(rs, (size_t) atm->np, 0);

  const int np = atm->np;
#ifdef _OPENACC
#pragma acc data present(ctl, met0, met1, atm, dt, rs)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < np; ip++)
    if (dt[ip] != 0) {

      double cape, cin, pel, ps;

      /* Interpolate CAPE... */
      INTPOL_INIT;
      INTPOL_2D(cape, 1);

      /* Check threshold... */
      if (isfinite(cape) && cape >= ctl->conv_cape) {

	/* Check CIN... */
	if (ctl->conv_cin > 0) {
	  INTPOL_2D(cin, 0);
	  if (isfinite(cin) && cin >= ctl->conv_cin)
	    continue;
	}

	/* Interpolate equilibrium level... */
	INTPOL_2D(pel, 0);

	/* Check whether particle is above cloud top... */
	if (!isfinite(pel) || atm->p[ip] < pel)
	  continue;

	/* Set pressure range for mixing... */
	double pbot = atm->p[ip];
	double ptop = atm->p[ip];
	if (ctl->conv_mix_bot == 1) {
	  INTPOL_2D(ps, 0);
	  pbot = ps;
	}
	if (ctl->conv_mix_top == 1)
	  ptop = pel;

	/* Vertical mixing based on pressure... */
	if (ctl->conv_mix == 0)
	  atm->p[ip] = pbot + (ptop - pbot) * rs[ip];

	/* Vertical mixing based on density... */
	else if (ctl->conv_mix == 1) {

	  /* Get density range... */
	  double tbot, ttop;
	  intpol_met_time_3d(met0, met0->t, met1, met1->t, atm->time[ip],
			     pbot, atm->lon[ip], atm->lat[ip], &tbot,
			     ci, cw, 1);
	  intpol_met_time_3d(met0, met0->t, met1, met1->t, atm->time[ip],
			     ptop, atm->lon[ip], atm->lat[ip], &ttop,
			     ci, cw, 1);
	  double rhobot = pbot / tbot;
	  double rhotop = ptop / ttop;

	  /* Get new density... */
	  double lrho = log(rhobot + (rhotop - rhobot) * rs[ip]);

	  /* Find pressure... */
	  double lrhobot = log(rhobot);
	  double lrhotop = log(rhotop);
	  double lpbot = log(pbot);
	  double lptop = log(ptop);
	  atm->p[ip] = exp(LIN(lrhobot, lpbot, lrhotop, lptop, lrho));
	}
      }
    }
}

/*****************************************************************************/

void module_decay(
  ctl_t * ctl,
  clim_t * clim,
  atm_t * atm,
  double *dt) {

  /* Set timer... */
  SELECT_TIMER("MODULE_DECAY", "PHYSICS", NVTX_GPU);

  /* Check quantity flags... */
  if (ctl->qnt_m < 0 && ctl->qnt_vmr < 0)
    ERRMSG("Module needs quantity mass or volume mixing ratio!");

  const int np = atm->np;
#ifdef _OPENACC
#pragma acc data present(ctl,clim,atm,dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < np; ip++)
    if (dt[ip] != 0) {

      /* Get weighting factor... */
      double w = tropo_weight(clim, atm->time[ip], atm->lat[ip], atm->p[ip]);

      /* Set lifetime... */
      double tdec = w * ctl->tdec_trop + (1 - w) * ctl->tdec_strat;

      /* Calculate exponential decay... */
      double aux = exp(-dt[ip] / tdec);
      if (ctl->qnt_m >= 0) {
	if (ctl->qnt_mloss_decay >= 0)
	  atm->q[ctl->qnt_mloss_decay][ip]
	    += atm->q[ctl->qnt_m][ip] * (1 - aux);
	atm->q[ctl->qnt_m][ip] *= aux;
      }
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
  SELECT_TIMER("MODULE_TURBMESO", "PHYSICS", NVTX_GPU);

  /* Create random numbers... */
  module_rng(rs, 3 * (size_t) atm->np, 1);

  const int np = atm->np;
#ifdef _OPENACC
#pragma acc data present(ctl, met0, met1, atm, cache, dt, rs)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < np; ip++)
    if (dt[ip] != 0) {

      /* Get indices... */
      int ix = locate_reg(met0->lon, met0->nx, atm->lon[ip]);
      int iy = locate_reg(met0->lat, met0->ny, atm->lat[ip]);
      int iz = locate_irr(met0->p, met0->np, atm->p[ip]);

      /* Get standard deviations of local wind data... */
      float umean = 0, usig = 0, vmean = 0, vsig = 0, wmean = 0, wsig = 0;
      for (int i = 0; i < 2; i++)
	for (int j = 0; j < 2; j++)
	  for (int k = 0; k < 2; k++) {
#ifdef UVW
	    umean += met0->uvw[ix + i][iy + j][iz + k][0];
	    usig += SQR(met0->uvw[ix + i][iy + j][iz + k][0]);
	    vmean += met0->uvw[ix + i][iy + j][iz + k][1];
	    vsig += SQR(met0->uvw[ix + i][iy + j][iz + k][1]);
	    wmean += met0->uvw[ix + i][iy + j][iz + k][2];
	    wsig += SQR(met0->uvw[ix + i][iy + j][iz + k][2]);

	    umean += met1->uvw[ix + i][iy + j][iz + k][0];
	    usig += SQR(met1->uvw[ix + i][iy + j][iz + k][0]);
	    vmean += met1->uvw[ix + i][iy + j][iz + k][1];
	    vsig += SQR(met1->uvw[ix + i][iy + j][iz + k][1]);
	    wmean += met1->uvw[ix + i][iy + j][iz + k][2];
	    wsig += SQR(met1->uvw[ix + i][iy + j][iz + k][2]);
#else
	    umean += met0->u[ix + i][iy + j][iz + k];
	    usig += SQR(met0->u[ix + i][iy + j][iz + k]);
	    vmean += met0->v[ix + i][iy + j][iz + k];
	    vsig += SQR(met0->v[ix + i][iy + j][iz + k]);
	    wmean += met0->w[ix + i][iy + j][iz + k];
	    wsig += SQR(met0->w[ix + i][iy + j][iz + k]);

	    umean += met1->u[ix + i][iy + j][iz + k];
	    usig += SQR(met1->u[ix + i][iy + j][iz + k]);
	    vmean += met1->v[ix + i][iy + j][iz + k];
	    vsig += SQR(met1->v[ix + i][iy + j][iz + k]);
	    wmean += met1->w[ix + i][iy + j][iz + k];
	    wsig += SQR(met1->w[ix + i][iy + j][iz + k]);
#endif
	  }
      usig = usig / 16.f - SQR(umean / 16.f);
      usig = (usig > 0 ? sqrtf(usig) : 0);
      vsig = vsig / 16.f - SQR(vmean / 16.f);
      vsig = (vsig > 0 ? sqrtf(vsig) : 0);
      wsig = wsig / 16.f - SQR(wmean / 16.f);
      wsig = (wsig > 0 ? sqrtf(wsig) : 0);

      /* Set temporal correlations for mesoscale fluctuations... */
      double r = 1 - 2 * fabs(dt[ip]) / ctl->dt_met;
      double r2 = sqrt(1 - r * r);

      /* Calculate horizontal mesoscale wind fluctuations... */
      if (ctl->turb_mesox > 0) {
	cache->uvwp[ip][0] =
	  (float) (r * cache->uvwp[ip][0] +
		   r2 * rs[3 * ip] * ctl->turb_mesox * usig);
	atm->lon[ip] +=
	  DX2DEG(cache->uvwp[ip][0] * dt[ip] / 1000., atm->lat[ip]);

	cache->uvwp[ip][1] =
	  (float) (r * cache->uvwp[ip][1] +
		   r2 * rs[3 * ip + 1] * ctl->turb_mesox * vsig);
	atm->lat[ip] += DY2DEG(cache->uvwp[ip][1] * dt[ip] / 1000.);
      }

      /* Calculate vertical mesoscale wind fluctuations... */
      if (ctl->turb_mesoz > 0) {
	cache->uvwp[ip][2] =
	  (float) (r * cache->uvwp[ip][2] +
		   r2 * rs[3 * ip + 2] * ctl->turb_mesoz * wsig);
	atm->p[ip] += cache->uvwp[ip][2] * dt[ip];
      }
    }
}

/*****************************************************************************/

void module_diffusion_turb(
  ctl_t * ctl,
  clim_t * clim,
  atm_t * atm,
  double *dt,
  double *rs) {

  /* Set timer... */
  SELECT_TIMER("MODULE_TURBDIFF", "PHYSICS", NVTX_GPU);

  /* Create random numbers... */
  module_rng(rs, 3 * (size_t) atm->np, 1);

  const int np = atm->np;
#ifdef _OPENACC
#pragma acc data present(ctl,clim,atm,dt,rs)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < np; ip++)
    if (dt[ip] != 0) {

      /* Get weighting factor... */
      double w = tropo_weight(clim, atm->time[ip], atm->lat[ip], atm->p[ip]);

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
	atm->p[ip] += DZ2DP(rs[3 * ip + 2] * sigma / 1000., atm->p[ip]);
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
  SELECT_TIMER("MODULE_DRYDEPO", "PHYSICS", NVTX_GPU);

  /* Check quantity flags... */
  if (ctl->qnt_m < 0 && ctl->qnt_vmr < 0)
    ERRMSG("Module needs quantity mass or volume mixing ratio!");

  const int np = atm->np;
#ifdef _OPENACC
#pragma acc data present(ctl, met0, met1, atm, dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < np; ip++)
    if (dt[ip] != 0) {

      double ps, t, v_dep;

      /* Get surface pressure... */
      INTPOL_INIT;
      INTPOL_2D(ps, 1);

      /* Check whether particle is above the surface layer... */
      if (atm->p[ip] < ps - ctl->dry_depo_dp)
	continue;

      /* Set depth of surface layer... */
      double dz = 1000. * (Z(ps - ctl->dry_depo_dp) - Z(ps));

      /* Calculate sedimentation velocity for particles... */
      if (ctl->qnt_rp > 0 && ctl->qnt_rhop > 0) {

	/* Get temperature... */
	INTPOL_3D(t, 1);

	/* Set deposition velocity... */
	v_dep = sedi(atm->p[ip], t, atm->q[ctl->qnt_rp][ip],
		     atm->q[ctl->qnt_rhop][ip]);
      }

      /* Use explicit sedimentation velocity for gases... */
      else
	v_dep = ctl->dry_depo_vdep;

      /* Calculate loss of mass based on deposition velocity... */
      double aux = exp(-dt[ip] * v_dep / dz);
      if (ctl->qnt_m >= 0) {
	if (ctl->qnt_mloss_dry >= 0)
	  atm->q[ctl->qnt_mloss_dry][ip]
	    += atm->q[ctl->qnt_m][ip] * (1 - aux);
	atm->q[ctl->qnt_m][ip] *= aux;
      }
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
  SELECT_TIMER("MODULE_ISOSURF", "PHYSICS", NVTX_GPU);

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
  SELECT_TIMER("MODULE_ISOSURF", "PHYSICS", NVTX_GPU);

  const int np = atm->np;
#ifdef _OPENACC
#pragma acc data present(ctl, met0, met1, atm, cache)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < np; ip++) {

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
  clim_t * clim,
  met_t * met0,
  met_t * met1,
  atm_t * atm) {

  /* Set timer... */
  SELECT_TIMER("MODULE_METEO", "PHYSICS", NVTX_GPU);

  /* Check quantity flags... */
  if (ctl->qnt_tsts >= 0)
    if (ctl->qnt_tice < 0 || ctl->qnt_tnat < 0)
      ERRMSG("Need T_ice and T_NAT to calculate T_STS!");

  const int np = atm->np;
#ifdef _OPENACC
#pragma acc data present(ctl, clim, met0, met1, atm)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < np; ip++) {

    double ps, ts, zs, us, vs, lsm, sst, pbl, pt, pct, pcb, cl, plcl, plfc,
      pel, cape, cin, pv, t, tt, u, v, w, h2o, h2ot, o3, lwc, iwc, cc, z, zt;

    /* Interpolate meteo data... */
    INTPOL_INIT;
    INTPOL_TIME_ALL(atm->time[ip], atm->p[ip], atm->lon[ip], atm->lat[ip]);

    /* Set quantities... */
    SET_ATM(qnt_ps, ps);
    SET_ATM(qnt_ts, ts);
    SET_ATM(qnt_zs, zs);
    SET_ATM(qnt_us, us);
    SET_ATM(qnt_vs, vs);
    SET_ATM(qnt_lsm, lsm);
    SET_ATM(qnt_sst, sst);
    SET_ATM(qnt_pbl, pbl);
    SET_ATM(qnt_pt, pt);
    SET_ATM(qnt_tt, tt);
    SET_ATM(qnt_zt, zt);
    SET_ATM(qnt_h2ot, h2ot);
    SET_ATM(qnt_z, z);
    SET_ATM(qnt_p, atm->p[ip]);
    SET_ATM(qnt_t, t);
    SET_ATM(qnt_rho, RHO(atm->p[ip], t));
    SET_ATM(qnt_u, u);
    SET_ATM(qnt_v, v);
    SET_ATM(qnt_w, w);
    SET_ATM(qnt_h2o, h2o);
    SET_ATM(qnt_o3, o3);
    SET_ATM(qnt_lwc, lwc);
    SET_ATM(qnt_iwc, iwc);
    SET_ATM(qnt_cc, cc);
    SET_ATM(qnt_pct, pct);
    SET_ATM(qnt_pcb, pcb);
    SET_ATM(qnt_cl, cl);
    SET_ATM(qnt_plcl, plcl);
    SET_ATM(qnt_plfc, plfc);
    SET_ATM(qnt_pel, pel);
    SET_ATM(qnt_cape, cape);
    SET_ATM(qnt_cin, cin);
    SET_ATM(qnt_hno3,
	    clim_hno3(clim, atm->time[ip], atm->lat[ip], atm->p[ip]));
    SET_ATM(qnt_oh, clim_oh(ctl, clim, atm->time[ip],
			    atm->lon[ip], atm->lat[ip], atm->p[ip]));
    SET_ATM(qnt_h2o2, clim_var(&clim->h2o2, atm->time[ip],
			       atm->lat[ip], atm->p[ip]));
    SET_ATM(qnt_ho2, clim_var(&clim->ho2, atm->time[ip],
			      atm->lat[ip], atm->p[ip]));
    SET_ATM(qnt_o1d, clim_var(&clim->o1d, atm->time[ip],
			      atm->lat[ip], atm->p[ip]));
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
			    clim_hno3(clim, atm->time[ip], atm->lat[ip],
				      atm->p[ip])));
    SET_ATM(qnt_tsts,
	    0.5 * (atm->q[ctl->qnt_tice][ip] + atm->q[ctl->qnt_tnat][ip]));
  }
}

/*****************************************************************************/

void module_oh_chem(
  ctl_t * ctl,
  clim_t * clim,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt) {

  /* Set timer... */
  SELECT_TIMER("MODULE_OHCHEM", "PHYSICS", NVTX_GPU);

  /* Check quantity flags... */
  if (ctl->qnt_m < 0 && ctl->qnt_vmr < 0)
    ERRMSG("Module needs quantity mass or volume mixing ratio!");

  const int np = atm->np;
#ifdef _OPENACC
#pragma acc data present(ctl, clim, met0, met1, atm, dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < np; ip++)
    if (dt[ip] != 0) {

      /* Get temperature... */
      double t;
      INTPOL_INIT;
      INTPOL_3D(t, 1);

      /* Use constant reaction rate... */
      double k = GSL_NAN;
      if (ctl->oh_chem_reaction == 1)
	k = ctl->oh_chem[0];

      /* Calculate bimolecular reaction rate... */
      else if (ctl->oh_chem_reaction == 2)
	k = ctl->oh_chem[0] * exp(-ctl->oh_chem[1] / t);

      /* Calculate termolecular reaction rate... */
      if (ctl->oh_chem_reaction == 3) {

	/* Calculate molecular density (IUPAC Data Sheet I.A4.86 SOx15)... */
	double M = MOLEC_DEN_AIR(atm->p[ip], t);

	/* Calculate rate coefficient for X + OH + M -> XOH + M
	   (JPL Publication 19-05) ... */
	double k0 =
	  ctl->oh_chem[0] * (ctl->oh_chem[1] !=
			     0 ? pow(298. / t, ctl->oh_chem[1]) : 1.);
	double ki =
	  ctl->oh_chem[2] * (ctl->oh_chem[3] !=
			     0 ? pow(298. / t, ctl->oh_chem[3]) : 1.);
	double c = log10(k0 * M / ki);
	k = k0 * M / (1. + k0 * M / ki) * pow(0.6, 1. / (1. + c * c));
      }

      /* Calculate exponential decay... */
      double rate_coef = k * clim_oh(ctl, clim, atm->time[ip], atm->lon[ip],
				     atm->lat[ip], atm->p[ip]);
      double aux = exp(-dt[ip] * rate_coef);
      if (ctl->qnt_m >= 0) {
	if (ctl->qnt_mloss_oh >= 0)
	  atm->q[ctl->qnt_mloss_oh][ip]
	    += atm->q[ctl->qnt_m][ip] * (1 - aux);
	atm->q[ctl->qnt_m][ip] *= aux;
      }
      if (ctl->qnt_vmr >= 0)
	atm->q[ctl->qnt_vmr][ip] *= aux;
    }
}

/*****************************************************************************/

void module_h2o2_chem(
  ctl_t * ctl,
  clim_t * clim,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt) {

  /* Set timer... */
  SELECT_TIMER("MODULE_H2O2CHEM", "PHYSICS", NVTX_GPU);

  /* Check quantity flags... */
  if (ctl->qnt_m < 0 && ctl->qnt_vmr < 0)
    ERRMSG("Module needs quantity mass or volume mixing ratio!");
  if (ctl->qnt_vmrimpl < 0)
    ERRMSG("Module needs quantity implicit volume mixing ratio!");

  const int np = atm->np;
#ifdef _OPENACC
#pragma acc data present(clim,ctl,met0,met1,atm,dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < np; ip++)
    if (dt[ip] != 0) {

      /* Check whether particle is inside cloud... */
      double lwc;
      INTPOL_INIT;
      INTPOL_3D(lwc, 1);
      if (!(lwc > 0))
	continue;

      /* Check implicit volume mixing ratio... */
      if (atm->q[ctl->qnt_vmrimpl][ip] == 0)
	continue;

      /* Get temperature... */
      double t;
      INTPOL_3D(t, 0);

      /* Reaction rate (Berglen et al., 2004)... */
      double k = 9.1e7 * exp(-29700 / RI * (1. / t - 1. / 298.15));	// Maass  1999 unit: M^(-2)

      /* Henry constant of SO2... */
      double H_SO2 = 1.3e-2 * exp(2900 * (1. / t - 1. / 298.15)) * RI * t;
      double K_1S = 1.23e-2 * exp(2.01e3 * (1. / t - 1. / 298.15));	// unit: M

      /* Henry constant of H2O2... */
      double H_h2o2 = 8.3e2 * exp(7600 * (1 / t - 1 / 298.15)) * RI * t;

      /* Concentration of H2O2 (Barth et al., 1989)... */
      double SO2 = atm->q[ctl->qnt_vmrimpl][ip] * 1e9;	// vmr unit: ppbv
      double h2o2 = H_h2o2
	* clim_var(&clim->h2o2, atm->time[ip], atm->lat[ip], atm->p[ip])
	* 0.59 * exp(-0.687 * SO2) * 1000 / AVO;	// unit: M

      /* Volume water content in cloud [m^3 m^(-3)]... */
      double rho_air = 100 * atm->p[ip] / (RI * t) * MA / 1000;
      //MA: Molar mass of dry air; RI: Ideal gas constant 8.314 [J/(mol K)]
      double CWC = lwc * rho_air / 1000;

      /* Calculate exponential decay (Rolph et al., 1992)... */
      double rate_coef = k * K_1S * h2o2 * H_SO2 * CWC;
      double aux = exp(-dt[ip] * rate_coef);
      if (ctl->qnt_m >= 0) {
	if (ctl->qnt_mloss_h2o2 >= 0)
	  atm->q[ctl->qnt_mloss_h2o2][ip] +=
	    atm->q[ctl->qnt_m][ip] * (1 - aux);
	atm->q[ctl->qnt_m][ip] *= aux;
      }
      if (ctl->qnt_vmr >= 0)
	atm->q[ctl->qnt_vmr][ip] *= aux;
    }
}

/*****************************************************************************/

void module_h2o2_chemgrid(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double t) {

  double *mass, *z, *lon, *lat, *press, *area;

  int *ixs, *iys, *izs;

  /* Update host... */
#ifdef _OPENACC
  SELECT_TIMER("UPDATE_HOST", "MEMORY", NVTX_D2H);
#pragma acc update host(atm[:1])
#endif

  /* Set timer... */
  SELECT_TIMER("MODULE_H2O2_CHEMGRID", "PHYSICS", NVTX_GPU);

  /* Check quantity flags... */
  if (ctl->qnt_m < 0)
    ERRMSG("Module needs quantity mass!");
  if (ctl->molmass <= 0)
    ERRMSG("Specify molar mass!");

  /* Allocate... */
  ALLOC(mass, double,
	ctl->chemgrid_nx * ctl->chemgrid_ny * ctl->chemgrid_nz);
  ALLOC(z, double,
	ctl->chemgrid_nz);
  ALLOC(lon, double,
	ctl->chemgrid_nx);
  ALLOC(lat, double,
	ctl->chemgrid_ny);
  ALLOC(area, double,
	ctl->chemgrid_ny);
  ALLOC(press, double,
	ctl->chemgrid_nz);
  ALLOC(ixs, int,
	atm->np);
  ALLOC(iys, int,
	atm->np);
  ALLOC(izs, int,
	atm->np);

  /* Set grid box size... */
  double dz = (ctl->chemgrid_z1 - ctl->chemgrid_z0) / ctl->chemgrid_nz;
  double dlon = (ctl->chemgrid_lon1 - ctl->chemgrid_lon0) / ctl->chemgrid_nx;
  double dlat = (ctl->chemgrid_lat1 - ctl->chemgrid_lat0) / ctl->chemgrid_ny;

  /* Set vertical coordinates... */
#pragma omp parallel for default(shared)
  for (int iz = 0; iz < ctl->chemgrid_nz; iz++) {
    z[iz] = ctl->chemgrid_z0 + dz * (iz + 0.5);
    press[iz] = P(z[iz]);
  }

  /* Set horizontal coordinates... */
  for (int ix = 0; ix < ctl->chemgrid_nx; ix++)
    lon[ix] = ctl->chemgrid_lon0 + dlon * (ix + 0.5);
#pragma omp parallel for default(shared)
  for (int iy = 0; iy < ctl->chemgrid_ny; iy++) {
    lat[iy] = ctl->chemgrid_lat0 + dlat * (iy + 0.5);
    area[iy] = dlat * dlon * SQR(RE * M_PI / 180.)
      * cos(lat[iy] * M_PI / 180.);
  }

  /* Set time interval for output... */
  double t0 = t - 0.5 * ctl->dt_mod;
  double t1 = t + 0.5 * ctl->dt_mod;

  /* Get indices... */
#pragma omp parallel for default(shared)
  for (int ip = 0; ip < atm->np; ip++) {
    ixs[ip] = (int) ((atm->lon[ip] - ctl->chemgrid_lon0) / dlon);
    iys[ip] = (int) ((atm->lat[ip] - ctl->chemgrid_lat0) / dlat);
    izs[ip] = (int) ((Z(atm->p[ip]) - ctl->chemgrid_z0) / dz);
    if (atm->time[ip] < t0 || atm->time[ip] > t1
	|| ixs[ip] < 0 || ixs[ip] >= ctl->chemgrid_nx
	|| iys[ip] < 0 || iys[ip] >= ctl->chemgrid_ny
	|| izs[ip] < 0 || izs[ip] >= ctl->chemgrid_nz)
      izs[ip] = -1;
  }

  /* Average data... */
  for (int ip = 0; ip < atm->np; ip++)
    if (izs[ip] >= 0)
      mass[ARRAY_3D
	   (ixs[ip], iys[ip], ctl->chemgrid_ny, izs[ip], ctl->chemgrid_nz)]
	+= atm->q[ctl->qnt_m][ip];

  /* Assign the grid data to air parcels ... */
#pragma omp parallel for default(shared)
  for (int ip = 0; ip < atm->np; ip++)
    if (izs[ip] >= 0) {

      if (ctl->qnt_vmrimpl > 0) {

	/* Inteprolate temperature... */
	double temp;
	INTPOL_INIT;
	intpol_met_time_3d(met0, met0->t, met1, met1->t, t, press[izs[ip]],
			   lon[ixs[ip]], lat[iys[ip]], &temp, ci, cw, 1);

	atm->q[ctl->qnt_vmrimpl][ip] = MA / ctl->molmass *
	  mass[ARRAY_3D
	       (ixs[ip], iys[ip], ctl->chemgrid_ny, izs[ip],
		ctl->chemgrid_nz)]
	  / (RHO(press[izs[ip]], temp) * 1e6 * area[iys[ip]] * 1e3 * dz);
      }
    }

  /* Free... */
  free(mass);
  free(z);
  free(lon);
  free(lat);
  free(area);
  free(press);
  free(ixs);
  free(iys);
  free(izs);

  /* Update device... */
#ifdef _OPENACC
  SELECT_TIMER("UPDATE_DEVICE", "MEMORY", NVTX_H2D);
#pragma acc update device(atm[:1])
#endif
}

/*****************************************************************************/

#ifdef KPP

void INTEGRATE(
  double TIN,
  double TOUT);

/*****************************************************************************/

void kpp_chemgrid_mass2concen(
  atm_t * atm,
  ctl_t * ctl,
  double *mass,
  double *area,
  int *ixs,
  int *iys,
  int *izs,
  double dz,
  int ip,
  int qnt_index) {

  if (qnt_index > 0)		// TODO: skip check of qnt_index in this function?
    atm->q[qnt_index][ip]
      = AVO * mass[ARRAY_3D(ixs[ip], iys[ip], ctl->chemgrid_ny,
			    izs[ip], ctl->chemgrid_nz)]
      / (1e18 * area[iys[ip]] * dz * ctl->molmass);	//Unit: molec/cm3
  else
    ERRMSG("Some quantity variables are not defined!");
}

/*****************************************************************************/

double param_mixing_calc(
  ctl_t * ctl,
  clim_t * clim,
  atm_t * atm,
  int ip) {

  /* Check control parameters... */// TODO: add check for mixparam_* < 0 ?
  if (ctl->chemgrid_mixparam_trop < 1 || ctl->chemgrid_mixparam_strat < 1) {

    /* Get weighting factor... */
    double w = tropo_weight(clim, atm->time[ip], atm->lat[ip], atm->p[ip]);

    /* Set interparcel exchange parameter (Collins et al. 1997)... */
    return w * ctl->chemgrid_mixparam_trop
      + (1 - w) * ctl->chemgrid_mixparam_strat;
  }

  /* Complete mixing... */
  else
    return 1;
}

/*****************************************************************************/

// TODO: it is a bit confusing that we have a function interparc_mixing_help() here, but not interparc_mixing() ?

void interparc_mixing_help(
  ctl_t * ctl,
  atm_t * atm,
  clim_t * clim,
  int *ixs,
  int *iys,
  int *izs,
  int qnt_idx) {

  /* Check quantity flag... */
  if (qnt_idx < 0)
    return;

  /* Allocate... */
  double *cmean;
  int *count;
  ALLOC(cmean, double,
	ctl->chemgrid_nx * ctl->chemgrid_ny * ctl->chemgrid_nz);
  ALLOC(count, int,
	ctl->chemgrid_nx * ctl->chemgrid_ny * ctl->chemgrid_nz);

  /* Loop over particles... */
  for (int ip = 0; ip < atm->np; ip++)
    if (izs[ip] >= 0) {
      cmean[ARRAY_3D
	    (ixs[ip], iys[ip], ctl->chemgrid_ny, izs[ip], ctl->chemgrid_nz)]
	+= atm->q[qnt_idx][ip];
      count[ARRAY_3D
	    (ixs[ip], iys[ip], ctl->chemgrid_ny, izs[ip],
	     ctl->chemgrid_nz)] += 1;
    }
  for (int i = 0;
       i < ctl->chemgrid_nx * ctl->chemgrid_ny * ctl->chemgrid_nz; i++)
    if (count[i] > 0)
      cmean[i] /= count[i];

#pragma omp parallel for
  for (int ip = 0; ip < atm->np; ip++)
    if (izs[ip] >= 0) {

      atm->q[qnt_idx][ip] +=
	(cmean
	 [ARRAY_3D
	  (ixs[ip], iys[ip], ctl->chemgrid_ny, izs[ip], ctl->chemgrid_nz)]
	 - atm->q[qnt_idx][ip]) * param_mixing_calc(ctl, clim, atm, ip);
    }

  /* Free... */
  free(cmean);
  free(count);
}

/*****************************************************************************/

void module_kpp_chemgrid(
  ctl_t * ctl,
  clim_t * clim,
  atm_t * atm,
  met_t * met0,
  met_t * met1,
  double t) {

  int *ixs, *iys, *izs;

  /* Update host... */
#ifdef _OPENACC
  SELECT_TIMER("UPDATE_HOST", "MEMORY", NVTX_D2H);
#pragma acc update host(atm[:1])
#endif

  /* Set timer... */
  SELECT_TIMER("MODULE_KPP_CHEMGRID", "PHYSICS", NVTX_GPU);

  /* Allocate... */
  ALLOC(ixs, int,
	atm->np);
  ALLOC(iys, int,
	atm->np);
  ALLOC(izs, int,
	atm->np);

  /* Set grid box size... */
  double dz = (ctl->chemgrid_z1 - ctl->chemgrid_z0) / ctl->chemgrid_nz;
  double dlon = (ctl->chemgrid_lon1 - ctl->chemgrid_lon0) / ctl->chemgrid_nx;
  double dlat = (ctl->chemgrid_lat1 - ctl->chemgrid_lat0) / ctl->chemgrid_ny;

  /* Set time interval for output... */
  double t0 = t - 0.5 * ctl->dt_mod;
  double t1 = t + 0.5 * ctl->dt_mod;

  /* Get indices... */
#pragma omp parallel for default(shared)
  for (int ip = 0; ip < atm->np; ip++) {
    ixs[ip] = (int) ((atm->lon[ip] - ctl->chemgrid_lon0) / dlon);
    iys[ip] = (int) ((atm->lat[ip] - ctl->chemgrid_lat0) / dlat);
    izs[ip] = (int) ((Z(atm->p[ip]) - ctl->chemgrid_z0) / dz);
    if (atm->time[ip] < t0 || atm->time[ip] > t1
	|| ixs[ip] < 0 || ixs[ip] >= ctl->chemgrid_nx
	|| iys[ip] < 0 || iys[ip] >= ctl->chemgrid_ny
	|| izs[ip] < 0 || izs[ip] >= ctl->chemgrid_nz)
      izs[ip] = -1;
  }

  /* Initialize quantity concentration variables... */
  static int init[NP] = { 0 };	// TODO: use short rather than int? static variables are initilaized to zero by default
#pragma omp parallel for default(shared)
  for (int ip = 0; ip < atm->np; ip++)
    if (izs[ip] >= 0)
      if (!init[ip]) {
	init[ip] += 1;		// TODO: why += 1 ? why not init[ip] = 1 ?
	kpp_chem_init_cqnt(ctl, atm, clim, met0, met1, ip);
      }

  /* Calculate trace species concentration according to mass data... */
  if (ctl->qnt_m >= 0) {	// TODO: move these checks close to the begin of the function?
    if (ctl->molmass <= 0)
      ERRMSG("Specify molar mass!");

    /* Allocate... */
    double *mass, *area, *lat;
    ALLOC(mass, double,
	  ctl->chemgrid_nx * ctl->chemgrid_ny * ctl->chemgrid_nz);
    ALLOC(lat, double,
	  ctl->chemgrid_ny);
    ALLOC(area, double,
	  ctl->chemgrid_ny);

    /* Calculate latitude and area... */
#pragma omp parallel for default(shared)
    for (int iy = 0; iy < ctl->chemgrid_ny; iy++) {
      lat[iy] = ctl->chemgrid_lat0 + dlat * (iy + 0.5);
      area[iy] = dlat * dlon * SQR(RE * M_PI / 180.)
	* cos(lat[iy] * M_PI / 180.);
    }

    /* Calculate mass per grid box... */
    for (int ip = 0; ip < atm->np; ip++)
      if (izs[ip] >= 0)
	mass[ARRAY_3D
	     (ixs[ip], iys[ip], ctl->chemgrid_ny, izs[ip], ctl->chemgrid_nz)]
	  += atm->q[ctl->qnt_m][ip];

    /* Convert mass to concentration... */
#pragma omp parallel for default(shared)
    for (int ip = 0; ip < atm->np; ip++)
      if (izs[ip] >= 0)
	if (ctl->qnt_m >= 0)
	  kpp_chemgrid_mass2concen(atm, ctl, mass, area, ixs, iys, izs,
				   dz, ip, ctl->qnt_Cx);

    /* Free... */
    free(mass);
    free(area);
    free(lat);
  }

  /* Calculate the inter-parcel exchange between parcel and backgroud... */
  interparc_mixing(ctl, atm, clim, ixs, iys, izs);

  /* Free... */
  free(ixs);
  free(iys);
  free(izs);

  /* Update device... */
#ifdef _OPENACC
  SELECT_TIMER("UPDATE_DEVICE", "MEMORY", NVTX_H2D);
#pragma acc update device(atm[:1])
#endif
}

/*****************************************************************************/

void module_kpp_chem(
  ctl_t * ctl,
  clim_t * clim,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt) {

  /* Set timer... */
  SELECT_TIMER("MODULE_KPP_CHEM", "PHYSICS", NVTX_GPU);

  /* Loop over particles... */
  const int np = atm->np;
#pragma omp parallel for
  for (int ip = 0; ip < np; ip++) {
    if (dt[ip] > 0) {

      /* Set species... */
      ALLOC(VAR, double,
	    NVAR);
      ALLOC(FIX, double,
	    NFIX);

      STEPMIN = 0;		// TODO: make STEPMIN and STEPMAX control parameters (ctl->kpp_stepmin, ctl->kpp_stepmax)?
      STEPMAX = 900.0;

      /* Set relative & absolute tolerances... */
      for (int i = 0; i < NVAR; i++) {
	RTOL[i] = 1.0e-4;	// TODO: make ATOL and RTOL control parameters  (ctl->kpp_atol, ctl->kpp_rtol)?
	ATOL[i] = 1.0e-3;
      }

      /* Initialize... */
      kpp_chem_initialize(ctl, clim, met0, met1, atm, ip);

      /* Integrate... */
      INTEGRATE(atm->time[ip] - dt[ip], atm->time[ip]);

      /* Output to air parcel.. */
      kpp_chem_output2atm(atm, ctl, ip);

      /* Free... */
      free(VAR);
      free(FIX);
    }
  }
}

#endif

/*****************************************************************************/

void module_position(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt) {

  /* Set timer... */
  SELECT_TIMER("MODULE_POSITION", "PHYSICS", NVTX_GPU);

  const int np = atm->np;
#ifdef _OPENACC
#pragma acc data present(met0, met1, atm, dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < np; ip++)
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
      if (atm->p[ip] < met0->p[met0->np - 1]) {
	if (ctl->reflect)
	  atm->p[ip] = 2. * met0->p[met0->np - 1] - atm->p[ip];
	else
	  atm->p[ip] = met0->p[met0->np - 1];
      } else if (atm->p[ip] > 300.) {
	INTPOL_2D(ps, 1);
	if (atm->p[ip] > ps) {
	  if (ctl->reflect)
	    atm->p[ip] = 2. * ps - atm->p[ip];
	  else
	    atm->p[ip] = ps;
	}
      }
    }
}

/*****************************************************************************/

void module_rng_init(
  int ntask) {

  /* Initialize random number generator... */
#ifdef _OPENACC

  if (curandCreateGenerator(&rng, CURAND_RNG_PSEUDO_DEFAULT) !=
      CURAND_STATUS_SUCCESS)
    ERRMSG("Cannot create random number generator!");
  if (curandSetPseudoRandomGeneratorSeed(rng, ntask) != CURAND_STATUS_SUCCESS)
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
    gsl_rng_set(rng[i],
		gsl_rng_default_seed + (long unsigned) (ntask * NTHREADS +
							i));
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
      if (curandGenerateUniformDouble(rng, rs, (n < 4 ? 4 : n)) !=
	  CURAND_STATUS_SUCCESS)
	ERRMSG("Cannot create random numbers!");
    }

    /* Normal distribution... */
    else if (method == 1) {
      if (curandGenerateNormalDouble(rng, rs, (n < 4 ? 4 : n), 0.0, 1.0) !=
	  CURAND_STATUS_SUCCESS)
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
  SELECT_TIMER("MODULE_SEDI", "PHYSICS", NVTX_GPU);

  const int np = atm->np;
#ifdef _OPENACC
#pragma acc data present(ctl, met0, met1, atm, dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < np; ip++)
    if (dt[ip] != 0) {

      /* Get temperature... */
      double t;
      INTPOL_INIT;
      INTPOL_3D(t, 1);

      /* Sedimentation velocity... */
      double v_s = sedi(atm->p[ip], t, atm->q[ctl->qnt_rp][ip],
			atm->q[ctl->qnt_rhop][ip]);

      /* Calculate pressure change... */
      atm->p[ip] += DZ2DP(v_s * dt[ip] / 1000., atm->p[ip]);
    }
}

/*****************************************************************************/

void module_sort(
  ctl_t * ctl,
  met_t * met0,
  atm_t * atm) {

  /* Set timer... */
  SELECT_TIMER("MODULE_SORT", "MEMORY", NVTX_GPU);

  /* Allocate... */
  const int np = atm->np;
  double *restrict const a = (double *) malloc((size_t) np * sizeof(double));
  int *restrict const p = (int *) malloc((size_t) np * sizeof(int));

#ifdef _OPENACC
#pragma acc enter data create(a[0:np],p[0:np])
#pragma acc data present(ctl,met0,atm,a,p)
#endif

  /* Get box index... */
#ifdef _OPENACC
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < np; ip++) {
    a[ip] =
      (double) ((locate_reg(met0->lon, met0->nx, atm->lon[ip]) * met0->ny +
		 locate_reg(met0->lat, met0->ny, atm->lat[ip]))
		* met0->np + locate_irr(met0->p, met0->np, atm->p[ip]));
    p[ip] = ip;
  }

  /* Sorting... */
#ifdef _OPENACC
  {
#ifdef THRUST
    {
#pragma acc host_data use_device(a, p)
      thrustSortWrapper(a, np, p);
    }
#else
    {
#pragma acc update host(a[0:np], p[0:np])
#pragma omp parallel
      {
#pragma omp single nowait
	quicksort(a, p, 0, np - 1);
      }
#pragma acc update device(a[0:np], p[0:np])
    }
#endif
  }
#else
  {
#ifdef THRUST
    {
      thrustSortWrapper(a, np, p);
    }
#else
    {
#pragma omp parallel
      {
#pragma omp single nowait
	quicksort(a, p, 0, np - 1);
      }
    }
#endif
  }
#endif

  /* Sort data... */
  module_sort_help(atm->time, p, np);
  module_sort_help(atm->p, p, np);
  module_sort_help(atm->lon, p, np);
  module_sort_help(atm->lat, p, np);
  for (int iq = 0; iq < ctl->nq; iq++)
    module_sort_help(atm->q[iq], p, np);

  /* Free... */
#ifdef _OPENACC
#pragma acc exit data delete(a,p)
#endif
  free(a);
  free(p);
}

/*****************************************************************************/

void module_sort_help(
  double *a,
  int *p,
  int np) {

  /* Allocate... */
  double *restrict const help =
    (double *) malloc((size_t) np * sizeof(double));

  /* Reordering of array... */
#ifdef _OPENACC
#pragma acc enter data create(help[0:np])
#pragma acc data present(a,p,help)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < np; ip++)
    help[ip] = a[p[ip]];
#ifdef _OPENACC
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < np; ip++)
    a[ip] = help[ip];

  /* Free... */
#ifdef _OPENACC
#pragma acc exit data delete(help)
#endif
  free(help);
}

/*****************************************************************************/

void module_timesteps(
  ctl_t * ctl,
  atm_t * atm,
  met_t * met0,
  double *dt,
  double t) {

  /* Set timer... */
  SELECT_TIMER("MODULE_TIMESTEPS", "PHYSICS", NVTX_GPU);

  const double latmin = gsl_stats_min(met0->lat, 1, (size_t) met0->ny),
    latmax = gsl_stats_max(met0->lat, 1, (size_t) met0->ny);

  const int np = atm->np,
    local = (fabs(met0->lon[met0->nx - 1] - met0->lon[0] - 360.0) >= 0.01);

#ifdef _OPENACC
#pragma acc data present(ctl, atm, dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < np; ip++) {

    /* Set time step for each air parcel... */
    if ((ctl->direction * (atm->time[ip] - ctl->t_start) >= 0
	 && ctl->direction * (atm->time[ip] - ctl->t_stop) <= 0
	 && ctl->direction * (atm->time[ip] - t) < 0))
      dt[ip] = t - atm->time[ip];
    else
      dt[ip] = 0.0;

    /* Check horizontal boundaries of local meteo data... */
    if (local && (atm->lon[ip] <= met0->lon[0]
		  || atm->lon[ip] >= met0->lon[met0->nx - 1]
		  || atm->lat[ip] <= latmin || atm->lat[ip] >= latmax))
      dt[ip] = 0.0;
  }
}

/*****************************************************************************/

void module_timesteps_init(
  ctl_t * ctl,
  atm_t * atm) {

  /* Set timer... */
  SELECT_TIMER("MODULE_TIMESTEPS", "PHYSICS", NVTX_GPU);

  /* Set start time... */
  if (ctl->direction == 1) {
    ctl->t_start = gsl_stats_min(atm->time, 1, (size_t) atm->np);
    if (ctl->t_stop > 1e99)
      ctl->t_stop = gsl_stats_max(atm->time, 1, (size_t) atm->np);
  } else {
    ctl->t_start = gsl_stats_max(atm->time, 1, (size_t) atm->np);
    if (ctl->t_stop > 1e99)
      ctl->t_stop = gsl_stats_min(atm->time, 1, (size_t) atm->np);
  }

  /* Check time interval... */
  if (ctl->direction * (ctl->t_stop - ctl->t_start) <= 0)
    ERRMSG("Nothing to do! Check T_STOP and DIRECTION!");

  /* Round start time... */
  if (ctl->direction == 1)
    ctl->t_start = floor(ctl->t_start / ctl->dt_mod) * ctl->dt_mod;
  else
    ctl->t_start = ceil(ctl->t_start / ctl->dt_mod) * ctl->dt_mod;
}

/*****************************************************************************/

void module_wet_deposition(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt) {

  /* Set timer... */
  SELECT_TIMER("MODULE_WETDEPO", "PHYSICS", NVTX_GPU);

  /* Check quantity flags... */
  if (ctl->qnt_m < 0 && ctl->qnt_vmr < 0)
    ERRMSG("Module needs quantity mass or volume mixing ratio!");

  const int np = atm->np;
#ifdef _OPENACC
#pragma acc data present(ctl, met0, met1, atm, dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < np; ip++)
    if (dt[ip] != 0) {

      double cl, dz, h, lambda = 0, t, iwc, lwc, pct, pcb;

      /* Check whether particle is below cloud top... */
      INTPOL_INIT;
      INTPOL_2D(pct, 1);
      if (!isfinite(pct) || atm->p[ip] <= pct)
	continue;

      /* Get cloud bottom pressure... */
      INTPOL_2D(pcb, 0);

      /* Estimate precipitation rate (Pisso et al., 2019)... */
      INTPOL_2D(cl, 0);
      double Is =
	pow(1. / ctl->wet_depo_pre[0] * cl, 1. / ctl->wet_depo_pre[1]);
      if (Is < 0.01)
	continue;

      /* Check whether particle is inside or below cloud... */
      INTPOL_3D(lwc, 1);
      INTPOL_3D(iwc, 0);
      int inside = (iwc > 0 || lwc > 0);

      /* Get temperature... */
      INTPOL_3D(t, 0);

      /* Calculate in-cloud scavenging coefficient... */
      if (inside) {

	/* Calculate retention factor... */
	double eta;
	if (t > 273.15)
	  eta = 1;
	else if (t <= 238.15)
	  eta = ctl->wet_depo_ic_ret_ratio;
	else
	  eta = LIN(273.15, 1, 238.15, ctl->wet_depo_ic_ret_ratio, t);

	/* Use exponential dependency for particles ... */
	if (ctl->wet_depo_ic_a > 0)
	  lambda = ctl->wet_depo_ic_a * pow(Is, ctl->wet_depo_ic_b) * eta;

	/* Use Henry's law for gases... */
	else if (ctl->wet_depo_ic_h[0] > 0) {

	  /* Get Henry's constant (Sander, 2015)... */
	  h = ctl->wet_depo_ic_h[0]
	    * exp(ctl->wet_depo_ic_h[1] * (1. / t - 1. / 298.15));

	  /* Use effective Henry's constant for SO2
	     (Berglen, 2004; Simpson, 2012)... */
	  if (ctl->wet_depo_ic_h[2] > 0) {
	    double H_ion = pow(10, ctl->wet_depo_ic_h[2] * (-1));
	    double K_1 = 1.23e-2 * exp(2.01e3 * (1. / t - 1. / 298.15));
	    double K_2 = 6e-8 * exp(1.12e3 * (1. / t - 1. / 298.15));
	    h *= (1 + K_1 / H_ion + K_1 * K_2 / pow(H_ion, 2));
	  }

	  /* Estimate depth of cloud layer... */
	  dz = 1e3 * (Z(pct) - Z(pcb));

	  /* Calculate scavenging coefficient (Draxler and Hess, 1997)... */
	  lambda = h * RI * t * Is / 3.6e6 / dz * eta;
	}
      }

      /* Calculate below-cloud scavenging coefficient... */
      else {

	/* Calculate retention factor... */
	double eta;
	if (t > 270)
	  eta = 1;
	else
	  eta = ctl->wet_depo_bc_ret_ratio;

	/* Use exponential dependency for particles... */
	if (ctl->wet_depo_bc_a > 0)
	  lambda = ctl->wet_depo_bc_a * pow(Is, ctl->wet_depo_bc_b) * eta;

	/* Use Henry's law for gases... */
	else if (ctl->wet_depo_bc_h[0] > 0) {

	  /* Get Henry's constant (Sander, 2015)... */
	  h = ctl->wet_depo_bc_h[0]
	    * exp(ctl->wet_depo_bc_h[1] * (1. / t - 1. / 298.15));

	  /* Estimate depth of cloud layer... */
	  dz = 1e3 * (Z(pct) - Z(pcb));

	  /* Calculate scavenging coefficient (Draxler and Hess, 1997)... */
	  lambda = h * RI * t * Is / 3.6e6 / dz * eta;
	}
      }

      /* Calculate exponential decay of mass... */
      double aux = exp(-dt[ip] * lambda);
      if (ctl->qnt_m >= 0) {
	if (ctl->qnt_mloss_wet >= 0)
	  atm->q[ctl->qnt_mloss_wet][ip]
	    += atm->q[ctl->qnt_m][ip] * (1 - aux);
	atm->q[ctl->qnt_m][ip] *= aux;
      }
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

  char ext[10], filename[2 * LEN];

  double r;

  int year, mon, day, hour, min, sec;

  /* Get time... */
  jsec2time(t, &year, &mon, &day, &hour, &min, &sec, &r);

  /* Update host... */
#ifdef _OPENACC
  if ((ctl->atm_basename[0] != '-' && fmod(t, ctl->atm_dt_out) == 0)
      || (ctl->grid_basename[0] != '-' && fmod(t, ctl->grid_dt_out) == 0)
      || (ctl->ens_basename[0] != '-' && fmod(t, ctl->ens_dt_out) == 0)
      || ctl->csi_basename[0] != '-' || ctl->prof_basename[0] != '-'
      || ctl->sample_basename[0] != '-' || ctl->stat_basename[0] != '-'
      || (ctl->vtk_basename[0] != '-' && fmod(t, ctl->vtk_dt_out) == 0)) {
    SELECT_TIMER("UPDATE_HOST", "MEMORY", NVTX_D2H);
#pragma acc update host(atm[:1])
  }
#endif

  /* Write atmospheric data... */
  if (ctl->atm_basename[0] != '-' && fmod(t, ctl->atm_dt_out) == 0) {
    if (ctl->atm_type == 0)
      sprintf(ext, "tab");
    else if (ctl->atm_type == 1)
      sprintf(ext, "bin");
    else if (ctl->atm_type == 2)
      sprintf(ext, "nc");
    sprintf(filename, "%s/%s_%04d_%02d_%02d_%02d_%02d.%s",
	    dirname, ctl->atm_basename, year, mon, day, hour, min, ext);
    write_atm(filename, ctl, atm, t);
  }

  /* Write gridded data... */
  if (ctl->grid_basename[0] != '-' && fmod(t, ctl->grid_dt_out) == 0) {
    sprintf(filename, "%s/%s_%04d_%02d_%02d_%02d_%02d.%s",
	    dirname, ctl->grid_basename, year, mon, day, hour, min,
	    ctl->grid_type == 0 ? "tab" : "nc");
    write_grid(filename, ctl, met0, met1, atm, t);
  }

  /* Write CSI data... */
  if (ctl->csi_basename[0] != '-') {
    sprintf(filename, "%s/%s.tab", dirname, ctl->csi_basename);
    write_csi(filename, ctl, atm, t);
  }

  /* Write ensemble data... */
  if (ctl->ens_basename[0] != '-' && fmod(t, ctl->ens_dt_out) == 0) {
    sprintf(filename, "%s/%s_%04d_%02d_%02d_%02d_%02d.tab",
	    dirname, ctl->ens_basename, year, mon, day, hour, min);
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

  /* Write VTK data... */
  if (ctl->vtk_basename[0] != '-' && fmod(t, ctl->vtk_dt_out) == 0) {
    static int nvtk;
    sprintf(filename, "%s/%s_%05d.vtk", dirname, ctl->vtk_basename, ++nvtk);
    write_vtk(filename, ctl, atm, t);
  }
}
