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
  Lagrangian particle dispersion model.
*/

#include "mptrac.h"

#ifdef KPP
#include "kpp_chem.h"
#endif

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

  double *dt, *rs;

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

    /* Initialize pressure heights consistent with zeta... */
    if (ctl.vert_coord_ap == 1) {
#pragma omp parallel for default(shared)
      for (int ip = 0; ip < atm->np; ip++) {
	INTPOL_INIT;
	intpol_met_4d_coord(met0, met0->zetal, met0->pl, met1, met1->zetal,
			    met1->pl, atm->time[ip], atm->q[ctl.qnt_zeta][ip],
			    atm->lon[ip], atm->lat[ip], &atm->p[ip], ci, cw,
			    1);
      }
    }

    /* Initialize species quantity values according to meteorological data or climatology... */
    module_quan_init(&ctl, clim, met0, met1, atm);

    /* Update GPU... */
#ifdef _OPENACC
    SELECT_TIMER("UPDATE_DEVICE", "MEMORY", NVTX_H2D);
#pragma acc update device(atm[:1], cache[:1], clim[:1], ctl)
#endif

    /* ------------------------------------------------------------
       Loop over timesteps...
       ------------------------------------------------------------ */

    /* Loop over timesteps... */
#ifdef ASYNCIO
    omp_set_nested(1);
    int ompTrdnum = omp_get_max_threads();
#endif
    for (double t = ctl.t_start;
	 ctl.direction * (t - ctl.t_stop) < ctl.dt_mod;
	 t += ctl.direction * ctl.dt_mod) {
#ifdef ASYNCIO
#pragma omp parallel num_threads(2)
      {
#endif

	/* Adjust length of final time step... */
	if (ctl.direction * (t - ctl.t_stop) > 0)
	  t = ctl.t_stop;

	/* Set time steps of air parcels... */
	module_timesteps(&ctl, met0, atm, dt, t);

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
	  if (ctl.advect > 0) {
	    if (ctl.vert_coord_ap == 0)
	      module_advect(&ctl, met0, met1, atm, dt);
	    else
	      module_advect_diabatic(&ctl, met0, met1, atm, dt);
	  }

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
	    module_isosurf(&ctl, met0, met1, atm, cache, dt);

	  /* Check positions (final)... */
	  module_position(&ctl, met0, met1, atm, dt);

	  /* Interpolate meteo data... */
	  if (ctl.met_dt_out > 0
	      && (ctl.met_dt_out < ctl.dt_mod
		  || fmod(t, ctl.met_dt_out) == 0))
	    module_meteo(&ctl, clim, met0, met1, atm, dt);

	  /* Check boundary conditions (initial)... */
	  if ((ctl.bound_lat0 < ctl.bound_lat1)
	      && (ctl.bound_p0 > ctl.bound_p1))
	    module_bound_cond(&ctl, clim, met0, met1, atm, dt);

	  /* Initialize quantity of total loss rate... */
	  if (ctl.qnt_loss_rate >= 0) {
	    PARTICLE_LOOP(0, atm->np, 1, "acc data present(ctl, atm)") {
	      atm->q[ctl.qnt_loss_rate][ip] = 0;
	    }
	  }

	  /* Decay of particle mass... */
	  if (ctl.tdec_trop > 0 && ctl.tdec_strat > 0)
	    module_decay(&ctl, clim, atm, dt);

	  /* Interparcel mixing... */
	  if (ctl.mixing_trop >= 0 && ctl.mixing_strat >= 0
	      && (ctl.mixing_dt <= 0 || fmod(t, ctl.mixing_dt) == 0))
	    module_mixing(&ctl, clim, atm, t);

	  /* Calculate the tracer vmr in the chemistry grid... */
	  if (ctl.qnt_Cx > 0
	      && (ctl.oh_chem_reaction != 0 || ctl.h2o2_chem_reaction != 0
		  || (ctl.kpp_chem && fmod(t, ctl.dt_kpp) == 0)))
	    module_chemgrid(&ctl, met0, met1, atm, t);

	  /* OH chemistry... */
	  if (ctl.oh_chem_reaction != 0)
	    module_oh_chem(&ctl, clim, met0, met1, atm, dt);

	  /* H2O2 chemistry (for SO2 aqueous phase oxidation)... */
	  if (ctl.h2o2_chem_reaction != 0)
	    module_h2o2_chem(&ctl, clim, met0, met1, atm, dt);

	  /* First-order tracer chemistry... */
	  if (ctl.tracer_chem)
	    module_tracer_chem(&ctl, clim, met0, met1, atm, dt);

	  /* KPP chemistry... */
	  if (ctl.kpp_chem && fmod(t, ctl.dt_kpp) == 0) {
#ifdef KPP
	    module_kpp_chem(&ctl, clim, met0, met1, atm, dt);
#else
	    ERRMSG("Code was compiled without KPP!");
#endif
	  }

	  /* Wet deposition... */
	  if ((ctl.wet_depo_ic_a > 0 || ctl.wet_depo_ic_h[0] > 0)
	      && (ctl.wet_depo_bc_a > 0 || ctl.wet_depo_bc_h[0] > 0))
	    module_wet_deposition(&ctl, met0, met1, atm, dt);

	  /* Dry deposition... */
	  if (ctl.dry_depo_vdep > 0)
	    module_dry_deposition(&ctl, met0, met1, atm, dt);

	  /* Check boundary conditions (final)... */
	  if ((ctl.bound_lat0 < ctl.bound_lat1)
	      && (ctl.bound_p0 > ctl.bound_p1))
	    module_bound_cond(&ctl, clim, met0, met1, atm, dt);

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

    /* Flush output buffer... */
    fflush(NULL);

    /* Report problem size... */
    LOG(1, "SIZE_NP = %d", atm->np);
    LOG(1, "SIZE_MPI_TASKS = %d", size);
    LOG(1, "SIZE_OMP_THREADS = %d", omp_get_max_threads());
#ifdef _OPENACC
    LOG(1, "SIZE_ACC_DEVICES = %d", acc_get_num_devices(acc_device_nvidia));
#else
    LOG(1, "SIZE_ACC_DEVICES = %d", 0);
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
