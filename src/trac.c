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

  Copyright (C) 2013-2025 Forschungszentrum Juelich GmbH
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

  ctl_t *ctl;

  atm_t *atm;

  cache_t *cache;

  clim_t *clim;

  met_t *met0, *met1;

  FILE *dirlist;

  char dirname[LEN], filename[2 * LEN];

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

#ifdef _OPENACC
    LOG(1, "Parallelization: ntask= %d | rank= %d | size= %d | acc_dev= %d",
	ntask, rank, size, acc_get_device_num(device_type));
#else
    LOG(1, "Parallelization: ntask= %d | rank= %d | size= %d | acc_dev= nan",
	ntask, rank, size);
#endif

    /* Allocate... */
    SELECT_TIMER("ALLOC", "MEMORY", NVTX_CPU);
    ALLOC(ctl, ctl_t, 1);
    ALLOC(atm, atm_t, 1);
    ALLOC(cache, cache_t, 1);
    ALLOC(clim, clim_t, 1);
    ALLOC(met0, met_t, 1);
    ALLOC(met1, met_t, 1);

    /* Create data region on GPUs... */
#ifdef _OPENACC
    SELECT_TIMER("CREATE_DATA_REGION", "MEMORY", NVTX_GPU);
#pragma acc enter data create(atm[:1], cache[:1], clim[:1], ctl[:1], met0[:1], met1[:1])
#endif

    /* Read control parameters... */
    sprintf(filename, "%s/%s", dirname, argv[2]);
    read_ctl(filename, argc, argv, ctl);

    /* Read climatological data... */
    read_clim(ctl, clim);

    /* Read atmospheric data... */
    sprintf(filename, "%s/%s", dirname, argv[3]);
    if (!read_atm(filename, ctl, atm))
      ERRMSG("Cannot open file!");

    /* Initialize timesteps... */
    module_timesteps_init(ctl, atm);

    /* Initialize random number generator... */
    module_rng_init(ntask);

    /* Initialize meteo data... */
    get_met(ctl, clim, ctl->t_start, &met0, &met1);

    /* Check time step... */
    if (ctl->dt_mod > fabs(met0->lon[1] - met0->lon[0]) * 111132. / 150.)
      WARN("Violation of CFL criterion! Check DT_MOD!");

    /* Initialize isosurface data... */
    if (ctl->isosurf >= 1 && ctl->isosurf <= 4)
      module_isosurf_init(ctl, cache, met0, met1, atm);

    /* Initialize advection... */
    module_advect_init(ctl, met0, met1, atm);

    /* Initialize chemistry... */
    module_chem_init(ctl, clim, met0, met1, atm);

    /* Update GPU... */
#ifdef _OPENACC
    SELECT_TIMER("UPDATE_DEVICE", "MEMORY", NVTX_H2D);
#pragma acc update device(atm[:1], cache[:1], clim[:1], ctl[:1])
#endif

    /* ------------------------------------------------------------
       Loop over timesteps...
       ------------------------------------------------------------ */

    /* Loop over timesteps... */
    for (double t = ctl->t_start;
	 ctl->direction * (t - ctl->t_stop) < ctl->dt_mod;
	 t += ctl->direction * ctl->dt_mod) {

      /* Adjust length of final time step... */
      if (ctl->direction * (t - ctl->t_stop) > 0)
	t = ctl->t_stop;

      /* Run a single time step... */
      mptrac_run_timestep(ctl, cache, clim, &met0, &met1, atm, t);

      /* Write output... */
      write_output(dirname, ctl, met0, met1, atm, t);
    }

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
#pragma acc exit data delete (ctl, atm, cache, clim, met0, met1)
#endif

    /* Free... */
    SELECT_TIMER("FREE", "MEMORY", NVTX_CPU);
    free(atm);
    free(ctl);
    free(cache);
    free(clim);
    free(met0);
    free(met1);
  }

  /* Report timers... */
  PRINT_TIMERS;

  /* Finalize MPI... */
#ifdef MPI
  MPI_Finalize();
#endif

  /* Stop timers... */
  STOP_TIMERS;

  return EXIT_SUCCESS;
}
