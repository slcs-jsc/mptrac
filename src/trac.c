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

  /* Initialize MPI... */
#ifdef MPI
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

    /* Write info... */
    LOG(1, "Parallelization: ntask= %d | rank= %d | size= %d",
	ntask, rank, size);

    /* ------------------------------------------------------------
       Initialize model run...
       ------------------------------------------------------------ */

    /* Start timers... */
    START_TIMERS;

    /* Allocate memory... */
    mptrac_alloc(&ctl, &cache, &clim, &met0, &met1, &atm);
    
    /* Register the MPI_Particle data type... */
    MPI_Datatype MPI_Particle;
    dd_reg_MPI_type_particle(&MPI_Particle);
  
    /* Define communication destinations ... */
    int destinations[8];
    dd_get_destinations(ctl, destinations, rank, size); 
  
    /* Check if particles are in domain. */
    dd_assign_rect_domains_atm( atm, met0, ctl, rank, destinations, 1);
   
    /* Initialize particles locally... */
    particle_t* particles;
    ALLOC(particles, particle_t, atm->np);

    /* Read control parameters... */
    sprintf(filename, "%s/%s", dirname, argv[2]);
    mptrac_read_ctl(filename, argc, argv, ctl);

    /* Read climatological data... */
    mptrac_read_clim(ctl, clim);

    /* Read atmospheric data... */
    sprintf(filename, "%s/%s", dirname, argv[3]);
    if (!mptrac_read_atm(filename, ctl, atm))
      ERRMSG("Cannot open file!");

    /* Initialize MPTRAC... */
    mptrac_init(ctl, cache, clim, atm, ntask);

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

      /* Get meteo data... */
      mptrac_get_met(ctl, clim, t, &met0, &met1);

      /* Check time step... */
      if (ctl->dt_mod > fabs(met0->lon[1] - met0->lon[0]) * 111132. / 150.)
	WARN("Violation of CFL criterion! Check DT_MOD!");

      /* Run a single time step... */
      mptrac_run_timestep(ctl, cache, clim, &met0, &met1, atm, t);

      /* Write output... */
      mptrac_write_output(dirname, ctl, met0, met1, atm, t);
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

    /* Report memory usage... */
    LOG(1, "MEMORY_ATM = %g MByte", sizeof(atm_t) / 1024. / 1024.);
    LOG(1, "MEMORY_CACHE = %g MByte", sizeof(cache_t) / 1024. / 1024.);
    LOG(1, "MEMORY_CLIM = %g MByte", sizeof(clim_t) / 1024. / 1024.);
    LOG(1, "MEMORY_METEO = %g MByte", sizeof(met_t) / 1024. / 1024.);

    /* Free memory... */
    mptrac_free(ctl, cache, clim, met0, met1, atm);
    /* Free MPI datatype... */
    MPI_Type_free(&MPI_Particle);
    
    /* Free local particle array... */
    free(particles);

    /* Report timers... */
    PRINT_TIMERS;
    STOP_TIMERS;
  }

  /* Finalize MPI... */
#ifdef MPI
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}
