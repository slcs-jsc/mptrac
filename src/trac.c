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

  Copyright (C) 2013-2026 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Lagrangian particle dispersion model.
*/

#include "mptrac.h"

#ifdef KPP
#include "kpp_chem.h"
#endif

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Print command-line help. */
void usage(
  void);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  ctl_t *ctl;

  atm_t *atm;

  cache_t *cache;

  clim_t *clim;

  met_t *met0, *met1;

  dd_t *dd;

  FILE *dirlist;

  char dirname[LEN], filename[2 * LEN];

  int ntask = -1, rank = 0, size = 1;

  /* Print usage information... */
  USAGE;

  /* Initialize MPI... */
#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: trac <dirlist> <ctl> <atm_in> [KEY VALUE ...]\n\n"
	   "Use -h for full help.");

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

    /* Allocate memory... */
    mptrac_alloc(&ctl, &cache, &clim, &met0, &met1, &atm, &dd);

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
      mptrac_get_met(ctl, clim, t, &met0, &met1, dd);

      /* Check time step... */
      if (ctl->dt_mod > fabs(met0->lon[1] - met0->lon[0]) * 111132. / 150.)
	WARN("Violation of CFL criterion! Check DT_MOD!");

#ifdef DD
      /* Set-up domain decomposition... */
      if ((t == ctl->t_start) && (ctl->dd == 1))
	dd_init(ctl, dd, atm);
#endif

      /* Run a single time step... */
      mptrac_run_timestep(ctl, cache, clim, &met0, &met1, atm, t, dd);

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
    mptrac_free(ctl, cache, clim, met0, met1, atm, dd);

    /* Report timers... */
    PRINT_TIMERS;
  }

  /* Finalize MPI... */
#ifdef MPI
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}

/*****************************************************************************/

/*! Print command-line help. */
void usage(
  void) {

  printf("\nMPTRAC trac tool.\n\n");
  printf("Run forward or backward trajectory calculations.\n");
  printf("\n");
  printf("Usage:\n");
  printf("  trac <dirlist> <ctl> <atm_in> [KEY VALUE ...]\n");
  printf("\n");
  printf("Arguments:\n");
  printf("  <dirlist>  Text file containing work directories to process.\n");
  printf("  <ctl>      Control file name relative to each work directory.\n");
  printf
    ("  <atm_in>   Atmospheric input file name relative to each work directory.\n");
  printf("  [KEY VALUE]  Optional control parameters.\n");
  printf("\nCommon control parameters:\n");
  printf("  T_STOP, DT_MOD                   Simulation end time and model time step.\n");
  printf("  DIRECTION                        Time-integration direction.\n");
  printf("  METBASE, DT_MET                  Meteorological input basename and step.\n");
  printf("  ATM_BASENAME, ATM_DT_OUT         Atmospheric output basename and step.\n");
  printf("  GRID_BASENAME, GRID_DT_OUT       Gridded output basename and step.\n");
  printf("  ENS_BASENAME, STAT_BASENAME      Ensemble and statistics output basenames.\n");
  printf("  SAMPLE_BASENAME, PROF_BASENAME   Sample and profile output basenames.\n");
  printf("\nFurther information:\n");
  printf("  Manual: https://slcs-jsc.github.io/mptrac/\n");
}
