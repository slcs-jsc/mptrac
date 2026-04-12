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
  Convert atmospheric data file to grid data file.
*/

#include "mptrac.h"

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

  ctl_t ctl;

  atm_t *atm;

  met_t *met0 = NULL, *met1 = NULL;

  /* Allocate... */
  ALLOC(atm, atm_t, 1);

  /* Print usage information... */
  USAGE;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: atm2grid <ctl> <atm_in> [KEY VALUE ...]\n\n"
	   "Use -h for full help.");

  /* Read control parameters... */
  mptrac_read_ctl(argv[1], argc, argv, &ctl);

  /* Check grid basename... */
  if (ctl.grid_basename[0] == '-')
    ERRMSG("You need to specify GRID_BASENAME!");

  /* Read atmospheric data... */
  if (!mptrac_read_atm(argv[2], &ctl, atm))
    ERRMSG("Cannot open file!");

  /* Get time from filename... */
  int year, mon, day, hour, min, sec;
  double r, t = time_from_filename(argv[2], ctl.atm_type < 2 ? 20 : 19);
  jsec2time(t, &year, &mon, &day, &hour, &min, &sec, &r);

  /* Set output filename... */
  char filename[3 * LEN];
  sprintf(filename, "%s_%04d_%02d_%02d_%02d_%02d.%s",
	  ctl.grid_basename, year, mon, day, hour, min,
	  ctl.grid_type == 0 ? "tab" : "nc");

  /* Write grid data... */
  write_grid(filename, &ctl, met0, met1, atm, t);

  /* Free... */
  free(atm);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

/*! Print command-line help. */
void usage(
  void) {

  printf("\nMPTRAC atm2grid tool.\n\n");
  printf
    ("Convert an atmospheric particle data file to gridded output data.\n");
  printf("\n");
  printf("Usage:\n");
  printf("  atm2grid <ctl> <atm_in> [KEY VALUE ...]\n");
  printf("\n");
  printf("Arguments:\n");
  printf("  <ctl>     Control file.\n");
  printf("  <atm_in>  Atmospheric input file.\n");
  printf("  [KEY VALUE]  Optional control parameters.\n");
  printf("\n");
  printf("Notes:\n");
  printf("  GRID_BASENAME must be set in the control parameters.\n");
  printf("\nFurther information:\n");
  printf("  Manual: https://slcs-jsc.github.io/mptrac/\n");
}
