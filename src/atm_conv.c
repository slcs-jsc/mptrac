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
  Convert file format of air parcel data files.
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

  /* Print usage information... */
  USAGE;

  /* Check arguments... */
  if (argc < 6)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: atm_conv <ctl> <atm_in> <atm_in_type> <atm_out> <atm_out_type>\n\n"
	   "Use -h for full help.");

  /* Allocate... */
  ALLOC(atm, atm_t, 1);

  /* Read control parameters... */
  mptrac_read_ctl(argv[1], argc, argv, &ctl);

  /* Read atmospheric data... */
  ctl.atm_type = atoi(argv[3]);
  if (!mptrac_read_atm(argv[2], &ctl, atm))
    ERRMSG("Cannot open file!");

  /* Write atmospheric data... */
  if (ctl.atm_type_out == 3) {

    /* For CLaMS trajectory files... */
    ctl.t_start = ctl.t_stop;
    ctl.atm_type_out = atoi(argv[5]);
    mptrac_write_atm(argv[4], &ctl, atm, ctl.t_stop);

  } else {

    /* Otherwise... */
    ctl.atm_type_out = atoi(argv[5]);
    mptrac_write_atm(argv[4], &ctl, atm, 0);
  }

  /* Free... */
  free(atm);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

/*! Print command-line help. */
void usage(
  void) {

  printf("\nMPTRAC atm_conv tool.\n\n");
  printf("Convert atmospheric particle data between file formats.\n");
  printf("\n");
  printf("Usage:\n");
  printf
    ("  atm_conv <ctl> <atm_in> <atm_in_type> <atm_out> <atm_out_type>\n");
  printf("\n");
  printf("Arguments:\n");
  printf("  <ctl>           Control file.\n");
  printf("  <atm_in>        Atmospheric input file.\n");
  printf("  <atm_in_type>   Input format: 0=ASCII, 1=binary, 2=netCDF,\n");
  printf
    ("                  3=CLaMS trajectory/position netCDF, 4=CLaMS position netCDF.\n");
  printf("  <atm_out>       Atmospheric output file.\n");
  printf("  <atm_out_type>  Output format: 0=ASCII, 1=binary, 2=netCDF,\n");
  printf
    ("                  3=CLaMS trajectory/position netCDF, 4=CLaMS position netCDF.\n");
  printf("\nFurther information:\n");
  printf("  Manual: https://slcs-jsc.github.io/mptrac/\n");
}
