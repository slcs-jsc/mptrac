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
  Convert file format of meteo data files.
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

  clim_t *clim;

  met_t *met;

  dd_t *dd;

  /* Print usage information... */
  USAGE;

  /* Check arguments... */
  if (argc < 6)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: met_conv <ctl> <met_in> <met_in_type> <met_out> <met_out_type> [KEY VALUE ...]\n\n"
	   "Use -h for full help.");

  /* Allocate... */
  ALLOC(clim, clim_t, 1);
  ALLOC(met, met_t, 1);
  ALLOC(dd, dd_t, 1);

  /* Read control parameters... */
  mptrac_read_ctl(argv[1], argc, argv, &ctl);

  /* Read climatological data... */
  mptrac_read_clim(&ctl, clim);

  /* Read meteo data... */
  ctl.met_type = atoi(argv[3]);
  if (!mptrac_read_met(argv[2], &ctl, clim, met, dd))
    ERRMSG("Cannot open file!");

  /* Write meteo data... */
  ctl.met_type = atoi(argv[5]);
  mptrac_write_met(argv[4], &ctl, met);

  /* Report timers... */
  PRINT_TIMERS;

  /* Free... */
  free(clim);
  free(met);
  free(dd);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

/*! Print command-line help. */
void usage(
  void) {

  printf("\nMPTRAC met_conv tool.\n\n");
  printf("Convert meteorological data between file formats.\n");
  printf("\n");
  printf("Usage:\n");
  printf
    ("  met_conv <ctl> <met_in> <met_in_type> <met_out> <met_out_type> [KEY VALUE ...]\n");
  printf("\n");
  printf("Arguments:\n");
  printf("  <ctl>           Control file.\n");
  printf("  <met_in>        Meteorological input file.\n");
  printf
    ("  <met_in_type>   Input format: 0=netCDF, 1=binary, 2=pck, 3=zfp,\n");
  printf("                  4=zstd, 5=cms, 6=GRIB, 7=SZ3.\n");
  printf("  <met_out>       Meteorological output file.\n");
  printf
    ("  <met_out_type>  Output format: 0=netCDF, 1=binary, 2=pck, 3=zfp,\n");
  printf("                  4=zstd, 5=cms, 6=GRIB, 7=SZ3.\n");
  printf("  [KEY VALUE]     Optional control parameters.\n");
  printf("\nFurther information:\n");
  printf("  Manual: https://slcs-jsc.github.io/mptrac/\n");
}
