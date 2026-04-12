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
  Convert Julian seconds to date.
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

  double remain;

  int day, hour, min, mon, sec, year;

  /* Print usage information... */
  USAGE;

  /* Check arguments... */
  if (argc < 2)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: jsec2time <jsec>\n\n" "Use -h for full help.");

  /* Read arguments... */
  const double jsec = atof(argv[1]);

  /* Convert time... */
  jsec2time(jsec, &year, &mon, &day, &hour, &min, &sec, &remain);
  printf("%d %d %d %d %d %d %g\n", year, mon, day, hour, min, sec, remain);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

/*! Print command-line help. */
void usage(
  void) {

  printf("\nMPTRAC jsec2time tool.\n\n");
  printf("Convert seconds since 2000-01-01 00:00 UTC to calendar time.\n");
  printf("\n");
  printf("Usage:\n");
  printf("  jsec2time <jsec>\n");
  printf("\n");
  printf("Arguments:\n");
  printf("  <jsec>  Seconds since 2000-01-01 00:00 UTC.\n");
  printf("\nFurther information:\n");
  printf("  Manual: https://slcs-jsc.github.io/mptrac/\n");
}
