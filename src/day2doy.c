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
  Convert date to day of year.
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

  /* Print usage information... */
  USAGE;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: day2doy <year> <mon> <day>\n\n" "Use -h for full help.");

  /* Read arguments... */
  const int year = atoi(argv[1]);
  const int mon = atoi(argv[2]);
  const int day = atoi(argv[3]);

  /* Convert... */
  int doy;
  day2doy(year, mon, day, &doy);
  printf("%d %d\n", year, doy);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

/*! Print command-line help. */
void usage(
  void) {

  printf("\nMPTRAC day2doy tool.\n\n");
  printf("Convert a calendar date to day of year.\n");
  printf("\n");
  printf("Usage:\n");
  printf("  day2doy <year> <mon> <day>\n");
  printf("\n");
  printf("Arguments:\n");
  printf("  <year>  Year.\n");
  printf("  <mon>   Month.\n");
  printf("  <day>   Day of month.\n");
  printf("\nFurther information:\n");
  printf("  Manual: https://slcs-jsc.github.io/mptrac/\n");
}
