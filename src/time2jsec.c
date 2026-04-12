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
  Convert date to Julian seconds.
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

  double jsec;

  /* Print usage information... */
  USAGE;

  /* Check arguments... */
  if (argc < 8)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: time2jsec <year> <mon> <day> <hour> <min> <sec> <remain>\n\n"
	   "Use -h for full help.");

  /* Read arguments... */
  const int year = atoi(argv[1]);
  const int mon = atoi(argv[2]);
  const int day = atoi(argv[3]);
  const int hour = atoi(argv[4]);
  const int min = atoi(argv[5]);
  const int sec = atoi(argv[6]);
  const double remain = atof(argv[7]);

  /* Convert... */
  time2jsec(year, mon, day, hour, min, sec, remain, &jsec);
  printf("%.2f\n", jsec);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

/*! Print command-line help. */
void usage(
  void) {

  printf("\nMPTRAC time2jsec tool.\n\n");
  printf("Convert calendar time to seconds since 2000-01-01 00:00 UTC.\n");
  printf("\n");
  printf("Usage:\n");
  printf("  time2jsec <year> <mon> <day> <hour> <min> <sec> <remain>\n");
  printf("\n");
  printf("Arguments:\n");
  printf("  <year>    Year.\n");
  printf("  <mon>     Month.\n");
  printf("  <day>     Day of month.\n");
  printf("  <hour>    Hour.\n");
  printf("  <min>     Minute.\n");
  printf("  <sec>     Second.\n");
  printf("  <remain>  Fractional-second remainder.\n");
  printf("\nFurther information:\n");
  printf("  Manual: https://slcs-jsc.github.io/mptrac/\n");
}
