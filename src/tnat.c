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
  Calculate PSC temperatures.
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
  if (argc < 3)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: tnat <p> <h2o> <hno3>\n\n" "Use -h for full help.");

  /* Get varibles... */
  const double p = atof(argv[1]);
  const double h2o = atof(argv[2]);
  const double hno3 = atof(argv[3]);

  /* Write output... */
  printf("     p= %g hPa\n", p);
  printf(" q_H2O= %g ppv\n", h2o);
  printf("q_HNO3= %g ppv\n", hno3);
  printf(" T_dew= %g K\n", TDEW(p, h2o));
  printf(" T_ice= %g K\n", TICE(p, h2o));
  printf(" T_NAT= %g K\n", nat_temperature(p, h2o, hno3));

  return EXIT_SUCCESS;
}

/*****************************************************************************/

/*! Print command-line help. */
void usage(
  void) {

  printf("\nMPTRAC tnat tool.\n\n");
  printf("Calculate PSC formation temperatures.\n");
  printf("\n");
  printf("Usage:\n");
  printf("  tnat <p> <h2o> <hno3>\n");
  printf("\n");
  printf("Arguments:\n");
  printf("  <p>     Pressure [hPa].\n");
  printf("  <h2o>   Water vapor volume mixing ratio [ppv].\n");
  printf("  <hno3>  Nitric acid volume mixing ratio [ppv].\n");
  printf("\nFurther information:\n");
  printf("  Manual: https://slcs-jsc.github.io/mptrac/\n");
}
