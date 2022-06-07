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
  
  Copyright (C) 2013-2021 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Calculate PSC temperatures.
*/

#include "libtrac.h"

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <p> <h2o> <hno3>");

  /* Get varibles... */
  double p = atof(argv[1]);
  double h2o = atof(argv[2]);
  double hno3 = atof(argv[3]);

  /* Calculate T_ice and T_NAT... */
  double tice = TICE(p, h2o);
  double tnat = nat_temperature(p, h2o, hno3);

  /* Write output... */
  printf("p= %g hPa\n", p);
  printf("q_H2O= %g ppv\n", h2o);
  printf("q_HNO3= %g ppv\n", hno3);
  printf("T_ice= %g K\n", tice);
  printf("T_NAT= %g K\n", tnat);

  return EXIT_SUCCESS;
}
