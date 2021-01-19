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
  
  Copyright (C) 2013-2019 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Calculate humidity parameters.
*/

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <p> <T> <h2o>");

  /* Read arguments... */
  double p = atof(argv[1]);
  double t = atof(argv[2]);
  double h2o = atof(argv[3]);

  /* Convert time... */
  printf("p= %g T= %g H2O= %g SH= %g PS= %g RH= %g Tdew= %g"
	 " PSice= %g RHice= %g Tice= %g lapse= %g\n",
	 p, t, h2o, SH(h2o), PSAT(t), RH(p, t, h2o), TDEW(p, h2o),
	 PSICE(t), RHICE(p, t, h2o), TICE(p, h2o), lapse_rate(t, h2o));

  return EXIT_SUCCESS;
}
