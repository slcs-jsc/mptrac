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
  Calculate sedimentation velocity.
*/

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {

  double p, T, r_p, rho_p;

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <p> <T> <r_p> <rho_p>");

  /* Read arguments... */
  p = atof(argv[1]);
  T = atof(argv[2]);
  r_p = atof(argv[3]);
  rho_p = atof(argv[4]);

  /* Convert... */
  printf("%g\n", sedi(p, T, r_p, rho_p));

  return EXIT_SUCCESS;
}
