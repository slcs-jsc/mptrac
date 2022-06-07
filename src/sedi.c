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

  double eta, p, T, r_p, rho, Re, rho_p, vs;

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <p> <T> <r_p> <rho_p>");

  /* Read arguments... */
  p = atof(argv[1]);
  T = atof(argv[2]);
  r_p = atof(argv[3]);
  rho_p = atof(argv[4]);

  /* Calculate sedimentation velocity... */
  vs = sedi(p, T, r_p, rho_p);

  /* Density of dry air [kg / m^3]... */
  rho = 100. * p / (RA * T);

  /* Dynamic viscosity of air [kg / (m s)]... */
  eta = 1.8325e-5 * (416.16 / (T + 120.)) * pow(T / 296.16, 1.5);

  /* Particle Reynolds number... */
  Re = 2e-6 * r_p * vs * rho / eta;

  /* Convert... */
  printf("%g %g\n", vs, Re);

  return EXIT_SUCCESS;
}
