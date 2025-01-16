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
  
  Copyright (C) 2013-2025 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Calculate sedimentation velocity.
*/

#include "mptrac.h"

int main(
  int argc,
  char *argv[]) {

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <p> <T> <r_p> <rho_p>");

  /* Read arguments... */
  const double p = atof(argv[1]);
  const double T = atof(argv[2]);
  const double r_p = atof(argv[3]);
  const double rho_p = atof(argv[4]);

  /* Calculate sedimentation velocity... */
  const double vs = sedi(p, T, r_p, rho_p);

  /* Density of dry air [kg / m^3]... */
  const double rho = 100. * p / (RA * T);

  /* Dynamic viscosity of air [kg / (m s)]... */
  const double eta = 1.8325e-5 * (416.16 / (T + 120.)) * pow(T / 296.16, 1.5);

  /* Particle Reynolds number... */
  const double Re = 2e-6 * r_p * vs * rho / eta;

  /* Write output... */
  printf("    p= %g hPa\n", p);
  printf("    T= %g K\n", T);
  printf("  r_p= %g microns\n", r_p);
  printf("rho_p= %g kg/m^3\n", rho_p);
  printf("rho_a= %g kg/m^3\n", RHO(p, T));
  printf("  v_s= %g m/s\n", vs);
  printf("   Re= %g\n", Re);

  return EXIT_SUCCESS;
}
