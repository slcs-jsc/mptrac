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
  Check model time step based on given meeorological data.
*/

#include "mptrac.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  clim_t *clim;

  met_t *met;

  dd_t *dd;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <met>");

  /* Allocate... */
  ALLOC(clim, clim_t, 1);
  ALLOC(met, met_t, 1);
  ALLOC(dd, dd_t, 1);

  /* Read control parameters... */
  mptrac_read_ctl(argv[1], argc, argv, &ctl);
  double kx = scan_ctl(argv[1], argc, argv, "KX", -1, "50.0", NULL);
  double kz = scan_ctl(argv[1], argc, argv, "KZ", -1, "0.1", NULL);

  /* Read climatological data... */
  mptrac_read_clim(&ctl, clim);

  /* Read meteo data... */
  if (!mptrac_read_met(argv[2], &ctl, clim, met, dd))
    ERRMSG("Cannot open file!");

  /* Set parameters... */
  const double c_max = 1.0;
  const double n = 0.3;

  /* Loop over pressure levels... */
  for (int ip = 1; ip < met->np - 1; ip++) {

    /* Init... */
    double dt_x_min = 1e100;
    double dt_y_min = 1e100;
    double dt_p_min = 1e100;
    double dt_dx_min = 1e100;
    double dt_dy_min = 1e100;
    double dt_dp_min = 1e100;

    /* Loop over columns... */
#pragma omp parallel for default(shared) collapse(2) reduction(min:dt_x_min,dt_y_min,dt_p_min,dt_dx_min,dt_dy_min,dt_dp_min)
    for (int ix = 0; ix < met->nx; ix++)
      for (int iy = 1; iy < met->ny - 1; iy++) {

	/* Check advection... */
	const double dx =
	  1e3 * fabs(DEG2DX(met->lon[1] - met->lon[0], met->lat[iy]));
	const double dt_x = fabs(c_max * dx / met->u[ix][iy][ip]);
	if (met->u[ix][iy][ip] != 0 && dt_x < dt_x_min)
	  dt_x_min = dt_x;

	const double dy =
	  0.5 * 1e3 * fabs(DEG2DY(met->lat[iy + 1] - met->lat[iy - 1]));
	const double dt_y = fabs(c_max * dy / met->v[ix][iy][ip]);
	if (met->v[ix][iy][ip] != 0 && dt_y < dt_y_min)
	  dt_y_min = dt_y;

	const double dp = 0.5 * fabs(met->p[ip + 1] - met->p[ip - 1]);
	const double dt_p = fabs(c_max * dp / met->w[ix][iy][ip]);
	if (met->w[ix][iy][ip] != 0 && dt_p < dt_p_min)
	  dt_p_min = dt_p;

	/* Check diffusion... */
	const double dt_dx = 0.5 * SQR(n * dx) / kx;
	if (dt_dx < dt_dx)
	  dt_dx_min = dt_dx;

	const double dt_dy = 0.5 * SQR(n * dy) / kx;
	if (dt_dy < dt_dy_min)
	  dt_dy_min = dt_dy;

	const double dt_dp = 0.5 * SQR(n * dp) / (SQR(met->p[ip] / H0 * 1e-3) * kz);
	if (dt_dp < dt_dp_min)
	  dt_dp_min = dt_dp;
      }

    /* Write output... */
    printf("%g %g %g %g %g %g %g\n", met->p[ip], dt_x_min, dt_y_min, dt_p_min, dt_dx_min, dt_dy_min, dt_dp_min);
  }

  /* Free... */
  free(clim);
  free(met);
  free(dd);

  return EXIT_SUCCESS;
}
