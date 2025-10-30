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
  Check model time step based on given meteorological data.
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
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <dt_file> <met>");

  /* Allocate... */
  ALLOC(clim, clim_t, 1);
  ALLOC(met, met_t, 1);
  ALLOC(dd, dd_t, 1);

  /* Read control parameters... */
  mptrac_read_ctl(argv[1], argc, argv, &ctl);
  const double kx = scan_ctl(argv[1], argc, argv, "KX", -1, "50.0", NULL);
  const double kz = scan_ctl(argv[1], argc, argv, "KZ", -1, "0.1", NULL);
  const double dx = 1e3 * scan_ctl(argv[1], argc, argv, "DX", -1, "", NULL);
  const double c_max = scan_ctl(argv[1], argc, argv, "CMAX", -1, "0.5", NULL);
  const double n_max = scan_ctl(argv[1], argc, argv, "NMAX", -1, "0.3", NULL);

  /* Read climatological data... */
  mptrac_read_clim(&ctl, clim);

  /* Read meteo data... */
  if (!mptrac_read_met(argv[3], &ctl, clim, met, dd))
    ERRMSG("Cannot open file!");

  /* Create output file... */
  FILE *out;
  LOG(1, "Write time step data file: %s", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = height [km]\n"
	  "# $2 = time step for horizontal advection [s]\n"
	  "# $3 = time step for vertical advection [s]\n"
	  "# $4 = time step for horizontal diffusion [s]\n"
	  "# $5 = time step for vertical diffusion [s]\n\n");

  /* Loop over pressure levels... */
  for (int ip = 1; ip < met->np - 1; ip++) {

    /* Init... */
    double dt_x_min = 1e100;
    double dt_p_min = 1e100;
    double dt_dx_min = 1e100;
    double dt_dp_min = 1e100;

    /* Loop over columns... */
#pragma omp parallel for default(shared) collapse(2) reduction(min:dt_x_min,dt_p_min,dt_dx_min,dt_dp_min)
    for (int ix = 0; ix < met->nx; ix++)
      for (int iy = 1; iy < met->ny - 1; iy++) {

	/* Check advection... */
	const double vh =
	  sqrt(SQR(met->u[ix][iy][ip]) + SQR(met->v[ix][iy][ip]));
	const double dt_x = fabs(c_max * dx / vh);
	if (vh != 0)
	  dt_x_min = MIN(dt_x, dt_x_min);

	const double dp = 0.5 * fabs(met->p[ip + 1] - met->p[ip - 1]);
	const double dt_p = fabs(c_max * dp / met->w[ix][iy][ip]);
	if (met->w[ix][iy][ip])
	  dt_p_min = MIN(dt_p, dt_p_min);

	/* Check diffusion... */
	const double dt_dx = 0.5 * SQR(n_max * dx) / kx;
	dt_dx_min = MIN(dt_dx, dt_dx_min);

	const double dt_dp =
	  0.5 * SQR(n_max * dp) / (SQR(met->p[ip] / (100. * H0)) * kz);
	dt_dp_min = MIN(dt_dp, dt_dp_min);
      }

    /* Write output... */
    fprintf(out, "%g %g %g %g %g\n", Z(met->p[ip]), dt_x_min, dt_p_min,
	    dt_dx_min, dt_dp_min);
  }

  /* Close output file... */
  fclose(out);

  /* Free... */
  free(clim);
  free(met);
  free(dd);

  return EXIT_SUCCESS;
}
