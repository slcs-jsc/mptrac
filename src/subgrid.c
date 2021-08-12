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
  Calculate standard deviations of horizontal wind and vertical velocity.
*/

#include "libtrac.h"

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  met_t *met0, *met1;

  FILE *out;

  static double usig[EP][EY], vsig[EP][EY], wsig[EP][EY], u[16], v[16], w[16];

  static int i, ix, iy, iz, n[EP][EY];

  /* Allocate... */
  ALLOC(met0, met_t, 1);
  ALLOC(met1, met_t, 1);

  /* Check arguments... */
  if (argc < 4 && argc % 2 != 0)
    ERRMSG
      ("Give parameters: <ctl> <zm.tab> <met0> <met1> [ <met0> <met1> ... ]");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);

  /* Loop over data files... */
  for (i = 3; i < argc - 1; i += 2) {

    /* Read meteorological data... */
    if (!read_met(&ctl, argv[i], met0))
      ERRMSG("Cannot open file!");
    if (!read_met(&ctl, argv[i + 1], met1))
      ERRMSG("Cannot open file!");

    /* Loop over grid boxes... */
    for (ix = 0; ix < met0->nx - 1; ix++)
      for (iy = 0; iy < met0->ny - 1; iy++)
	for (iz = 0; iz < met0->np - 1; iz++) {

	  /* Collect local wind data... */
	  u[0] = met0->u[ix][iy][iz];
	  u[1] = met0->u[ix + 1][iy][iz];
	  u[2] = met0->u[ix][iy + 1][iz];
	  u[3] = met0->u[ix + 1][iy + 1][iz];
	  u[4] = met0->u[ix][iy][iz + 1];
	  u[5] = met0->u[ix + 1][iy][iz + 1];
	  u[6] = met0->u[ix][iy + 1][iz + 1];
	  u[7] = met0->u[ix + 1][iy + 1][iz + 1];

	  v[0] = met0->v[ix][iy][iz];
	  v[1] = met0->v[ix + 1][iy][iz];
	  v[2] = met0->v[ix][iy + 1][iz];
	  v[3] = met0->v[ix + 1][iy + 1][iz];
	  v[4] = met0->v[ix][iy][iz + 1];
	  v[5] = met0->v[ix + 1][iy][iz + 1];
	  v[6] = met0->v[ix][iy + 1][iz + 1];
	  v[7] = met0->v[ix + 1][iy + 1][iz + 1];

	  w[0] = 1e3 * DP2DZ(met0->w[ix][iy][iz], met0->p[iz]);
	  w[1] = 1e3 * DP2DZ(met0->w[ix + 1][iy][iz], met0->p[iz]);
	  w[2] = 1e3 * DP2DZ(met0->w[ix][iy + 1][iz], met0->p[iz]);
	  w[3] = 1e3 * DP2DZ(met0->w[ix + 1][iy + 1][iz], met0->p[iz]);
	  w[4] = 1e3 * DP2DZ(met0->w[ix][iy][iz + 1], met0->p[iz + 1]);
	  w[5] = 1e3 * DP2DZ(met0->w[ix + 1][iy][iz + 1], met0->p[iz + 1]);
	  w[6] = 1e3 * DP2DZ(met0->w[ix][iy + 1][iz + 1], met0->p[iz + 1]);
	  w[7] =
	    1e3 * DP2DZ(met0->w[ix + 1][iy + 1][iz + 1], met0->p[iz + 1]);

	  /* Collect local wind data... */
	  u[8] = met1->u[ix][iy][iz];
	  u[9] = met1->u[ix + 1][iy][iz];
	  u[10] = met1->u[ix][iy + 1][iz];
	  u[11] = met1->u[ix + 1][iy + 1][iz];
	  u[12] = met1->u[ix][iy][iz + 1];
	  u[13] = met1->u[ix + 1][iy][iz + 1];
	  u[14] = met1->u[ix][iy + 1][iz + 1];
	  u[15] = met1->u[ix + 1][iy + 1][iz + 1];

	  v[8] = met1->v[ix][iy][iz];
	  v[9] = met1->v[ix + 1][iy][iz];
	  v[10] = met1->v[ix][iy + 1][iz];
	  v[11] = met1->v[ix + 1][iy + 1][iz];
	  v[12] = met1->v[ix][iy][iz + 1];
	  v[13] = met1->v[ix + 1][iy][iz + 1];
	  v[14] = met1->v[ix][iy + 1][iz + 1];
	  v[15] = met1->v[ix + 1][iy + 1][iz + 1];

	  w[8] = 1e3 * DP2DZ(met1->w[ix][iy][iz], met1->p[iz]);
	  w[9] = 1e3 * DP2DZ(met1->w[ix + 1][iy][iz], met1->p[iz]);
	  w[10] = 1e3 * DP2DZ(met1->w[ix][iy + 1][iz], met1->p[iz]);
	  w[11] = 1e3 * DP2DZ(met1->w[ix + 1][iy + 1][iz], met1->p[iz]);
	  w[12] = 1e3 * DP2DZ(met1->w[ix][iy][iz + 1], met1->p[iz + 1]);
	  w[13] = 1e3 * DP2DZ(met1->w[ix + 1][iy][iz + 1], met1->p[iz + 1]);
	  w[14] = 1e3 * DP2DZ(met1->w[ix][iy + 1][iz + 1], met1->p[iz + 1]);
	  w[15] =
	    1e3 * DP2DZ(met1->w[ix + 1][iy + 1][iz + 1], met1->p[iz + 1]);

	  /* Get standard deviations of local wind data... */
	  usig[iz][iy] += stddev(u, 16);
	  vsig[iz][iy] += stddev(v, 16);
	  wsig[iz][iy] += stddev(w, 16);
	  n[iz][iy]++;

	  /* Check surface pressure... */
	  if (met0->p[iz] > met0->ps[ix][iy]
	      || met1->p[iz] > met1->ps[ix][iy]) {
	    usig[iz][iy] = GSL_NAN;
	    vsig[iz][iy] = GSL_NAN;
	    wsig[iz][iy] = GSL_NAN;
	    n[iz][iy] = 0;
	  }
	}
  }

  /* Create output file... */
  LOG(1, "Write subgrid data file: %s", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = altitude [km]\n"
	  "# $3 = longitude [deg]\n"
	  "# $4 = latitude [deg]\n"
	  "# $5 = zonal wind standard deviation [m/s]\n"
	  "# $6 = meridional wind standard deviation [m/s]\n"
	  "# $7 = vertical velocity standard deviation [m/s]\n"
	  "# $8 = number of data points\n");

  /* Write output... */
  for (iy = 0; iy < met0->ny - 1; iy++) {
    fprintf(out, "\n");
    for (iz = 0; iz < met0->np - 1; iz++)
      fprintf(out, "%.2f %g %g %g %g %g %g %d\n",
	      0.5 * (met0->time + met1->time),
	      0.5 * (Z(met0->p[iz]) + Z(met1->p[iz + 1])),
	      0.0, 0.5 * (met0->lat[iy] + met1->lat[iy + 1]),
	      usig[iz][iy] / n[iz][iy], vsig[iz][iy] / n[iz][iy],
	      wsig[iz][iy] / n[iz][iy], n[iz][iy]);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(met0);
  free(met1);

  return EXIT_SUCCESS;
}
