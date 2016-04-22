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
  
  Copright (C) 2013-2015 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Extract global map from meteorological data.
*/

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  met_t *met;

  FILE *in, *out;

  static double dz, dzmin = 1e10, z, timem[EX][EY], psm[EX][EY], tm[EX][EY],
    um[EX][EY], vm[EX][EY], wm[EX][EY], h2om[EX][EY], o3m[EX][EY];

  static int i, ip, ip2, ix, iy, np[EX][EY];

  /* Allocate... */
  ALLOC(met, met_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <map.tab> <met0> [ <met1> ... ]");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  z = scan_ctl(argv[1], argc, argv, "Z", -1, "", NULL);

  /* Loop over files... */
  for (i = 3; i < argc; i++) {

    /* Read meteorological data... */
    if (!(in = fopen(argv[i], "r")))
      continue;
    else
      fclose(in);
    read_met(argv[i], met);

    /* Find nearest pressure level... */
    for (ip2 = 0; ip2 < met->np; ip2++) {
      dz = fabs(Z(met->p[ip2]) - z);
      if (dz < dzmin) {
	dzmin = dz;
	ip = ip2;
      }
    }

    /* Average data... */
    for (ix = 0; ix < met->nx; ix++)
      for (iy = 0; iy < met->ny; iy++) {
	timem[ix][iy] += met->time;
	tm[ix][iy] += met->t[ix][iy][ip];
	um[ix][iy] += met->u[ix][iy][ip];
	vm[ix][iy] += met->v[ix][iy][ip];
	wm[ix][iy] += met->w[ix][iy][ip];
	h2om[ix][iy] += met->h2o[ix][iy][ip];
	o3m[ix][iy] += met->o3[ix][iy][ip];
	psm[ix][iy] += met->ps[ix][iy];
	np[ix][iy]++;
      }
  }

  /* Create output file... */
  printf("Write meteorological data file: %s\n", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1  = time [s]\n"
	  "# $2  = altitude [km]\n"
	  "# $3  = longitude [deg]\n"
	  "# $4  = latitude [deg]\n"
	  "# $5  = pressure [hPa]\n"
	  "# $6  = temperature [K]\n"
	  "# $7  = zonal wind [m/s]\n"
	  "# $8  = meridional wind [m/s]\n"
	  "# $9  = vertical wind [hPa/s]\n"
	  "# $10 = H2O volume mixing ratio [1]\n"
	  "# $11 = O3 volume mixing ratio [1]\n"
	  "# $12 = surface pressure [hPa]\n");

  /* Write data... */
  for (iy = 0; iy < met->ny; iy++) {
    fprintf(out, "\n");
    for (ix = 0; ix < met->nx; ix++)
      if (met->lon[ix] >= 180)
	fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %g %g\n",
		timem[ix][iy] / np[ix][iy], Z(met->p[ip]),
		met->lon[ix] - 360.0, met->lat[iy], met->p[ip],
		tm[ix][iy] / np[ix][iy], um[ix][iy] / np[ix][iy],
		vm[ix][iy] / np[ix][iy], wm[ix][iy] / np[ix][iy],
		h2om[ix][iy] / np[ix][iy], o3m[ix][iy] / np[ix][iy],
		psm[ix][iy] / np[ix][iy]);
    for (ix = 0; ix < met->nx; ix++)
      if (met->lon[ix] <= 180)
	fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %g %g\n",
		timem[ix][iy] / np[ix][iy], Z(met->p[ip]),
		met->lon[ix], met->lat[iy], met->p[ip],
		tm[ix][iy] / np[ix][iy], um[ix][iy] / np[ix][iy],
		vm[ix][iy] / np[ix][iy], wm[ix][iy] / np[ix][iy],
		h2om[ix][iy] / np[ix][iy], o3m[ix][iy] / np[ix][iy],
		psm[ix][iy] / np[ix][iy]);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(met);

  return EXIT_SUCCESS;
}
