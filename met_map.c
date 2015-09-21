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

  static double dz, dzmin = 1e10, z, timem[EX][EY], tm[EX][EY],
    um[EX][EY], vm[EX][EY], wm[EX][EY];

  static int epol, i, ip, ip2, ix, iy, np[EX][EY], redx, redy, redp, wrap;

  /* Allocate... */
  ALLOC(met, met_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <map.tab> <met0> [ <met1> ... ]");

  /* Read control parameters... */
  read_ctl(NULL, NULL, argc, argv, &ctl);
  z = scan_ctl(NULL, NULL, argc, argv, "Z", -1, "", NULL);
  epol = (int) scan_ctl(NULL, NULL, argc, argv, "EPOL", -1, "0", NULL);
  redx = (int) scan_ctl(NULL, NULL, argc, argv, "REDX", -1, "1", NULL);
  redy = (int) scan_ctl(NULL, NULL, argc, argv, "REDY", -1, "1", NULL);
  redp = (int) scan_ctl(NULL, NULL, argc, argv, "REDP", -1, "1", NULL);
  wrap = (int) scan_ctl(NULL, NULL, argc, argv, "WRAP", -1, "0", NULL);

  /* Loop over files... */
  for (i = 3; i < argc; i++) {

    /* Read meteorological data... */
    if (!(in = fopen(argv[i], "r")))
      continue;
    else
      fclose(in);
    read_met(argv[i], met);

    /* Extrapolate... */
    if (epol)
      extrapolate_met(met);

    /* Reduce resolution... */
    reduce_met(met, redx, redy, redp);

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
	  "# $9  = vertical wind [hPa/s]\n");

  /* Write data... */
  for (iy = 0; iy < met->ny; iy++) {
    fprintf(out, "\n");
    if (wrap) {
      for (ix = 0; ix < met->nx; ix++)
	if (met->lon[ix] >= 180)
	  fprintf(out, "%.2f %g %g %g %g %g %g %g %g\n",
		  timem[ix][iy] / np[ix][iy], Z(met->p[ip]),
		  met->lon[ix] - 360.0, met->lat[iy], met->p[ip],
		  tm[ix][iy] / np[ix][iy], um[ix][iy] / np[ix][iy],
		  vm[ix][iy] / np[ix][iy], wm[ix][iy] / np[ix][iy]);
      for (ix = 0; ix < met->nx; ix++)
	if (met->lon[ix] <= 180)
	  fprintf(out, "%.2f %g %g %g %g %g %g %g %g\n",
		  timem[ix][iy] / np[ix][iy], Z(met->p[ip]),
		  met->lon[ix], met->lat[iy], met->p[ip],
		  tm[ix][iy] / np[ix][iy], um[ix][iy] / np[ix][iy],
		  vm[ix][iy] / np[ix][iy], wm[ix][iy] / np[ix][iy]);
    } else {
      for (ix = 0; ix < met->nx; ix++)
	fprintf(out, "%.2f %g %g %g %g %g %g %g %g\n",
		timem[ix][iy] / np[ix][iy], Z(met->p[ip]),
		met->lon[ix], met->lat[iy], met->p[ip],
		tm[ix][iy] / np[ix][iy], um[ix][iy] / np[ix][iy],
		vm[ix][iy] / np[ix][iy], wm[ix][iy] / np[ix][iy]);
    }
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(met);

  return EXIT_SUCCESS;
}
