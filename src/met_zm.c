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
  
  Copyright (C) 2013-2018 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Extract zonal mean from meteorological data.
*/

#include "libtrac.h"

int main(
  int const argc,
  char const *argv[]) {

  ctl_t ctl;

  met_t *met;

  FILE *out;

  static double timem[EP][EY], psm[EP][EY], ptm[EP][EY], ttm[EP][EY],
    ztm[EP][EY], Tm[EP][EY], um[EP][EY], vm[EP][EY], wm[EP][EY], h2om[EP][EY],
    pvm[EP][EY], o3m[EP][EY], zm[EP][EY], zt, tt;

  static int i, ip, ix, iy, np[EP][EY], npt[EP][EY];

  /* Allocate... */
  ALLOC(met, met_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <zm.tab> <met0> [ <met1> ... ]");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);

  /* Loop over files... */
  for (i = 3; i < argc; i++) {

    /* Read meteorological data... */
    read_met(&ctl, argv[i], met);

    /* Average data... */
    for (ix = 0; ix < met->nx; ix++)
      for (iy = 0; iy < met->ny; iy++)
	for (ip = 0; ip < met->np; ip++) {
	  intpol_met_space(met, met->pt[ix][iy], met->lon[ix], met->lat[iy],
			   NULL, NULL, &zt, &tt, NULL, NULL, NULL, NULL, NULL,
			   NULL);
	  timem[ip][iy] += met->time;
	  zm[ip][iy] += met->z[ix][iy][ip];
	  Tm[ip][iy] += met->T[ix][iy][ip];
	  um[ip][iy] += met->u[ix][iy][ip];
	  vm[ip][iy] += met->v[ix][iy][ip];
	  wm[ip][iy] += met->w[ix][iy][ip];
	  pvm[ip][iy] += met->pv[ix][iy][ip];
	  h2om[ip][iy] += met->h2o[ix][iy][ip];
	  o3m[ip][iy] += met->o3[ix][iy][ip];
	  psm[ip][iy] += met->ps[ix][iy];
	  if (gsl_finite(met->pt[ix][iy])) {
	    ptm[ip][iy] += met->pt[ix][iy];
	    ztm[ip][iy] += zt;
	    ttm[ip][iy] += tt;
	    npt[ip][iy]++;
	  }
	  np[ip][iy]++;
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
	  "# $10 = H2O volume mixing ratio [1]\n");
  fprintf(out,
	  "# $11 = O3 volume mixing ratio [1]\n"
	  "# $12 = geopotential height [km]\n"
	  "# $13 = potential vorticity [PVU]\n"
	  "# $14 = surface pressure [hPa]\n"
	  "# $15 = tropopause pressure [hPa]\n"
	  "# $16 = tropopause geopotential height [km]\n"
	  "# $17 = tropopause temperature [K]\n");

  /* Write data... */
  for (ip = 0; ip < met->np; ip++) {
    fprintf(out, "\n");
    for (iy = 0; iy < met->ny; iy++)
      fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	      timem[ip][iy] / np[ip][iy], Z(met->p[ip]), 0.0, met->lat[iy],
	      met->p[ip], Tm[ip][iy] / np[ip][iy], um[ip][iy] / np[ip][iy],
	      vm[ip][iy] / np[ip][iy], wm[ip][iy] / np[ip][iy],
	      h2om[ip][iy] / np[ip][iy], o3m[ip][iy] / np[ip][iy],
	      zm[ip][iy] / np[ip][iy], pvm[ip][iy] / np[ip][iy],
	      psm[ip][iy] / np[ip][iy], ptm[ip][iy] / npt[ip][iy],
	      ztm[ip][iy] / npt[ip][iy], ttm[ip][iy] / npt[ip][iy]);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(met);

  return EXIT_SUCCESS;
}
