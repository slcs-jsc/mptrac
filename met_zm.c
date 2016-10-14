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
  Extract zonal mean from meteorological data.
*/

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  met_t *met;

  FILE *in, *out;

  static double timem[EP][EY], psm[EP][EY], tm[EP][EY],
    um[EP][EY], vm[EP][EY], wm[EP][EY], h2om[EP][EY], o3m[EP][EY],
    psm2[EP][EY], tm2[EP][EY], um2[EP][EY], vm2[EP][EY], wm2[EP][EY],
    h2om2[EP][EY], o3m2[EP][EY];

  static int i, ip, ix, iy, np[EP][EY];

  /* Allocate... */
  ALLOC(met, met_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <map.tab> <met0> [ <met1> ... ]");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);

  /* Loop over files... */
  for (i = 3; i < argc; i++) {

    /* Read meteorological data... */
    if (!(in = fopen(argv[i], "r")))
      continue;
    else
      fclose(in);
    read_met(argv[i], met);

    /* Average data... */
    for (ix = 0; ix < met->nx; ix++)
      for (iy = 0; iy < met->ny; iy++)
	for (ip = 0; ip < met->np; ip++) {
	  timem[ip][iy] += met->time;
	  tm[ip][iy] += met->t[ix][iy][ip];
	  um[ip][iy] += met->u[ix][iy][ip];
	  vm[ip][iy] += met->v[ix][iy][ip];
	  wm[ip][iy] += met->w[ix][iy][ip];
	  h2om[ip][iy] += met->h2o[ix][iy][ip];
	  o3m[ip][iy] += met->o3[ix][iy][ip];
	  psm[ip][iy] += met->ps[ix][iy];
	  tm2[ip][iy] += gsl_pow_2(met->t[ix][iy][ip]);
	  um2[ip][iy] += gsl_pow_2(met->u[ix][iy][ip]);
	  vm2[ip][iy] += gsl_pow_2(met->v[ix][iy][ip]);
	  wm2[ip][iy] += gsl_pow_2(met->w[ix][iy][ip]);
	  h2om2[ip][iy] += gsl_pow_2(met->h2o[ix][iy][ip]);
	  o3m2[ip][iy] += gsl_pow_2(met->o3[ix][iy][ip]);
	  psm2[ip][iy] += gsl_pow_2(met->ps[ix][iy]);
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
	  "# $3  = latitude [deg]\n"
	  "# $4  = temperature mean [K]\n"
	  "# $5  = temperature standard deviation [K]\n"
	  "# $6  = zonal wind mean [m/s]\n"
	  "# $7  = zonal wind standard deviation [m/s]\n"
	  "# $8  = meridional wind mean [m/s]\n"
	  "# $9  = meridional wind standard deviation [m/s]\n"
	  "# $10 = vertical wind mean [hPa/s]\n"
	  "# $11 = vertical wind standard deviation [hPa/s]\n"
	  "# $12 = H2O vmr mean [1]\n"
	  "# $13 = H2O vmr standard deviation [1]\n"
	  "# $14 = O3 vmr mean [1]\n"
	  "# $15 = O3 vmr standard deviation [1]\n"
	  "# $16 = surface pressure mean [hPa]\n"
	  "# $17 = surface pressure standard deviation [hPa]\n");

  /* Write data... */
  for (iy = 0; iy < met->ny; iy++) {
    fprintf(out, "\n");
    for (ip = 0; ip < met->np; ip++)
      fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	      timem[ip][iy] / np[ip][iy], Z(met->p[ip]), met->lat[iy],
	      tm[ip][iy] / np[ip][iy],
	      sqrt(tm2[ip][iy] / np[ip][iy] -
		   gsl_pow_2(tm[ip][iy] / np[ip][iy])),
	      um[ip][iy] / np[ip][iy],
	      sqrt(um2[ip][iy] / np[ip][iy] -
		   gsl_pow_2(um[ip][iy] / np[ip][iy])),
	      vm[ip][iy] / np[ip][iy],
	      sqrt(vm2[ip][iy] / np[ip][iy] -
		   gsl_pow_2(vm[ip][iy] / np[ip][iy])),
	      wm[ip][iy] / np[ip][iy],
	      sqrt(wm2[ip][iy] / np[ip][iy] -
		   gsl_pow_2(wm[ip][iy] / np[ip][iy])),
	      h2om[ip][iy] / np[ip][iy],
	      sqrt(h2om2[ip][iy] / np[ip][iy] -
		   gsl_pow_2(h2om[ip][iy] / np[ip][iy])),
	      o3m[ip][iy] / np[ip][iy],
	      sqrt(o3m2[ip][iy] / np[ip][iy] -
		   gsl_pow_2(o3m[ip][iy] / np[ip][iy])),
	      psm[ip][iy] / np[ip][iy],
	      sqrt(psm2[ip][iy] / np[ip][iy] -
		   gsl_pow_2(psm[ip][iy] / np[ip][iy])));
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(met);

  return EXIT_SUCCESS;
}
