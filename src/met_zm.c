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

/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

/*! Calculate zonal mean. */
#define MEAN(x)					\
  (x[ip][iy] / np[ip][iy])

/*! Calculate standard deviation. */
#define SIGMA(x, x2)							\
  (np[ix][iy] > 1 ? sqrt(x2[ip][iy] / np[ip][iy] -			\
			 gsl_pow_2(x[ip][iy] / np[ip][iy])) : 0)

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  met_t *met;

  FILE *out;

  static double timem[EP][EY], psm[EP][EY], ptm[EP][EY], tm[EP][EY],
    um[EP][EY], vm[EP][EY], vhm[EP][EY], wm[EP][EY], h2om[EP][EY],
    o3m[EP][EY], zm[EP][EY], psm2[EP][EY], ptm2[EP][EY], tm2[EP][EY],
    um2[EP][EY], vm2[EP][EY], vhm2[EP][EY], wm2[EP][EY], h2om2[EP][EY],
    o3m2[EP][EY], zm2[EP][EY];

  static int i, ip, ix, iy, np[EP][EY];

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
	  timem[ip][iy] += met->time;
	  psm[ip][iy] += met->ps[ix][iy];
	  psm2[ip][iy] += gsl_pow_2(met->ps[ix][iy]);
	  ptm[ip][iy] += met->pt[ix][iy];
	  ptm2[ip][iy] += gsl_pow_2(met->pt[ix][iy]);
	  zm[ip][iy] += met->z[ix][iy][ip];
	  zm2[ip][iy] += gsl_pow_2(met->z[ix][iy][ip]);
	  tm[ip][iy] += met->t[ix][iy][ip];
	  tm2[ip][iy] += gsl_pow_2(met->t[ix][iy][ip]);
	  um[ip][iy] += met->u[ix][iy][ip];
	  um2[ip][iy] += gsl_pow_2(met->u[ix][iy][ip]);
	  vm[ip][iy] += met->v[ix][iy][ip];
	  vm2[ip][iy] += gsl_pow_2(met->v[ix][iy][ip]);
	  vhm[ip][iy] += sqrt(gsl_pow_2(met->u[ix][iy][ip])
			      + gsl_pow_2(met->v[ix][iy][ip]));
	  vhm2[ip][iy] += gsl_pow_2(met->u[ix][iy][ip])
	    + gsl_pow_2(met->v[ix][iy][ip]);
	  wm[ip][iy] += met->w[ix][iy][ip];
	  wm2[ip][iy] += gsl_pow_2(met->w[ix][iy][ip]);
	  h2om[ip][iy] += met->h2o[ix][iy][ip];
	  h2om2[ip][iy] += gsl_pow_2(met->h2o[ix][iy][ip]);
	  o3m[ip][iy] += met->o3[ix][iy][ip];
	  o3m2[ip][iy] += gsl_pow_2(met->o3[ix][iy][ip]);
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
	  "# $9  = meridional wind standard deviation [m/s]\n");
  fprintf(out,
	  "# $10 = horizontal wind mean [m/s]\n"
	  "# $11 = horizontal wind standard deviation [m/s]\n"
	  "# $12 = vertical wind mean [hPa/s]\n"
	  "# $13 = vertical wind standard deviation [hPa/s]\n"
	  "# $14 = H2O vmr mean [1]\n"
	  "# $15 = H2O vmr standard deviation [1]\n"
	  "# $16 = O3 vmr mean [1]\n"
	  "# $17 = O3 vmr standard deviation [1]\n"
	  "# $18 = geopotential height mean [hPa]\n"
	  "# $19 = geopotential height standard deviation [hPa]\n");
  fprintf(out,
	  "# $20 = surface pressure mean [hPa]\n"
	  "# $21 = surface pressure standard deviation [hPa]\n"
	  "# $22 = tropopause pressure mean [hPa]\n"
	  "# $23 = tropopause pressure standard deviation [hPa]\n");

  /* Write data... */
  for (iy = 0; iy < met->ny; iy++) {
    fprintf(out, "\n");
    for (ip = 0; ip < met->np; ip++)
      fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %g"
	      " %g %g %g %g %g %g %g %g %g %g %g %g\n",
	      MEAN(timem), Z(met->p[ip]), met->lat[iy],
	      MEAN(tm), SIGMA(tm, tm2), MEAN(um), SIGMA(um, um2),
	      MEAN(vm), SIGMA(vm, vm2), MEAN(vhm), SIGMA(vhm, vhm2),
	      MEAN(wm), SIGMA(wm, wm2), MEAN(h2om), SIGMA(h2om, h2om2),
	      MEAN(o3m), SIGMA(o3m, o3m2), MEAN(zm), SIGMA(zm, zm2),
	      MEAN(psm), SIGMA(psm, psm2), MEAN(ptm), SIGMA(ptm, ptm2));
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(met);

  return EXIT_SUCCESS;
}
