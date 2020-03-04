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
  
  Copyright (C) 2013-2019 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Sample meteorological data at given geolocations.
*/

#include "libtrac.h"

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  atm_t *atm;

  met_t *met0, *met1;

  FILE *out;

  double h2o, h2ot, o3, p0, p1, pref, ps, pt, pv, t, tt, u, v, w,
    z, zm, zref, zt, cw[3];

  int geopot, ip, it, ci[3];

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <sample.tab> <metbase> <atm_in>");

  /* Allocate... */
  ALLOC(atm, atm_t, 1);
  ALLOC(met0, met_t, 1);
  ALLOC(met1, met_t, 1);

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  geopot =
    (int) scan_ctl(argv[1], argc, argv, "SAMPLE_GEOPOT", -1, "0", NULL);

  /* Read atmospheric data... */
  if (!read_atm(argv[4], &ctl, atm))
    ERRMSG("Cannot open file!");

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
  fprintf(out,
	  "# $10 = H2O volume mixing ratio [1]\n"
	  "# $11 = O3 volume mixing ratio [1]\n"
	  "# $12 = geopotential height [km]\n"
	  "# $13 = potential vorticity [PVU]\n"
	  "# $14 = surface pressure [hPa]\n"
	  "# $15 = tropopause pressure [hPa]\n"
	  "# $16 = tropopause geopotential height [km]\n"
	  "# $17 = tropopause temperature [K]\n"
	  "# $18 = tropopause water vapor [ppv]\n\n");

  /* Loop over air parcels... */
  for (ip = 0; ip < atm->np; ip++) {

    /* Get meteorological data... */
    get_met(&ctl, argv[3], atm->time[ip], &met0, &met1);

    /* Set reference pressure for interpolation... */
    pref = atm->p[ip];
    if (geopot) {
      zref = Z(pref);
      p0 = met0->p[0];
      p1 = met0->p[met0->np - 1];
      for (it = 0; it < 24; it++) {
	pref = 0.5 * (p0 + p1);
	intpol_met_time_3d(met0, met0->z, met1, met1->z, atm->time[ip], pref,
			   atm->lon[ip], atm->lat[ip], &zm, ci, cw, 1);
	if (zref > zm || !gsl_finite(zm))
	  p0 = pref;
	else
	  p1 = pref;
      }
      pref = 0.5 * (p0 + p1);
    }

    /* Interpolate meteo data... */
    intpol_met_time_3d(met0, met0->z, met1, met1->z, atm->time[ip], pref,
		       atm->lon[ip], atm->lat[ip], &z, ci, cw, 1);
    intpol_met_time_3d(met0, met0->t, met1, met1->t, atm->time[ip], pref,
		       atm->lon[ip], atm->lat[ip], &t, ci, cw, 0);
    intpol_met_time_3d(met0, met0->u, met1, met1->u, atm->time[ip], pref,
		       atm->lon[ip], atm->lat[ip], &u, ci, cw, 0);
    intpol_met_time_3d(met0, met0->v, met1, met1->v, atm->time[ip], pref,
		       atm->lon[ip], atm->lat[ip], &v, ci, cw, 0);
    intpol_met_time_3d(met0, met0->w, met1, met1->w, atm->time[ip], pref,
		       atm->lon[ip], atm->lat[ip], &w, ci, cw, 0);
    intpol_met_time_3d(met0, met0->pv, met1, met1->pv, atm->time[ip], pref,
		       atm->lon[ip], atm->lat[ip], &pv, ci, cw, 0);
    intpol_met_time_3d(met0, met0->h2o, met1, met1->h2o, atm->time[ip], pref,
		       atm->lon[ip], atm->lat[ip], &h2o, ci, cw, 0);
    intpol_met_time_3d(met0, met0->o3, met1, met1->o3, atm->time[ip], pref,
		       atm->lon[ip], atm->lat[ip], &o3, ci, cw, 0);
    intpol_met_time_2d(met0, met0->ps, met1, met1->ps, atm->time[ip],
		       atm->lon[ip], atm->lat[ip], &ps, ci, cw, 0);
    intpol_met_time_2d(met0, met0->pt, met1, met1->pt, atm->time[ip],
		       atm->lon[ip], atm->lat[ip], &pt, ci, cw, 0);

    /* Interpolate tropopause data... */
    intpol_met_time_3d(met0, met0->z, met1, met1->z, atm->time[ip], pt,
		       atm->lon[ip], atm->lat[ip], &zt, ci, cw, 1);
    intpol_met_time_3d(met0, met0->t, met1, met1->t, atm->time[ip], pt,
		       atm->lon[ip], atm->lat[ip], &tt, ci, cw, 0);
    intpol_met_time_3d(met0, met0->h2o, met1, met1->h2o, atm->time[ip], pt,
		       atm->lon[ip], atm->lat[ip], &h2ot, ci, cw, 0);

    /* Write data... */
    fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	    atm->time[ip], Z(atm->p[ip]), atm->lon[ip], atm->lat[ip],
	    atm->p[ip], t, u, v, w, h2o, o3, z, pv, ps, pt, zt, tt, h2ot);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(atm);
  free(met0);
  free(met1);

  return EXIT_SUCCESS;
}
