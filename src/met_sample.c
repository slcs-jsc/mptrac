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
  Sample meteorological data at given geolocations.
*/

#include "libtrac.h"

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int const argc,
  char const *argv[]) {

  ctl_t ctl;

  atm_t *atm;

  met_t *met0, *met1;

  FILE *out;

  double h2o, o3, p0, p1, pref, ps, pt, pv, t, tt, uvw[3], z, zm, zref, zt;

  int geopot, ip, it;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <metbase> <atm_in> <sample.tab>");

  /* Allocate... */
  ALLOC(atm, atm_t, 1);
  ALLOC(met0, met_t, 1);
  ALLOC(met1, met_t, 1);

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  geopot =
    (int) scan_ctl(argv[1], argc, argv, "MET_SAMPLE_GEOPOT", -1, "0", NULL);

  /* Read atmospheric data... */
  read_atm(argv[3], &ctl, atm);

  /* Create output file... */
  printf("Write meteorological data file: %s\n", argv[4]);
  if (!(out = fopen(argv[4], "w")))
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
	  "# $17 = tropopause temperature [K]\n\n");

  /* Loop over air parcels... */
  for (ip = 0; ip < atm->np; ip++) {

    /* Get meteorological data... */
    get_met(&ctl, argv[2], atm->time[ip], &met0, &met1);

    /* Set reference pressure for interpolation... */
    pref = atm->p[ip];
    if (geopot) {
      zref = Z(pref);
      p0 = met0->p[0];
      p1 = met0->p[met0->np - 1];
      for (it = 0; it < 24; it++) {
	pref = 0.5 * (p0 + p1);
	intpol_met_time(met0, met1, atm->time[ip], pref, atm->lon[ip], atm->lat[ip],
        NULL, NULL, &zm, NULL, NULL, NULL, NULL);
	if (zref > zm || !gsl_finite(zm))
	  p0 = pref;
	else
	  p1 = pref;
      }
      pref = 0.5 * (p0 + p1);
    }

    /* Interpolate meteorological data... */
    intpol_met_time(met0, met1, atm->time[ip], pref, atm->lon[ip], atm->lat[ip], 
        &ps, &pt, &z, &t, &pv, &h2o, &o3);
    intpol_met_time(met0, met1, atm->time[ip], pt, atm->lon[ip], atm->lat[ip],
        NULL, NULL, &zt, &tt, NULL, NULL, NULL);
    intpol_winds_time(met0, met1, atm->time[ip], pref, atm->lon[ip], atm->lat[ip], uvw);

    /* Write data... */
    fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	    atm->time[ip], Z(atm->p[ip]), atm->lon[ip], atm->lat[ip],
	    atm->p[ip], t, uvw[0], uvw[1], uvw[2], h2o, o3, z, pv, ps, pt, zt, tt);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(atm);
  free(met0);
  free(met1);

  return EXIT_SUCCESS;
}
