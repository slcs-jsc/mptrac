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
  
  Copyright (C) 2013-2024 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Sample meteorological data at given geolocations.
*/

#include "mptrac.h"

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  clim_t *clim;

  atm_t *atm;

  met_t *met0, *met1;

  FILE *out;

  double h2o, h2ot, o3, lwc, rwc, iwc, swc, cc, p0, p1, ps, ts, zs, us, vs,
    ess, nss, shf, lsm, sst, pbl, pt, pct, pcb, cl, plcl, plfc, pel,
    cape, cin, o3c, pv, t, tt, u, v, w, z, zm, zref, zt, time_old = -999,
    p_old = -999, lon_old = -999, lat_old = -999;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <sample.tab> <atm_in>");

  /* Allocate... */
  ALLOC(clim, clim_t, 1);
  ALLOC(atm, atm_t, 1);
  ALLOC(met0, met_t, 1);
  ALLOC(met1, met_t, 1);

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  int geopot =
    (int) scan_ctl(argv[1], argc, argv, "SAMPLE_GEOPOT", -1, "0", NULL);
  int grid_time =
    (int) scan_ctl(argv[1], argc, argv, "SAMPLE_GRID_TIME", -1, "0", NULL);
  int grid_z =
    (int) scan_ctl(argv[1], argc, argv, "SAMPLE_GRID_Z", -1, "0", NULL);
  int grid_lon =
    (int) scan_ctl(argv[1], argc, argv, "SAMPLE_GRID_LON", -1, "0", NULL);
  int grid_lat =
    (int) scan_ctl(argv[1], argc, argv, "SAMPLE_GRID_LAT", -1, "0", NULL);

  /* Read climatological data... */
  read_clim(&ctl, clim);

  /* Read atmospheric data... */
  if (!read_atm(argv[3], &ctl, atm))
    ERRMSG("Cannot open file!");

  /* Create output file... */
  LOG(1, "Write meteorological data file: %s", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  MET_HEADER;

  /* Loop over air parcels... */
  for (int ip = 0; ip < atm->np; ip++) {

    /* Get meteorological data... */
    get_met(&ctl, clim, atm->time[ip], &met0, &met1);

    /* Set reference pressure for interpolation... */
    INTPOL_INIT;
    double pref = atm->p[ip];
    if (geopot) {
      zref = Z(pref);
      p0 = met0->p[0];
      p1 = met0->p[met0->np - 1];
      for (int it = 0; it < 24; it++) {
	pref = 0.5 * (p0 + p1);
	intpol_met_time_3d(met0, met0->z, met1, met1->z, atm->time[ip], pref,
			   atm->lon[ip], atm->lat[ip], &zm, ci, cw, 1);
	if (zref > zm || !isfinite(zm))
	  p0 = pref;
	else
	  p1 = pref;
      }
      pref = 0.5 * (p0 + p1);
    }

    /* Interpolate meteo data... */
    INTPOL_TIME_ALL(atm->time[ip], pref, atm->lon[ip], atm->lat[ip]);

    /* Make blank lines... */
    if (ip == 0 || (grid_time && atm->time[ip] != time_old)
	|| (grid_z && atm->p[ip] != p_old)
	|| (grid_lon && atm->lon[ip] != lon_old)
	|| (grid_lat && atm->lat[ip] != lat_old))
      fprintf(out, "\n");
    time_old = atm->time[ip];
    p_old = atm->p[ip];
    lon_old = atm->lon[ip];
    lat_old = atm->lat[ip];

    /* Write data... */
    fprintf(out,
	    "%.2f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g"
	    " %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g"
	    " %g %g %g %g %g %g %g %g %g %g %g %g %g %g 1 1 1\n",
	    atm->time[ip], Z(atm->p[ip]), atm->lon[ip], atm->lat[ip],
	    atm->p[ip], t, u, v, w, h2o, o3, z, pv, ps, ts, zs, us, vs,
	    ess, nss, shf, lsm,
	    sst, pt, zt, tt, h2ot, lwc, rwc, iwc, swc, cc, cl, pct, pcb, plcl,
	    plfc, pel, cape, cin, RH(atm->p[ip], t, h2o), RHICE(atm->p[ip], t,
								h2o),
	    TDEW(atm->p[ip], h2o), TICE(atm->p[ip], h2o),
	    nat_temperature(atm->p[ip], h2o,
			    clim_zm(&clim->hno3, atm->time[ip], atm->lat[ip],
				    atm->p[ip])), clim_zm(&clim->hno3,
							  atm->time[ip],
							  atm->lat[ip],
							  atm->p[ip]),
	    clim_oh(&ctl, clim, atm->time[ip], atm->lon[ip], atm->lat[ip],
		    atm->p[ip]), clim_zm(&clim->h2o2, atm->time[ip],
					 atm->lat[ip], atm->p[ip]),
	    clim_zm(&clim->ho2, atm->time[ip], atm->lat[ip], atm->p[ip]),
	    clim_zm(&clim->o1d, atm->time[ip], atm->lat[ip], atm->p[ip]), pbl,
	    o3c);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(clim);
  free(atm);
  free(met0);
  free(met1);

  return EXIT_SUCCESS;
}
