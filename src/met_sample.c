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

  double h2o, h2ot, o3, lwc, iwc, p0, p1, ps, ts, zs, us, vs, pbl, pt,
    pct, pcb, cl, plcl, plfc, pel, cape, cin, pv, t, tt, u, v, w, z, zm, zref,
    zt, cw[3], time_old = -999, p_old = -999, lon_old = -999, lat_old = -999;

  int geopot, grid_time, grid_z, grid_lon, grid_lat, ip, it, ci[3];

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <sample.tab> <atm_in>");

  /* Allocate... */
  ALLOC(atm, atm_t, 1);
  ALLOC(met0, met_t, 1);
  ALLOC(met1, met_t, 1);

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  geopot =
    (int) scan_ctl(argv[1], argc, argv, "SAMPLE_GEOPOT", -1, "0", NULL);
  grid_time =
    (int) scan_ctl(argv[1], argc, argv, "SAMPLE_GRID_TIME", -1, "0", NULL);
  grid_z =
    (int) scan_ctl(argv[1], argc, argv, "SAMPLE_GRID_Z", -1, "0", NULL);
  grid_lon =
    (int) scan_ctl(argv[1], argc, argv, "SAMPLE_GRID_LON", -1, "0", NULL);
  grid_lat =
    (int) scan_ctl(argv[1], argc, argv, "SAMPLE_GRID_LAT", -1, "0", NULL);

  /* Read atmospheric data... */
  if (!read_atm(argv[3], &ctl, atm))
    ERRMSG("Cannot open file!");

  /* Create output file... */
  LOG(1, "Write meteorological data file: %s", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = altitude [km]\n"
	  "# $3 = longitude [deg]\n"
	  "# $4 = latitude [deg]\n"
	  "# $5 = pressure [hPa]\n"
	  "# $6 = temperature [K]\n"
	  "# $7 = zonal wind [m/s]\n"
	  "# $8 = meridional wind [m/s]\n"
	  "# $9 = vertical velocity [hPa/s]\n"
	  "# $10 = H2O volume mixing ratio [ppv]\n");
  fprintf(out,
	  "# $11 = O3 volume mixing ratio [ppv]\n"
	  "# $12 = geopotential height [km]\n"
	  "# $13 = potential vorticity [PVU]\n"
	  "# $14 = surface pressure [hPa]\n"
	  "# $15 = surface temperature [K]\n"
	  "# $16 = surface geopotential height [km]\n"
	  "# $17 = surface zonal wind [m/s]\n"
	  "# $18 = surface meridional wind [m/s]\n"
	  "# $19 = tropopause pressure [hPa]\n"
	  "# $20 = tropopause geopotential height [km]\n");
  fprintf(out,
	  "# $21 = tropopause temperature [K]\n"
	  "# $22 = tropopause water vapor [ppv]\n"
	  "# $23 = cloud liquid water content [kg/kg]\n"
	  "# $24 = cloud ice water content [kg/kg]\n"
	  "# $25 = total column cloud water [kg/m^2]\n"
	  "# $26 = cloud top pressure [hPa]\n"
	  "# $27 = cloud bottom pressure [hPa]\n"
	  "# $28 = pressure at lifted condensation level (LCL) [hPa]\n"
	  "# $29 = pressure at level of free convection (LFC) [hPa]\n"
	  "# $30 = pressure at equilibrium level (EL) [hPa]\n");
  fprintf(out,
	  "# $31 = convective available potential energy (CAPE) [J/kg]\n"
	  "# $32 = convective inhibition (CIN) [J/kg]\n"
	  "# $33 = relative humidity over water [%%]\n"
	  "# $34 = relative humidity over ice [%%]\n"
	  "# $35 = dew point temperature [K]\n"
	  "# $36 = frost point temperature [K]\n"
	  "# $37 = NAT temperature [K]\n"
	  "# $38 = HNO3 volume mixing ratio [ppv]\n"
	  "# $39 = OH concentration [molec/cm^3]\n"
	  "# $40 = boundary layer pressure [hPa]\n");

  /* Loop over air parcels... */
  for (ip = 0; ip < atm->np; ip++) {

    /* Get meteorological data... */
    get_met(&ctl, atm->time[ip], &met0, &met1);

    /* Set reference pressure for interpolation... */
    double pref = atm->p[ip];
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
	    "%.2f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g"
	    " %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	    atm->time[ip], Z(atm->p[ip]), atm->lon[ip], atm->lat[ip],
	    atm->p[ip], t, u, v, w, h2o, o3, z, pv, ps, ts, zs, us, vs,
	    pt, zt, tt, h2ot, lwc, iwc, cl, pct, pcb, plcl, plfc, pel, cape,
	    cin, RH(atm->p[ip], t, h2o), RHICE(atm->p[ip], t, h2o),
	    TDEW(atm->p[ip], h2o), TICE(atm->p[ip], h2o),
	    nat_temperature(atm->p[ip], h2o,
			    clim_hno3(atm->time[ip], atm->lat[ip],
				      atm->p[ip])), clim_hno3(atm->time[ip],
							      atm->lat[ip],
							      atm->p[ip]),
	    clim_oh(atm->time[ip], atm->lat[ip], atm->p[ip]), pbl);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(atm);
  free(met0);
  free(met1);

  return EXIT_SUCCESS;
}
