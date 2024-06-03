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
  Sample tropopause data set.
*/

#include "mptrac.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Maximum number of time steps. */
#define NT 744

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  atm_t *atm;

  static FILE *out;

  static char varname[LEN];

  static double times[NT], lons[EX], lats[EY], time0, time1, z0, z0sig, p0,
    p0sig, t0, t0sig, q0, q0sig, o30, o30sig;

  static float help[EX * EY], tropo_z0[EX][EY], tropo_z1[EX][EY],
    tropo_p0[EX][EY], tropo_p1[EX][EY], tropo_t0[EX][EY], tropo_t1[EX][EY],
    tropo_q0[EX][EY], tropo_q1[EX][EY], tropo_o30[EX][EY], tropo_o31[EX][EY];

  static int it_old =
    -999, ncid, varid, varid_z, varid_p, varid_t, varid_q, varid_o3, h2o, o3,
    ntime, nlon, nlat, ilon, ilat;

  static size_t count[10], start[10];

  /* Allocate... */
  ALLOC(atm, atm_t, 1);

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <sample.tab> <tropo.nc> <var> <atm_in>");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  int method =
    (int) scan_ctl(argv[1], argc, argv, "TROPO_SAMPLE_METHOD", -1, "1", NULL);

  /* Read atmospheric data... */
  if (!read_atm(argv[5], &ctl, atm))
    ERRMSG("Cannot open file!");

  /* Open tropopause file... */
  LOG(1, "Read tropopause data: %s", argv[3]);
  if (nc_open(argv[3], NC_NOWRITE, &ncid) != NC_NOERR)
    ERRMSG("Cannot open file!");

  /* Get dimensions... */
  NC_INQ_DIM("time", &ntime, 1, NT);
  NC_INQ_DIM("lat", &nlat, 1, EY);
  NC_INQ_DIM("lon", &nlon, 1, EX);

  /* Read coordinates... */
  NC_GET_DOUBLE("time", times, 1);
  NC_GET_DOUBLE("lat", lats, 1);
  NC_GET_DOUBLE("lon", lons, 1);

  /* Get variable indices... */
  sprintf(varname, "%s_z", argv[4]);
  NC(nc_inq_varid(ncid, varname, &varid_z));
  sprintf(varname, "%s_p", argv[4]);
  NC(nc_inq_varid(ncid, varname, &varid_p));
  sprintf(varname, "%s_t", argv[4]);
  NC(nc_inq_varid(ncid, varname, &varid_t));
  sprintf(varname, "%s_q", argv[4]);
  h2o = (nc_inq_varid(ncid, varname, &varid_q) == NC_NOERR);
  sprintf(varname, "%s_o3", argv[4]);
  o3 = (nc_inq_varid(ncid, varname, &varid_o3) == NC_NOERR);

  /* Set dimensions... */
  count[0] = 1;
  count[1] = (size_t) nlat;
  count[2] = (size_t) nlon;

  /* Create file... */
  LOG(1, "Write tropopause sample data: %s", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = altitude [km]\n"
	  "# $3 = longitude [deg]\n" "# $4 = latitude [deg]\n");
  for (int iq = 0; iq < ctl.nq; iq++)
    fprintf(out, "# $%i = %s [%s]\n", iq + 5, ctl.qnt_name[iq],
	    ctl.qnt_unit[iq]);
  fprintf(out, "# $%d = tropopause height (mean) [km]\n", 5 + ctl.nq);
  fprintf(out, "# $%d = tropopause pressure (mean) [hPa]\n", 6 + ctl.nq);
  fprintf(out, "# $%d = tropopause temperature (mean) [K]\n", 7 + ctl.nq);
  fprintf(out, "# $%d = tropopause water vapor (mean) [ppv]\n", 8 + ctl.nq);
  fprintf(out, "# $%d = tropopause ozone (mean) [ppv]\n", 9 + ctl.nq);
  fprintf(out, "# $%d = tropopause height (sigma) [km]\n", 10 + ctl.nq);
  fprintf(out, "# $%d = tropopause pressure (sigma) [hPa]\n", 11 + ctl.nq);
  fprintf(out, "# $%d = tropopause temperature (sigma) [K]\n", 12 + ctl.nq);
  fprintf(out, "# $%d = tropopause water vapor (sigma) [ppv]\n", 13 + ctl.nq);
  fprintf(out, "# $%d = tropopause ozone (sigma) [ppv]\n\n", 14 + ctl.nq);

  /* Loop over particles... */
  for (int ip = 0; ip < atm->np; ip++) {

    /* Check temporal ordering... */
    if (ip > 0 && atm->time[ip] < atm->time[ip - 1])
      ERRMSG("Time must be ascending!");

    /* Check range... */
    if (atm->time[ip] < times[0] || atm->time[ip] > times[ntime - 1])
      continue;

    /* Read data... */
    int it = locate_irr(times, (int) ntime, atm->time[ip]);
    if (it != it_old) {

      time0 = times[it];
      start[0] = (size_t) it;
      NC(nc_get_vara_float(ncid, varid_z, start, count, help));
      for (ilon = 0; ilon < nlon; ilon++)
	for (ilat = 0; ilat < nlat; ilat++)
	  tropo_z0[ilon][ilat] = help[ilat * nlon + ilon];
      NC(nc_get_vara_float(ncid, varid_p, start, count, help));
      for (ilon = 0; ilon < nlon; ilon++)
	for (ilat = 0; ilat < nlat; ilat++)
	  tropo_p0[ilon][ilat] = help[ilat * nlon + ilon];
      NC(nc_get_vara_float(ncid, varid_t, start, count, help));
      for (ilon = 0; ilon < nlon; ilon++)
	for (ilat = 0; ilat < nlat; ilat++)
	  tropo_t0[ilon][ilat] = help[ilat * nlon + ilon];
      if (h2o) {
	NC(nc_get_vara_float(ncid, varid_q, start, count, help));
	for (ilon = 0; ilon < nlon; ilon++)
	  for (ilat = 0; ilat < nlat; ilat++)
	    tropo_q0[ilon][ilat] = help[ilat * nlon + ilon];
      } else
	for (ilon = 0; ilon < nlon; ilon++)
	  for (ilat = 0; ilat < nlat; ilat++)
	    tropo_q0[ilon][ilat] = NAN;
      if (o3) {
	NC(nc_get_vara_float(ncid, varid_o3, start, count, help));
	for (ilon = 0; ilon < nlon; ilon++)
	  for (ilat = 0; ilat < nlat; ilat++)
	    tropo_o30[ilon][ilat] = help[ilat * nlon + ilon];
      } else
	for (ilon = 0; ilon < nlon; ilon++)
	  for (ilat = 0; ilat < nlat; ilat++)
	    tropo_o30[ilon][ilat] = NAN;

      time1 = times[it + 1];
      start[0] = (size_t) it + 1;
      NC(nc_get_vara_float(ncid, varid_z, start, count, help));
      for (ilon = 0; ilon < nlon; ilon++)
	for (ilat = 0; ilat < nlat; ilat++)
	  tropo_z1[ilon][ilat] = help[ilat * nlon + ilon];
      NC(nc_get_vara_float(ncid, varid_p, start, count, help));
      for (ilon = 0; ilon < nlon; ilon++)
	for (ilat = 0; ilat < nlat; ilat++)
	  tropo_p1[ilon][ilat] = help[ilat * nlon + ilon];
      NC(nc_get_vara_float(ncid, varid_t, start, count, help));
      for (ilon = 0; ilon < nlon; ilon++)
	for (ilat = 0; ilat < nlat; ilat++)
	  tropo_t1[ilon][ilat] = help[ilat * nlon + ilon];
      if (h2o) {
	NC(nc_get_vara_float(ncid, varid_q, start, count, help));
	for (ilon = 0; ilon < nlon; ilon++)
	  for (ilat = 0; ilat < nlat; ilat++)
	    tropo_q1[ilon][ilat] = help[ilat * nlon + ilon];
      } else
	for (ilon = 0; ilon < nlon; ilon++)
	  for (ilat = 0; ilat < nlat; ilat++)
	    tropo_q1[ilon][ilat] = NAN;
      if (o3) {
	NC(nc_get_vara_float(ncid, varid_o3, start, count, help));
	for (ilon = 0; ilon < nlon; ilon++)
	  for (ilat = 0; ilat < nlat; ilat++)
	    tropo_o31[ilon][ilat] = help[ilat * nlon + ilon];
      } else
	for (ilon = 0; ilon < nlon; ilon++)
	  for (ilat = 0; ilat < nlat; ilat++)
	    tropo_o31[ilon][ilat] = NAN;
    }
    it_old = it;

    /* Interpolate... */
    intpol_tropo_3d(time0, tropo_z0, time1, tropo_z1,
		    lons, lats, nlon, nlat, atm->time[ip], atm->lon[ip],
		    atm->lat[ip], method, &z0, &z0sig);
    intpol_tropo_3d(time0, tropo_p0, time1, tropo_p1,
		    lons, lats, nlon, nlat, atm->time[ip], atm->lon[ip],
		    atm->lat[ip], method, &p0, &p0sig);
    intpol_tropo_3d(time0, tropo_t0, time1, tropo_t1,
		    lons, lats, nlon, nlat, atm->time[ip], atm->lon[ip],
		    atm->lat[ip], method, &t0, &t0sig);
    intpol_tropo_3d(time0, tropo_q0, time1, tropo_q1,
		    lons, lats, nlon, nlat, atm->time[ip], atm->lon[ip],
		    atm->lat[ip], method, &q0, &q0sig);
    intpol_tropo_3d(time0, tropo_o30, time1, tropo_o31,
		    lons, lats, nlon, nlat, atm->time[ip], atm->lon[ip],
		    atm->lat[ip], method, &o30, &o30sig);

    /* Write output... */
    fprintf(out, "%.2f %g %g %g", atm->time[ip], Z(atm->p[ip]),
	    atm->lon[ip], atm->lat[ip]);
    for (int iq = 0; iq < ctl.nq; iq++) {
      fprintf(out, " ");
      fprintf(out, ctl.qnt_format[iq], atm->q[iq][ip]);
    }
    fprintf(out, " %g %g %g %g %g %g %g %g %g %g\n",
	    z0, p0, t0, q0, o30, z0sig, p0sig, t0sig, q0sig, o30sig);
  }

  /* Close files... */
  fclose(out);
  NC(nc_close(ncid));

  /* Free... */
  free(atm);

  return EXIT_SUCCESS;
}
