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
  Sample tropopause data set.
*/

#include "libtrac.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Maximum number of time steps. */
#define NT 744

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! 3-D linear interpolation of tropopause data. */
double intpol_tropo_3d(
  double time0,
  float array0[EX][EY],
  float arrayz0[EX][EY],
  double time1,
  float array1[EX][EY],
  float arrayz1[EX][EY],
  double lons[EX],
  double lats[EY],
  size_t nlon,
  size_t nlat,
  double time,
  double lon,
  double lat,
  double deltaz);

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

  static double times[NT], lons[EX], lats[EY], time0, time1, z0, p0, t0, q0,
    deltaz;

  static float help[EX * EY], tropo_z0[EX][EY], tropo_z1[EX][EY],
    tropo_p0[EX][EY], tropo_p1[EX][EY], tropo_t0[EX][EY],
    tropo_t1[EX][EY], tropo_q0[EX][EY], tropo_q1[EX][EY];

  static int ip, iq, it, it_old = -999, dimid[10], ncid,
    varid, varid_z, varid_p, varid_t, varid_q, h2o;

  static size_t count[10], start[10], ntime, nlon, nlat, ilon, ilat;

  /* Allocate... */
  ALLOC(atm, atm_t, 1);

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <sample.tab> <tropo.nc> <var> <atm_in>");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  deltaz =
    scan_ctl(argv[1], argc, argv, "TROPO_SAMPLE_DELTAZ", -1, "-999", NULL);

  /* Read atmospheric data... */
  if (!read_atm(argv[5], &ctl, atm))
    ERRMSG("Cannot open file!");

  /* Open tropopause file... */
  printf("Read tropopause data: %s\n", argv[3]);
  if (nc_open(argv[3], NC_NOWRITE, &ncid) != NC_NOERR)
    ERRMSG("Cannot open file!");

  /* Get dimensions... */
  NC(nc_inq_dimid(ncid, "time", &dimid[0]));
  NC(nc_inq_dimlen(ncid, dimid[0], &ntime));
  if (ntime > NT)
    ERRMSG("Too many times!");
  NC(nc_inq_dimid(ncid, "lat", &dimid[1]));
  NC(nc_inq_dimlen(ncid, dimid[1], &nlat));
  if (nlat > EY)
    ERRMSG("Too many latitudes!");
  NC(nc_inq_dimid(ncid, "lon", &dimid[2]));
  NC(nc_inq_dimlen(ncid, dimid[2], &nlon));
  if (nlon > EX)
    ERRMSG("Too many longitudes!");

  /* Read coordinates... */
  NC(nc_inq_varid(ncid, "time", &varid));
  NC(nc_get_var_double(ncid, varid, times));
  NC(nc_inq_varid(ncid, "lat", &varid));
  NC(nc_get_var_double(ncid, varid, lats));
  NC(nc_inq_varid(ncid, "lon", &varid));
  NC(nc_get_var_double(ncid, varid, lons));

  /* Get variable indices... */
  sprintf(varname, "%s_z", argv[4]);
  NC(nc_inq_varid(ncid, varname, &varid_z));
  sprintf(varname, "%s_p", argv[4]);
  NC(nc_inq_varid(ncid, varname, &varid_p));
  sprintf(varname, "%s_t", argv[4]);
  NC(nc_inq_varid(ncid, varname, &varid_t));
  sprintf(varname, "%s_q", argv[4]);
  h2o = (nc_inq_varid(ncid, varname, &varid_q) == NC_NOERR);

  /* Set dimensions... */
  count[0] = 1;
  count[1] = nlat;
  count[2] = nlon;

  /* Create file... */
  printf("Write tropopause sample data: %s\n", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = altitude [km]\n"
	  "# $3 = longitude [deg]\n" "# $4 = latitude [deg]\n");
  for (iq = 0; iq < ctl.nq; iq++)
    fprintf(out, "# $%i = %s [%s]\n", iq + 5, ctl.qnt_name[iq],
	    ctl.qnt_unit[iq]);
  fprintf(out, "# $%d = tropopause height [km]\n", 5 + ctl.nq);
  fprintf(out, "# $%d = tropopause pressure [hPa]\n", 6 + ctl.nq);
  fprintf(out, "# $%d = tropopause temperature [K]\n", 7 + ctl.nq);
  fprintf(out, "# $%d = tropopause water vapor [ppv]\n\n", 8 + ctl.nq);

  /* Loop over particles... */
  for (ip = 0; ip < atm->np; ip++) {

    /* Check temporal ordering... */
    if (ip > 0 && atm->time[ip] < atm->time[ip - 1])
      ERRMSG("Time must be ascending!");

    /* Check range... */
    if (atm->time[ip] < times[0] || atm->time[ip] > times[ntime - 1])
      continue;

    /* Read data... */
    it = locate_irr(times, (int) ntime, atm->time[ip]);
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
	    tropo_q0[ilon][ilat] = GSL_NAN;

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
	    tropo_q1[ilon][ilat] = GSL_NAN;;
    }
    it_old = it;

    /* Interpolate... */
    z0 =
      intpol_tropo_3d(time0, tropo_z0, tropo_z0, time1, tropo_z1, tropo_z1,
		      lons, lats, nlon, nlat, atm->time[ip], atm->lon[ip],
		      atm->lat[ip], deltaz);
    p0 =
      intpol_tropo_3d(time0, tropo_p0, tropo_z0, time1, tropo_p1, tropo_z1,
		      lons, lats, nlon, nlat, atm->time[ip], atm->lon[ip],
		      atm->lat[ip], deltaz);
    t0 =
      intpol_tropo_3d(time0, tropo_t0, tropo_z0, time1, tropo_t1, tropo_z1,
		      lons, lats, nlon, nlat, atm->time[ip], atm->lon[ip],
		      atm->lat[ip], deltaz);
    q0 =
      intpol_tropo_3d(time0, tropo_q0, tropo_z0, time1, tropo_q1, tropo_z1,
		      lons, lats, nlon, nlat, atm->time[ip], atm->lon[ip],
		      atm->lat[ip], deltaz);

    /* Write output... */
    fprintf(out, "%.2f %g %g %g", atm->time[ip], Z(atm->p[ip]),
	    atm->lon[ip], atm->lat[ip]);
    for (iq = 0; iq < ctl.nq; iq++) {
      fprintf(out, " ");
      fprintf(out, ctl.qnt_format[iq], atm->q[iq][ip]);
    }
    fprintf(out, " %g %g %g %g\n", z0, p0, t0, q0);
  }

  /* Close files... */
  fclose(out);
  NC(nc_close(ncid));

  /* Free... */
  free(atm);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

double intpol_tropo_3d(
  double time0,
  float array0[EX][EY],
  float arrayz0[EX][EY],
  double time1,
  float array1[EX][EY],
  float arrayz1[EX][EY],
  double lons[EX],
  double lats[EY],
  size_t nlon,
  size_t nlat,
  double time,
  double lon,
  double lat,
  double deltaz) {

  double sum, w, wsum;

  /* Adjust longitude... */
  if (lon < lons[0])
    lon += 360;
  else if (lon > lons[nlon - 1])
    lon -= 360;

  /* Get indices... */
  int ix = locate_reg(lons, (int) nlon, lon);
  int iy = locate_reg(lats, (int) nlat, lat);

  /* Check gradients... */
  if (deltaz > 0)
    if (fabs(arrayz0[ix + 1][iy] - arrayz0[ix][iy]) > deltaz
	|| fabs(arrayz0[ix + 1][iy + 1] - arrayz0[ix][iy + 1]) > deltaz
	|| fabs(arrayz0[ix][iy + 1] - arrayz0[ix][iy]) > deltaz
	|| fabs(arrayz0[ix + 1][iy + 1] - arrayz0[ix + 1][iy]) > deltaz
	|| fabs(arrayz1[ix + 1][iy] - arrayz1[ix][iy]) > deltaz
	|| fabs(arrayz1[ix + 1][iy + 1] - arrayz1[ix][iy + 1]) > deltaz
	|| fabs(arrayz1[ix][iy + 1] - arrayz1[ix][iy]) > deltaz
	|| fabs(arrayz1[ix + 1][iy + 1] - arrayz1[ix + 1][iy]) > deltaz)
      return GSL_NAN;

  /* Init... */
  sum = 0;
  wsum = 0;

  /* Interpolate... */
  if (isfinite(array0[ix][iy])) {
    w = (lons[ix + 1] - lon) * (lats[iy + 1] - lat) * (time1 - time);
    sum += w * array0[ix][iy];
    wsum += w;
  } else {
    if (fabs(lon - lons[ix]) < fabs(lon - lons[ix + 1])
	&& fabs(lat - lats[iy]) < fabs(lat - lats[iy + 1])
	&& fabs(time - time0) < fabs(time - time1))
      sum += GSL_NAN;
  }

  if (isfinite(array0[ix + 1][iy])) {
    w = (lon - lons[ix]) * (lats[iy + 1] - lat) * (time1 - time);
    sum += w * array0[ix + 1][iy];
    wsum += w;
  } else {
    if (fabs(lon - lons[ix + 1]) < fabs(lon - lons[ix])
	&& fabs(lat - lats[iy]) < fabs(lat - lats[iy + 1])
	&& fabs(time - time0) < fabs(time - time1))
      sum += GSL_NAN;
  }

  if (isfinite(array0[ix][iy + 1])) {
    w = (lons[ix + 1] - lon) * (lat - lats[iy]) * (time1 - time);
    sum += w * array0[ix][iy + 1];
    wsum += w;
  } else {
    if (fabs(lon - lons[ix]) < fabs(lon - lons[ix + 1])
	&& fabs(lat - lats[iy + 1]) < fabs(lat - lats[iy])
	&& fabs(time - time0) < fabs(time - time1))
      sum += GSL_NAN;
  }

  if (isfinite(array0[ix + 1][iy + 1])) {
    w = (lon - lons[ix]) * (lat - lats[iy]) * (time1 - time);
    sum += w * array0[ix + 1][iy + 1];
    wsum += w;
  } else {
    if (fabs(lon - lons[ix + 1]) < fabs(lon - lons[ix])
	&& fabs(lat - lats[iy + 1]) < fabs(lat - lats[iy])
	&& fabs(time - time0) < fabs(time - time1))
      sum += GSL_NAN;
  }

  if (isfinite(array1[ix][iy])) {
    w = (lons[ix + 1] - lon) * (lats[iy + 1] - lat) * (time - time0);
    sum += w * array1[ix][iy];
    wsum += w;
  } else {
    if (fabs(lon - lons[ix]) < fabs(lon - lons[ix + 1])
	&& fabs(lat - lats[iy]) < fabs(lat - lats[iy + 1])
	&& fabs(time - time1) < fabs(time - time0))
      sum += GSL_NAN;
  }

  if (isfinite(array1[ix + 1][iy])) {
    w = (lon - lons[ix]) * (lats[iy + 1] - lat) * (time - time0);
    sum += w * array1[ix + 1][iy];
    wsum += w;
  } else {
    if (fabs(lon - lons[ix + 1]) < fabs(lon - lons[ix])
	&& fabs(lat - lats[iy]) < fabs(lat - lats[iy + 1])
	&& fabs(time - time1) < fabs(time - time0))
      sum += GSL_NAN;
  }

  if (isfinite(array1[ix][iy + 1])) {
    w = (lons[ix + 1] - lon) * (lat - lats[iy]) * (time - time0);
    sum += w * array1[ix][iy + 1];
    wsum += w;
  } else {
    if (fabs(lon - lons[ix]) < fabs(lon - lons[ix + 1])
	&& fabs(lat - lats[iy + 1]) < fabs(lat - lats[iy])
	&& fabs(time - time1) < fabs(time - time0))
      sum += GSL_NAN;
  }

  if (isfinite(array1[ix + 1][iy + 1])) {
    w = (lon - lons[ix]) * (lat - lats[iy]) * (time - time0);
    sum += w * array1[ix + 1][iy + 1];
    wsum += w;
  } else {
    if (fabs(lon - lons[ix + 1]) < fabs(lon - lons[ix])
	&& fabs(lat - lats[iy + 1]) < fabs(lat - lats[iy])
	&& fabs(time - time1) < fabs(time - time0))
      sum += GSL_NAN;
  }

  /* Return value... */
  if (wsum > 0)
    return sum / wsum;
  else
    return GSL_NAN;
}
