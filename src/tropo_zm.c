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
  
  Copyright (C) 2013-2025 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Extract zonal mean of tropopause data set.
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

  static FILE *out;

  static char varname[LEN];

  static double time0, lons[EX], lats[EY], zm[EY], zs[EY], pm[EY],
    ps[EY], tm[EY], ts[EY], qm[EY], qs[EY], o3m[EY], o3s[EY];

  static float help[EX * EY], tropo_z0[EX][EY], tropo_p0[EX][EY],
    tropo_t0[EX][EY], tropo_q0[EX][EY], tropo_o30[EX][EY];

  static int ncid, varid, varid_z, varid_p, varid_t, varid_q, varid_o3, h2o,
    o3, n[EY], nt[EY], init, ntime, nlon, nlat, ilon, ilat;

  static size_t count[10], start[10];

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <zm.tab> <var> <tropo.nc>");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);

  /* Loop over tropopause files... */
  for (int iarg = 4; iarg < argc; iarg++) {

    /* Open tropopause file... */
    LOG(1, "Read tropopause data: %s", argv[iarg]);
    if (nc_open(argv[iarg], NC_NOWRITE, &ncid) != NC_NOERR)
      ERRMSG("Cannot open file!");

    /* Get dimensions... */
    NC_INQ_DIM("time", &ntime, 1, NT);
    NC_INQ_DIM("lat", &nlat, 1, EY);
    NC_INQ_DIM("lon", &nlon, 1, EX);

    /* Read coordinates... */
    NC_GET_DOUBLE("lat", lats, 1);
    NC_GET_DOUBLE("lon", lons, 1);

    /* Get variable indices... */
    sprintf(varname, "%s_z", argv[3]);
    NC(nc_inq_varid(ncid, varname, &varid_z));
    sprintf(varname, "%s_p", argv[3]);
    NC(nc_inq_varid(ncid, varname, &varid_p));
    sprintf(varname, "%s_t", argv[3]);
    NC(nc_inq_varid(ncid, varname, &varid_t));
    sprintf(varname, "%s_q", argv[3]);
    h2o = (nc_inq_varid(ncid, varname, &varid_q) == NC_NOERR);
    sprintf(varname, "%s_o3", argv[3]);
    o3 = (nc_inq_varid(ncid, varname, &varid_o3) == NC_NOERR);

    /* Set dimensions... */
    count[0] = 1;
    count[1] = (size_t) nlat;
    count[2] = (size_t) nlon;

    /* Loop over time steps... */
    for (int it = 0; it < ntime; it++) {

      /* Get time from filename... */
      if (!init) {
	init = 1;
	time0 = time_from_filename(argv[iarg], 13);
      }

      /* Read data... */
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

      /* Averaging... */
      for (ilat = 0; ilat < nlat; ilat++)
	for (ilon = 0; ilon < nlon; ilon++) {
	  nt[ilat]++;
	  if (isfinite(tropo_z0[ilon][ilat])
	      && isfinite(tropo_p0[ilon][ilat])
	      && isfinite(tropo_t0[ilon][ilat])
	      && (!h2o || isfinite(tropo_q0[ilon][ilat]))
	      && (!o3 || isfinite(tropo_o30[ilon][ilat]))) {
	    zm[ilat] += tropo_z0[ilon][ilat];
	    zs[ilat] += SQR(tropo_z0[ilon][ilat]);
	    pm[ilat] += tropo_p0[ilon][ilat];
	    ps[ilat] += SQR(tropo_p0[ilon][ilat]);
	    tm[ilat] += tropo_t0[ilon][ilat];
	    ts[ilat] += SQR(tropo_t0[ilon][ilat]);
	    qm[ilat] += tropo_q0[ilon][ilat];
	    qs[ilat] += SQR(tropo_q0[ilon][ilat]);
	    o3m[ilat] += tropo_o30[ilon][ilat];
	    o3s[ilat] += SQR(tropo_o30[ilon][ilat]);
	    n[ilat]++;
	  }
	}
    }

    /* Close files... */
    NC(nc_close(ncid));
  }

  /* Normalize... */
  for (ilat = 0; ilat < nlat; ilat++)
    if (n[ilat] > 0) {
      zm[ilat] /= n[ilat];
      pm[ilat] /= n[ilat];
      tm[ilat] /= n[ilat];
      qm[ilat] /= n[ilat];
      o3m[ilat] /= n[ilat];
      double aux = zs[ilat] / n[ilat] - SQR(zm[ilat]);
      zs[ilat] = aux > 0 ? sqrt(aux) : 0.0;
      aux = ps[ilat] / n[ilat] - SQR(pm[ilat]);
      ps[ilat] = aux > 0 ? sqrt(aux) : 0.0;
      aux = ts[ilat] / n[ilat] - SQR(tm[ilat]);
      ts[ilat] = aux > 0 ? sqrt(aux) : 0.0;
      aux = qs[ilat] / n[ilat] - SQR(qm[ilat]);
      qs[ilat] = aux > 0 ? sqrt(aux) : 0.0;
      aux = o3s[ilat] / n[ilat] - SQR(o3m[ilat]);
      o3s[ilat] = aux > 0 ? sqrt(aux) : 0.0;
    }

  /* Create file... */
  LOG(1, "Write tropopause zonal mean data: %s", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = latitude [deg]\n"
	  "# $3 = tropopause height (mean) [km]\n"
	  "# $4 = tropopause pressure (mean) [hPa]\n"
	  "# $5 = tropopause temperature (mean) [K]\n"
	  "# $6 = tropopause water vapor (mean) [ppv]\n"
	  "# $7 = tropopause ozone (mean) [ppv]\n"
	  "# $8 = tropopause height (sigma) [km]\n"
	  "# $9 = tropopause pressure (sigma) [hPa]\n"
	  "# $10 = tropopause temperature (sigma) [K]\n"
	  "# $11 = tropopause water vapor (sigma) [ppv]\n"
	  "# $12 = tropopause ozone (sigma) [ppv]\n"
	  "# $13 = number of data points\n"
	  "# $14 = occurrence frequency [%%]\n\n");

  /* Write output... */
  for (ilat = 0; ilat < nlat; ilat++)
    fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %g %g %d %g\n",
	    time0, lats[ilat], zm[ilat], pm[ilat], tm[ilat], qm[ilat],
	    o3m[ilat], zs[ilat], ps[ilat], ts[ilat], qs[ilat], o3s[ilat],
	    n[ilat], 100. * n[ilat] / nt[ilat]);

  /* Close files... */
  fclose(out);

  return EXIT_SUCCESS;
}
