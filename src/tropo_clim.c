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
  Calculate tropopause climatology.
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

  static double lons[EX], lats[EY], *zm, *zs, *pm, *ps, *tm, *ts, *qm, *qs,
    *o3m, *o3s;

  static float *tropo_z0, *tropo_p0, *tropo_t0, *tropo_q0, *tropo_o30;

  static int ncid, varid, varid_z, varid_p, varid_t, varid_q, varid_o3, h2o,
    o3, *n, *nt, ntime, nlon, nlat;

  static size_t count[10], start[10];

  /* Allocate... */
  ALLOC(zm, double,
	EX * EY);
  ALLOC(zs, double,
	EX * EY);
  ALLOC(pm, double,
	EX * EY);
  ALLOC(ps, double,
	EX * EY);
  ALLOC(tm, double,
	EX * EY);
  ALLOC(ts, double,
	EX * EY);
  ALLOC(qm, double,
	EX * EY);
  ALLOC(qs, double,
	EX * EY);
  ALLOC(o3m, double,
	EX * EY);
  ALLOC(o3s, double,
	EX * EY);

  ALLOC(tropo_z0, float,
	EX * EY);
  ALLOC(tropo_p0, float,
	EX * EY);
  ALLOC(tropo_t0, float,
	EX * EY);
  ALLOC(tropo_q0, float,
	EX * EY);
  ALLOC(tropo_o30, float,
	EX * EY);

  ALLOC(n, int,
	EX * EY);
  ALLOC(nt, int,
	EX * EY);

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <clim.tab> <var> <tropo.nc>");

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

      /* Read data... */
      start[0] = (size_t) it;
      NC(nc_get_vara_float(ncid, varid_z, start, count, tropo_z0));
      NC(nc_get_vara_float(ncid, varid_p, start, count, tropo_p0));
      NC(nc_get_vara_float(ncid, varid_t, start, count, tropo_t0));
      if (h2o) {
	NC(nc_get_vara_float(ncid, varid_q, start, count, tropo_q0));
      } else
	for (int i = 0; i < nlon * nlat; i++)
	  tropo_q0[i] = NAN;
      if (o3) {
	NC(nc_get_vara_float(ncid, varid_o3, start, count, tropo_o30));
      } else
	for (int i = 0; i < nlon * nlat; i++)
	  tropo_o30[i] = NAN;

      /* Averaging... */
      for (int i = 0; i < nlon * nlat; i++) {
	nt[i]++;
	if (isfinite(tropo_z0[i])
	    && isfinite(tropo_p0[i])
	    && isfinite(tropo_t0[i])
	    && (!h2o || isfinite(tropo_q0[i]))
	    && (!o3 || isfinite(tropo_o30[i]))) {
	  zm[i] += tropo_z0[i];
	  zs[i] += SQR(tropo_z0[i]);
	  pm[i] += tropo_p0[i];
	  ps[i] += SQR(tropo_p0[i]);
	  tm[i] += tropo_t0[i];
	  ts[i] += SQR(tropo_t0[i]);
	  qm[i] += tropo_q0[i];
	  qs[i] += SQR(tropo_q0[i]);
	  o3m[i] += tropo_o30[i];
	  o3s[i] += SQR(tropo_o30[i]);
	  n[i]++;
	}
      }
    }

    /* Close files... */
    NC(nc_close(ncid));
  }

  /* Normalize... */
  for (int i = 0; i < nlon * nlat; i++)
    if (n[i] > 0) {
      zm[i] /= n[i];
      pm[i] /= n[i];
      tm[i] /= n[i];
      qm[i] /= n[i];
      o3m[i] /= n[i];
      double aux = zs[i] / n[i] - SQR(zm[i]);
      zs[i] = aux > 0 ? sqrt(aux) : 0.0;
      aux = ps[i] / n[i] - SQR(pm[i]);
      ps[i] = aux > 0 ? sqrt(aux) : 0.0;
      aux = ts[i] / n[i] - SQR(tm[i]);
      ts[i] = aux > 0 ? sqrt(aux) : 0.0;
      aux = qs[i] / n[i] - SQR(qm[i]);
      qs[i] = aux > 0 ? sqrt(aux) : 0.0;
      aux = o3s[i] / n[i] - SQR(o3m[i]);
      o3s[i] = aux > 0 ? sqrt(aux) : 0.0;
    }

  /* Create file... */
  LOG(1, "Write tropopause climatological data: %s", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = longitude [deg]\n"
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
	  "# $14 = occurrence frequency [%%]\n");

  /* Write output... */
  for (int ilat = 0; ilat < nlat; ilat++) {
    fprintf(out, "\n");
    for (int ilon = 0; ilon < nlon; ilon++)
      fprintf(out, "%g %g %g %g %g %g %g %g %g %g %g %g %d %g\n",
	      lons[ilon], lats[ilat], zm[ilat * nlon + ilon],
	      pm[ilat * nlon + ilon], tm[ilat * nlon + ilon],
	      qm[ilat * nlon + ilon], o3m[ilat * nlon + ilon],
	      zs[ilat * nlon + ilon], ps[ilat * nlon + ilon],
	      ts[ilat * nlon + ilon], qs[ilat * nlon + ilon],
	      o3s[ilat * nlon + ilon], n[ilat * nlon + ilon],
	      100. * n[ilat * nlon + ilon] / nt[ilat * nlon + ilon]);
  }

  /* Close files... */
  fclose(out);

  /* Free... */
  free(zm);
  free(zs);
  free(pm);
  free(ps);
  free(tm);
  free(ts);
  free(qm);
  free(qs);
  free(o3m);
  free(o3s);

  free(tropo_z0);
  free(tropo_p0);
  free(tropo_t0);
  free(tropo_q0);
  free(tropo_o30);

  free(n);
  free(nt);

  return EXIT_SUCCESS;
}
