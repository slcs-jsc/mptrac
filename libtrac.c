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
  MPTRAC library definitions.
*/

#include "libtrac.h"

/*****************************************************************************/

void cart2geo(
  double *x,
  double *z,
  double *lon,
  double *lat) {

  double radius;

  radius = NORM(x);
  *lat = asin(x[2] / radius) * 180 / M_PI;
  *lon = atan2(x[1], x[0]) * 180 / M_PI;
  *z = radius - RE;
}

/*****************************************************************************/

double deg2dx(
  double dlon,
  double lat) {

  return dlon * M_PI * RE / 180. * cos(lat / 180. * M_PI);
}

/*****************************************************************************/

double deg2dy(
  double dlat) {

  return dlat * M_PI * RE / 180.;
}

/*****************************************************************************/

double dp2dz(
  double dp,
  double p) {

  return -dp * H0 / p;
}

/*****************************************************************************/

double dx2deg(
  double dx,
  double lat) {

  /* Avoid singularity at poles... */
  if (lat < -89.999 || lat > 89.999)
    return 0;
  else
    return dx * 180. / (M_PI * RE * cos(lat / 180. * M_PI));
}

/*****************************************************************************/

double dy2deg(
  double dy) {

  return dy * 180. / (M_PI * RE);
}

/*****************************************************************************/

double dz2dp(
  double dz,
  double p) {

  return -dz * p / H0;
}

/*****************************************************************************/

void extrapolate_met(
  met_t * met) {

  int ip, ip0, ix, iy;

  /* Loop over columns... */
  for (ix = 0; ix < met->nx; ix++)
    for (iy = 0; iy < met->ny; iy++) {

      /* Find lowest valid data point... */
      for (ip0 = met->np - 1; ip0 >= 0; ip0--)
	if (!gsl_finite(met->u[ix][iy][ip0])
	    || !gsl_finite(met->v[ix][iy][ip0])
	    || !gsl_finite(met->w[ix][iy][ip0])
	    || !gsl_finite(met->t[ix][iy][ip0]))
	  break;

      /* Extrapolate... */
      for (ip = ip0; ip >= 0; ip--) {
	met->t[ix][iy][ip]
	  = (float) LIN(met->p[ip + 1], met->t[ix][iy][ip + 1],
			met->p[ip + 2], met->t[ix][iy][ip + 2], met->p[ip]);
	met->u[ix][iy][ip]
	  = (float) LIN(met->p[ip + 1], met->u[ix][iy][ip + 1],
			met->p[ip + 2], met->u[ix][iy][ip + 2], met->p[ip]);
	met->v[ix][iy][ip]
	  = (float) LIN(met->p[ip + 1], met->v[ix][iy][ip + 1],
			met->p[ip + 2], met->v[ix][iy][ip + 2], met->p[ip]);
	met->w[ix][iy][ip]
	  = (float) LIN(met->p[ip + 1], met->w[ix][iy][ip + 1],
			met->p[ip + 2], met->w[ix][iy][ip + 2], met->p[ip]);
      }
    }
}

/*****************************************************************************/

void geo2cart(
  double z,
  double lon,
  double lat,
  double *x) {

  double radius;

  radius = z + RE;
  x[0] = radius * cos(lat / 180 * M_PI) * cos(lon / 180 * M_PI);
  x[1] = radius * cos(lat / 180 * M_PI) * sin(lon / 180 * M_PI);
  x[2] = radius * sin(lat / 180 * M_PI);
}

/*****************************************************************************/

void get_met(
  double t,
  int direct,
  char *metbase,
  double dt_met,
  int reduce,
  met_t * met0,
  met_t * met1) {

  char filename[LEN];

  static int init;

  /* Initialize at first call... */
  if (!init) {
    init = 1;

    get_met_help(t, -1, metbase, dt_met, filename);
    read_met(filename, met0);
    extrapolate_met(met0);
    reduce_met(met0, reduce, reduce, 1);

    get_met_help(t, 1, metbase, dt_met, filename);
    read_met(filename, met1);
    extrapolate_met(met1);
    reduce_met(met1, reduce, reduce, 1);
  }

  /* Read new data for forward trajectories... */
  if (t > met1->time && direct == 1) {
    memcpy(met0, met1, sizeof(met_t));
    get_met_help(t, 1, metbase, dt_met, filename);
    read_met(filename, met1);
    extrapolate_met(met1);
    reduce_met(met1, reduce, reduce, 1);
  }

  /* Read new data for backward trajectories... */
  if (t < met0->time && direct == -1) {
    memcpy(met1, met0, sizeof(met_t));
    get_met_help(t, -1, metbase, dt_met, filename);
    read_met(filename, met0);
    extrapolate_met(met0);
    reduce_met(met0, reduce, reduce, 1);
  }
}

/*****************************************************************************/

void get_met_help(
  double t,
  int direct,
  char *metbase,
  double dt_met,
  char *filename) {

  double t6, r;

  int year, mon, day, hour, min, sec;

  /* Round time to fixed intervals... */
  if (direct == -1)
    t6 = floor(t / dt_met) * dt_met;
  else
    t6 = ceil(t / dt_met) * dt_met;

  /* Decode time... */
  jsec2time(t6, &year, &mon, &day, &hour, &min, &sec, &r);

  /* Set filename... */
  sprintf(filename, "%s_%d_%02d_%02d_%02d.nc", metbase, year, mon, day, hour);
}

/*****************************************************************************/

void intpol_met_help(
  float array[EX][EY][EP],
  int ip,
  int ix,
  int iy,
  double wp,
  double wx,
  double wy,
  double *var) {

  double aux00, aux01, aux10, aux11;

  /* Interpolate vertically... */
  aux00 = wp * (array[ix][iy][ip] - array[ix][iy][ip + 1])
    + array[ix][iy][ip + 1];
  aux01 = wp * (array[ix][iy + 1][ip] - array[ix][iy + 1][ip + 1])
    + array[ix][iy + 1][ip + 1];
  aux10 = wp * (array[ix + 1][iy][ip] - array[ix + 1][iy][ip + 1])
    + array[ix + 1][iy][ip + 1];
  aux11 = wp * (array[ix + 1][iy + 1][ip] - array[ix + 1][iy + 1][ip + 1])
    + array[ix + 1][iy + 1][ip + 1];

  /* Interpolate horizontally... */
  aux00 = wy * (aux00 - aux01) + aux01;
  aux11 = wy * (aux10 - aux11) + aux11;
  *var = wx * (aux00 - aux11) + aux11;
}

/*****************************************************************************/

void intpol_met_space(
  met_t * met,
  double p,
  double lon,
  double lat,
  double *t,
  double *u,
  double *v,
  double *w) {

  double wp, wx, wy;

  int ip, ix, iy;

  /* Check longitude... */
  if (lon < 0)
    lon += 360;

  /* Get indices... */
  ip = locate(met->p, met->np, p);
  ix = locate(met->lon, met->nx, lon);
  iy = locate(met->lat, met->ny, lat);

  /* Get weights... */
  wp = (met->p[ip + 1] - p) / (met->p[ip + 1] - met->p[ip]);
  wx = (met->lon[ix + 1] - lon) / (met->lon[ix + 1] - met->lon[ix]);
  wy = (met->lat[iy + 1] - lat) / (met->lat[iy + 1] - met->lat[iy]);

  /* Interpolate... */
  intpol_met_help(met->t, ip, ix, iy, wp, wx, wy, t);
  intpol_met_help(met->u, ip, ix, iy, wp, wx, wy, u);
  intpol_met_help(met->v, ip, ix, iy, wp, wx, wy, v);
  intpol_met_help(met->w, ip, ix, iy, wp, wx, wy, w);
}

/*****************************************************************************/

void intpol_met_time(
  met_t * met0,
  met_t * met1,
  double ts,
  double p,
  double lon,
  double lat,
  double *t,
  double *u,
  double *v,
  double *w) {

  double t0, t1, u0, u1, v0, v1, w0, w1, wt;

  /* Spatial interpolation... */
  intpol_met_space(met0, p, lon, lat, &t0, &u0, &v0, &w0);
  intpol_met_space(met1, p, lon, lat, &t1, &u1, &v1, &w1);

  /* Get weighting factor... */
  wt = (met1->time - ts) / (met1->time - met0->time);

  /* Interpolate... */
  *t = wt * (t0 - t1) + t1;
  *u = wt * (u0 - u1) + u1;
  *v = wt * (v0 - v1) + v1;
  *w = wt * (w0 - w1) + w1;
}

/*****************************************************************************/

void jsec2time(
  double jsec,
  int *year,
  int *mon,
  int *day,
  int *hour,
  int *min,
  int *sec,
  double *remain) {

  struct tm t0, *t1;

  time_t jsec0;

  t0.tm_year = 100;
  t0.tm_mon = 0;
  t0.tm_mday = 1;
  t0.tm_hour = 0;
  t0.tm_min = 0;
  t0.tm_sec = 0;

  jsec0 = (time_t) jsec + timegm(&t0);
  t1 = gmtime(&jsec0);

  *year = t1->tm_year + 1900;
  *mon = t1->tm_mon + 1;
  *day = t1->tm_mday;
  *hour = t1->tm_hour;
  *min = t1->tm_min;
  *sec = t1->tm_sec;
  *remain = jsec - floor(jsec);
}

/*****************************************************************************/

int locate(
  double *xx,
  int n,
  double x) {

  int i, ilo, ihi;

  ilo = 0;
  ihi = n - 1;
  i = (ihi + ilo) >> 1;

  if (xx[i] < xx[i + 1])
    while (ihi > ilo + 1) {
      i = (ihi + ilo) >> 1;
      if (xx[i] > x)
	ihi = i;
      else
	ilo = i;
  } else
    while (ihi > ilo + 1) {
      i = (ihi + ilo) >> 1;
      if (xx[i] <= x)
	ihi = i;
      else
	ilo = i;
    }

  return ilo;
}

/*****************************************************************************/

void read_atm(
  const char *dirname,
  const char *filename,
  atm_t * atm,
  ctl_t * ctl) {

  char file[LEN];

  /* Set filename... */
  if (dirname != NULL)
    sprintf(file, "%s/%s", dirname, filename);
  else
    sprintf(file, "%s", filename);

  /* Write info... */
  printf("Read atmospheric data: %s\n", file);

  /* Call function for the selected file format. */
  if (ctl->atm_iformat == 0)
    read_atm_from_ascii(file, atm, ctl);
  else if (ctl->atm_iformat == 1)
    read_atm_from_netcdf(file, atm, ctl);
  else
    ERRMSG("File format undefined!");
}


/*****************************************************************************/

void read_atm_from_ascii(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl) {

  FILE *in;

  char line[LEN], *tok;

  int iq;

  atm->np = 0;

  /* Open file... */
  if (!(in = fopen(filename, "r"))) {
    printf("%s\n", filename);
    ERRMSG("Cannot open file!");
  }

  /* Read line... */
  while (fgets(line, LEN, in)) {

    /* Read data... */
    TOK(line, tok, "%lg", atm->time[atm->np]);
    TOK(NULL, tok, "%lg", atm->p[atm->np]);
    TOK(NULL, tok, "%lg", atm->lon[atm->np]);
    TOK(NULL, tok, "%lg", atm->lat[atm->np]);
    for (iq = 0; iq < ctl->nq; iq++)
      TOK(NULL, tok, "%lg", atm->q[iq][atm->np]);

    /* Convert altitude to pressure... */
    atm->p[atm->np] = P(atm->p[atm->np]);

    /* Increment data point counter... */
    if ((++atm->np) > NP)
      ERRMSG("Too many data points!");
  }

  /* Close file... */
  fclose(in);

  /* Check number of points... */
  if (atm->np < 1)
    ERRMSG("Could not read any data!");
}

/*****************************************************************************/

void read_atm_from_netcdf(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl) {

  size_t size, loadnum[1], startpos[1];

  int iq, varid, ncid, dimid;

  /* Open file... */
  NC(nc_open(filename, NC_NOWRITE, &ncid));

  /* Read dimensions... */
  NC(nc_inq_dimid(ncid, "np", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &size));
  atm->np = (int) size;

  if (atm->np > NP)
    ERRMSG("Too many data points!");
  startpos[0] = (size_t) 0;
  loadnum[0] = (size_t) atm->np;

  /* Read coordinate variables... */
  NC(nc_inq_varid(ncid, "time", &varid));
  NC(nc_get_vara_double(ncid, varid, startpos, loadnum, atm->time));
  NC(nc_inq_varid(ncid, "pressure", &varid));
  NC(nc_get_vara_double(ncid, varid, startpos, loadnum, atm->p));
  NC(nc_inq_varid(ncid, "longitude", &varid));
  NC(nc_get_vara_double(ncid, varid, startpos, loadnum, atm->lon));
  NC(nc_inq_varid(ncid, "latitude", &varid));
  NC(nc_get_vara_double(ncid, varid, startpos, loadnum, atm->lat));

  /* Read quantities... */
  for (iq = 0; iq < ctl->nq; iq++) {
    NC(nc_inq_varid(ncid, ctl->qnt_name[iq], &varid));
    NC(nc_get_vara_double(ncid, varid, startpos, loadnum, atm->q[iq]));
  }

  /* Close file... */
  NC(nc_close(ncid));
}

/*****************************************************************************/

void read_ctl(
  const char *dirname,
  const char *filename,
  int argc,
  char *argv[],
  ctl_t * ctl) {

  int iq;

  /* Write info... */
  printf("\nMassive-Parallel Trajectory Calculations (MPTRAC)\n"
	 "(executable: %s | compiled: %s, %s)\n\n",
	 argv[0], __DATE__, __TIME__);

  /* Initialize quantity indices... */
  ctl->qnt_mass = -1;
  ctl->qnt_r_p = -1;
  ctl->qnt_rho_p = -1;
  ctl->qnt_station = -1;
  ctl->qnt_temp = -1;

  /* Read quantities... */
  ctl->nq =
    (int) scan_ctl(dirname, filename, argc, argv, "NQ", -1, "0", NULL);
  for (iq = 0; iq < ctl->nq; iq++) {

    /* Read quantity names, units, and formats... */
    scan_ctl(dirname, filename, argc, argv, "QNT_NAME", iq, "",
	     ctl->qnt_name[iq]);
    scan_ctl(dirname, filename, argc, argv, "QNT_UNIT", iq, "",
	     ctl->qnt_unit[iq]);
    scan_ctl(dirname, filename, argc, argv, "QNT_FORMAT", iq, "%g",
	     ctl->qnt_format[iq]);

    /* Try to identify quantities... */
    if (strcmp(ctl->qnt_name[iq], "mass") == 0)
      ctl->qnt_mass = iq;
    else if (strcmp(ctl->qnt_name[iq], "part_radius") == 0)
      ctl->qnt_r_p = iq;
    else if (strcmp(ctl->qnt_name[iq], "part_density") == 0)
      ctl->qnt_rho_p = iq;
    else if (strcmp(ctl->qnt_name[iq], "station") == 0)
      ctl->qnt_station = iq;
    else if (strcmp(ctl->qnt_name[iq], "temperature") == 0)
      ctl->qnt_temp = iq;
  }

  /* Time steps of simulation... */
  ctl->direction =
    (int) scan_ctl(dirname, filename, argc, argv, "DIRECTION", -1, "1", NULL);
  if (ctl->direction != -1 && ctl->direction != 1)
    ERRMSG("Set DIRECTION to -1 or 1!");
  ctl->t_start =
    scan_ctl(dirname, filename, argc, argv, "T_START", -1, "-1e100", NULL);
  ctl->t_stop =
    scan_ctl(dirname, filename, argc, argv, "T_STOP", -1, "-1e100", NULL);
  ctl->dt_mod =
    scan_ctl(dirname, filename, argc, argv, "DT_MOD", -1, "600", NULL);

  /* Meteorological data... */
  ctl->dt_met =
    scan_ctl(dirname, filename, argc, argv, "DT_MET", -1, "21600", NULL);
  ctl->red_met =
    (int) scan_ctl(dirname, filename, argc, argv, "RED_MET", -1, "1", NULL);

  /* Dispersion parameters... */
  ctl->turb_dx =
    scan_ctl(dirname, filename, argc, argv, "TURB_DX", -1, "50.0", NULL);
  ctl->turb_dz =
    scan_ctl(dirname, filename, argc, argv, "TURB_DZ", -1, "0.1", NULL);
  ctl->turb_meso =
    scan_ctl(dirname, filename, argc, argv, "TURB_MESO", -1, "0.16", NULL);

  /* Half life time of particles... */
  ctl->t12 = scan_ctl(dirname, filename, argc, argv, "T12", -1, "0", NULL);

  /* Output of atmospheric data... */
  scan_ctl(dirname, filename, argc, argv, "ATM_BASENAME", -1, "atm",
	   ctl->atm_basename);
  ctl->atm_dt_out =
    scan_ctl(dirname, filename, argc, argv, "ATM_DT_OUT", -1, "0", NULL);
  ctl->atm_iformat =
    (int) scan_ctl(dirname, filename, argc, argv, "ATM_IFORMAT", -1, "0",
		   NULL);
  ctl->atm_oformat =
    (int) scan_ctl(dirname, filename, argc, argv, "ATM_OFORMAT", -1, "0",
		   NULL);

  /* Output of CSI data... */
  scan_ctl(dirname, filename, argc, argv, "CSI_BASENAME", -1, "csi",
	   ctl->csi_basename);
  ctl->csi_dt_out =
    scan_ctl(dirname, filename, argc, argv, "CSI_DT_OUT", -1, "0", NULL);
  ctl->csi_dt_update =
    scan_ctl(dirname, filename, argc, argv, "CSI_DT_UPDATE", -1, "0", NULL);
  scan_ctl(dirname, filename, argc, argv, "CSI_OBSFILE", -1, "obs.tab",
	   ctl->csi_obsfile);
  ctl->csi_obsmin =
    scan_ctl(dirname, filename, argc, argv, "CSI_OBSMIN", -1, "0", NULL);
  ctl->csi_modmin =
    scan_ctl(dirname, filename, argc, argv, "CSI_MODMIN", -1, "0", NULL);
  ctl->csi_z0 =
    scan_ctl(dirname, filename, argc, argv, "CSI_Z0", -1, "0", NULL);
  ctl->csi_z1 =
    scan_ctl(dirname, filename, argc, argv, "CSI_Z1", -1, "100", NULL);
  ctl->csi_nz =
    (int) scan_ctl(dirname, filename, argc, argv, "CSI_NZ", -1, "1", NULL);
  ctl->csi_lon0 =
    scan_ctl(dirname, filename, argc, argv, "CSI_LON0", -1, "-180", NULL);
  ctl->csi_lon1 =
    scan_ctl(dirname, filename, argc, argv, "CSI_LON1", -1, "180", NULL);
  ctl->csi_nx =
    (int) scan_ctl(dirname, filename, argc, argv, "CSI_NX", -1, "360", NULL);
  ctl->csi_lat0 =
    scan_ctl(dirname, filename, argc, argv, "CSI_LAT0", -1, "-90", NULL);
  ctl->csi_lat1 =
    scan_ctl(dirname, filename, argc, argv, "CSI_LAT1", -1, "90", NULL);
  ctl->csi_ny =
    (int) scan_ctl(dirname, filename, argc, argv, "CSI_NY", -1, "180", NULL);

  /* Spearman */
  ctl->spearman_modmin =
    scan_ctl(dirname, filename, argc, argv, "SPEARMAN_MODMIN", -1, "-9999999",
	     NULL);
  ctl->spearman_obsmin =
    scan_ctl(dirname, filename, argc, argv, "SPEARMAN_OBSMIN", -1, "-9999999",
	     NULL);

  /* Output of gridded data... */
  scan_ctl(dirname, filename, argc, argv, "GRID_BASENAME", -1, "grid",
	   ctl->grid_basename);
  ctl->grid_dt_out =
    scan_ctl(dirname, filename, argc, argv, "GRID_DT_OUT", -1, "0", NULL);
  ctl->grid_z0 =
    scan_ctl(dirname, filename, argc, argv, "GRID_Z0", -1, "0", NULL);
  ctl->grid_z1 =
    scan_ctl(dirname, filename, argc, argv, "GRID_Z1", -1, "100", NULL);
  ctl->grid_nz =
    (int) scan_ctl(dirname, filename, argc, argv, "GRID_NZ", -1, "1", NULL);
  ctl->grid_lon0 =
    scan_ctl(dirname, filename, argc, argv, "GRID_LON0", -1, "-180", NULL);
  ctl->grid_lon1 =
    scan_ctl(dirname, filename, argc, argv, "GRID_LON1", -1, "180", NULL);
  ctl->grid_nx =
    (int) scan_ctl(dirname, filename, argc, argv, "GRID_NX", -1, "360", NULL);
  ctl->grid_lat0 =
    scan_ctl(dirname, filename, argc, argv, "GRID_LAT0", -1, "-90", NULL);
  ctl->grid_lat1 =
    scan_ctl(dirname, filename, argc, argv, "GRID_LAT1", -1, "90", NULL);
  ctl->grid_ny =
    (int) scan_ctl(dirname, filename, argc, argv, "GRID_NY", -1, "180", NULL);

  /* Output of station data... */
  scan_ctl(dirname, filename, argc, argv, "STAT_BASENAME", -1, "stat",
	   ctl->stat_basename);
  ctl->stat_dt_out =
    scan_ctl(dirname, filename, argc, argv, "STAT_DT_OUT", -1, "0", NULL);
  ctl->stat_lon =
    scan_ctl(dirname, filename, argc, argv, "STAT_LON", -1, "0", NULL);
  ctl->stat_lat =
    scan_ctl(dirname, filename, argc, argv, "STAT_LAT", -1, "0", NULL);
  ctl->stat_r =
    scan_ctl(dirname, filename, argc, argv, "STAT_R", -1, "100", NULL);
}

/*****************************************************************************/

void read_met(
  char *filename,
  met_t * met) {

  int ip, dimid, ncid, varid, year, mon, day, hour;

  size_t np, np2 = 1, nx, ny;

  /* Write info... */
  printf("Read meteorological data: %s\n", filename);

  /* Open netCDF file... */
  NC(nc_open(filename, NC_NOWRITE, &ncid));

  /* Get dimensions... */
  NC(nc_inq_dimid(ncid, "lon", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &nx));
  if (nx > EX)
    ERRMSG("Too many longitudes!");

  NC(nc_inq_dimid(ncid, "lat", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &ny));
  if (ny > EY)
    ERRMSG("Too many latitudes!");

  NC(nc_inq_dimid(ncid, "lev", &dimid));
  NC(nc_inq_dimlen(ncid, dimid, &np));
  if (np > EP)
    ERRMSG("Too many pressure levels!");

  /* Check for different number of pressure levels... */
  if (nc_inq_dimid(ncid, "lev_2", &dimid) == NC_NOERR)
    NC(nc_inq_dimlen(ncid, dimid, &np2));
  if (np2 == 1)
    np2 = np;

  /* Store dimensions... */
  met->np = (int) np;
  met->nx = (int) nx;
  met->ny = (int) ny;

  /* Read geolocations... */
  NC(nc_inq_varid(ncid, "time", &varid));
  NC(nc_get_var_double(ncid, varid, &met->time));

  NC(nc_inq_varid(ncid, "lev", &varid));
  NC(nc_get_var_double(ncid, varid, met->p));

  NC(nc_inq_varid(ncid, "lon", &varid));
  NC(nc_get_var_double(ncid, varid, met->lon));

  NC(nc_inq_varid(ncid, "lat", &varid));
  NC(nc_get_var_double(ncid, varid, met->lat));

  /* Convert time... */
  year = (int) met->time / 10000;
  met->time -= year * 10000;
  mon = (int) met->time / 100;
  met->time -= mon * 100;
  day = (int) (met->time);
  met->time -= day;
  hour = (int) (met->time * 24.);
  time2jsec(year, mon, day, hour, 0, 0, 0, &met->time);

  /* Convert pressure to hPa... */
  for (ip = 0; ip < met->np; ip++)
    met->p[ip] /= 100.;

  /* Read wind... */
  read_met_help(ncid, "T", met, met->np, met->t, 1.0);
  read_met_help(ncid, "U", met, met->np, met->u, 1.0);
  read_met_help(ncid, "V", met, met->np, met->v, 1.0);
  read_met_help(ncid, "W", met, (int) np2, met->w, 0.01f);

  /* Close file... */
  NC(nc_close(ncid));
}

/*****************************************************************************/

void read_met_help(
  int ncid,
  char *varname,
  met_t * met,
  int np,
  float dest[EX][EY][EP],
  float scl) {

  static float help[EX * EY * EP];

  char varname2[LEN];

  int i, ip, ix, iy, n = 0, varid;

  /* Check if variable exists... */
  if (nc_inq_varid(ncid, varname, &varid) != NC_NOERR) {
    for (i = 0; i < (int) strlen(varname); i++) {
      varname2[i] = (char) tolower(varname[i]);
      varname2[i + 1] = '\0';
    }
    NC(nc_inq_varid(ncid, varname2, &varid));
  }
  NC(nc_get_var_float(ncid, varid, help));

  /* Copy data... */
  for (ip = 0; ip < np; ip++)
    for (iy = 0; iy < met->ny; iy++)
      for (ix = 0; ix < met->nx; ix++)
	dest[ix][iy][ip] = scl * help[n++];

  /* Check data... */
  for (ip = 0; ip < met->np; ip++)
    for (iy = 0; iy < met->ny; iy++)
      for (ix = 0; ix < met->nx; ix++)
	if (dest[ix][iy][ip] < -1e10 || dest[ix][iy][ip] > 1e10)
	  dest[ix][iy][ip] = GSL_NAN;
}

/*****************************************************************************/

void reduce_met(
  met_t * met,
  int dx,
  int dy,
  int dp) {

  met_t *met2;

  int ix, iy, ip;

  /* Check grid distances... */
  if (dx < 1 || dy < 1 || dp < 1 || dx * dy * dp <= 1)
    return;

  /* Allocate... */
  ALLOC(met2, met_t, 1);

  /* Reduce resolution... */
  for (ix = 0; ix < met->nx; ix++) {
    met2->lon[ix / dx] += met->lon[ix] / (double) (dx);
    for (iy = 0; iy < met->ny; iy++) {
      if (ix == 0)
	met2->lat[iy / dy] += met->lat[iy] / (double) (dy);
      for (ip = 0; ip < met->np; ip++) {
	if (ix == 0 && iy == 0)
	  met2->p[ip / dp] += met->p[ip] / (double) (dp);
	met2->t[ix / dx][iy / dy][ip / dp] +=
	  met->t[ix][iy][ip] / (float) (dx * dy * dp);
	met2->u[ix / dx][iy / dy][ip / dp] +=
	  met->u[ix][iy][ip] / (float) (dx * dy * dp);
	met2->v[ix / dx][iy / dy][ip / dp] +=
	  met->v[ix][iy][ip] / (float) (dx * dy * dp);
	met2->w[ix / dx][iy / dy][ip / dp] +=
	  met->w[ix][iy][ip] / (float) (dx * dy * dp);
      }
    }
  }
  met2->nx = met->nx / dx;
  met2->ny = met->ny / dy;
  met2->np = met->np / dp;

  /* Copy... */
  memcpy(met, met2, sizeof(met_t));

  /* Free... */
  free(met2);
}

/*****************************************************************************/

double scan_ctl(
  const char *dirname,
  const char *filename,
  int argc,
  char *argv[],
  const char *varname,
  int arridx,
  const char *defvalue,
  char *value) {

  FILE *in = NULL;

  char dummy[LEN], file[LEN], fullname1[LEN], fullname2[LEN], line[LEN],
    msg[LEN], rvarname[LEN], rval[LEN];

  int contain = 0, i;

  /* Set filename... */
  if (dirname != NULL && filename != NULL)
    sprintf(file, "%s/%s", dirname, filename);
  else if (filename != NULL)
    sprintf(file, "%s", filename);
  else
    sprintf(file, "%s", argv[1]);

  /* Open file... */
  if (file[strlen(file) - 1] != '-')
    if (!(in = fopen(file, "r")))
      ERRMSG("Cannot open file!");

  /* Set full variable name... */
  if (arridx >= 0) {
    sprintf(fullname1, "%s[%d]", varname, arridx);
    sprintf(fullname2, "%s[*]", varname);
  } else {
    sprintf(fullname1, "%s", varname);
    sprintf(fullname2, "%s", varname);
  }

  /* Read data... */
  if (in != NULL)
    while (fgets(line, LEN, in))
      if (sscanf(line, "%s %s %s", rvarname, dummy, rval) == 3)
	if (strcasecmp(rvarname, fullname1) == 0 ||
	    strcasecmp(rvarname, fullname2) == 0) {
	  contain = 1;
	  break;
	}
  for (i = 1; i < argc - 1; i++)
    if (strcasecmp(argv[i], fullname1) == 0 ||
	strcasecmp(argv[i], fullname2) == 0) {
      sprintf(rval, "%s", argv[i + 1]);
      contain = 1;
      break;
    }

  /* Close file... */
  if (in != NULL)
    fclose(in);

  /* Check for missing variables... */
  if (!contain) {
    if (strlen(defvalue) > 0)
      sprintf(rval, "%s", defvalue);
    else {
      sprintf(msg, "Missing variable %s!\n", fullname1);
      ERRMSG(msg);
    }
  }

  /* Write info... */
  printf("%s = %s\n", fullname1, rval);

  /* Return values... */
  if (value != NULL)
    sprintf(value, "%s", rval);
  return atof(rval);
}

/*****************************************************************************/

void time2jsec(
  int year,
  int mon,
  int day,
  int hour,
  int min,
  int sec,
  double remain,
  double *jsec) {

  struct tm t0, t1;

  t0.tm_year = 100;
  t0.tm_mon = 0;
  t0.tm_mday = 1;
  t0.tm_hour = 0;
  t0.tm_min = 0;
  t0.tm_sec = 0;

  t1.tm_year = year - 1900;
  t1.tm_mon = mon - 1;
  t1.tm_mday = day;
  t1.tm_hour = hour;
  t1.tm_min = min;
  t1.tm_sec = sec;

  *jsec = (double) timegm(&t1) - (double) timegm(&t0) + remain;
}

/*****************************************************************************/

int timer(
  const char *name,
  int id,
  int mode) {

  static struct timeval tim;

  static char names[NTIMER][LEN];

  static double starttime[NTIMER], runtime[NTIMER];

  static int n;

  /* Create new timer... */
  if (name != NULL) {
    if (n >= NTIMER)
      ERRMSG("Too many timers!");
    strcpy(names[n], name);
    id = n;
    n++;
  } else {

    /* Start timer... */
    if (mode == 1) {
      if (starttime[id] <= 0) {
	gettimeofday(&tim, NULL);
	starttime[id] = (double) tim.tv_sec + (double) tim.tv_usec / 1e6;
      } else
	ERRMSG("Timer already started!");
    }

    /* Stop timer... */
    else if (mode == 2) {
      if (starttime[id] > 0) {
	gettimeofday(&tim, NULL);
	runtime[id] =
	  runtime[id] + (double) tim.tv_sec + (double) tim.tv_usec / 1e6 -
	  starttime[id];
	starttime[id] = -1;
      } else
	ERRMSG("Timer not started!");
    }

    /* Print timer... */
    else if (mode == 3) {
      if (starttime[id] > 0) {
	gettimeofday(&tim, NULL);
	printf("TIMER: %s %.4f s\n", names[id],
	       runtime[id] + (double) tim.tv_sec +
	       (double) tim.tv_usec / 1e6 - starttime[id]);
      } else
	printf("TIMER: %s %.4f s\n", names[id], runtime[id]);
    }
  }

  return id;
}

/*****************************************************************************/

void write_atm(
  const char *dirname,
  const char *filename,
  atm_t * atm,
  ctl_t * ctl) {

  char file[LEN];

  /* Set filename... */
  if (dirname != NULL)
    sprintf(file, "%s/%s", dirname, filename);
  else
    sprintf(file, "%s", filename);

  /* Write info... */
  printf("Write atmospheric data: %s\n", file);

  /* Select file format... */
  if (ctl->atm_oformat == 0)
    write_atm_to_ascii(file, atm, ctl);
  else if (ctl->atm_oformat == 1)
    write_atm_to_netcdf(file, atm, ctl);
  else
    ERRMSG("Unknown file format!");
}

/*****************************************************************************/

void write_atm_to_ascii(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl) {

  static FILE *out;

  int ip, iq;

  /* Create file... */
  if (!(out = fopen(filename, "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = altitude [km]\n"
	  "# $3 = longitude [deg]\n" "# $4 = latitude [deg]\n");
  for (iq = 0; iq < ctl->nq; iq++)
    fprintf(out, "# $%i = %s [%s]\n", iq + 5, ctl->qnt_name[iq],
	    ctl->qnt_unit[iq]);
  fprintf(out, "\n");

  /* Write data... */
  for (ip = 0; ip < atm->np; ip++) {
    fprintf(out, "%.2f %g %g %g", atm->time[ip], Z(atm->p[ip]),
	    atm->lon[ip], atm->lat[ip]);
    for (iq = 0; iq < ctl->nq; iq++) {
      fprintf(out, " ");
      fprintf(out, ctl->qnt_format[iq], atm->q[iq][ip]);
    }
    fprintf(out, "\n");
  }

  /* Close file... */
  fclose(out);
}

/*****************************************************************************/

void write_atm_to_netcdf(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl) {

  static int n_dimid, dimids[1], iq, varid[4 + NQ], ncid = -1;

  size_t start[1], count[1];

  char unit_str[LEN];

  /* Create file... */
  NC(nc_create(filename, NC_CLOBBER, &ncid));

  /* Create dimensions... */
  NC(nc_def_dim(ncid, "np", NC_UNLIMITED, &n_dimid));
  dimids[0] = n_dimid;

  /* Create variables... */
  NC(nc_def_var(ncid, "time", NC_DOUBLE, 1, dimids, &varid[0]));
  NC(nc_put_att_text(ncid, varid[0], "units", 1, "s"));

  NC(nc_def_var(ncid, "pressure", NC_DOUBLE, 1, dimids, &varid[1]));
  NC(nc_put_att_text(ncid, varid[1], "units", 3, "hPa"));

  NC(nc_def_var(ncid, "longitude", NC_DOUBLE, 1, dimids, &varid[2]));
  NC(nc_put_att_text(ncid, varid[2], "units", 3, "deg"));

  NC(nc_def_var(ncid, "latitude", NC_DOUBLE, 1, dimids, &varid[3]));
  NC(nc_put_att_text(ncid, varid[3], "units", 3, "deg"));

  /* Loop over quantities... */
  for (iq = 0; iq < ctl->nq; iq++) {
    NC(nc_def_var
       (ncid, ctl->qnt_name[iq], NC_DOUBLE, 1, dimids, &varid[4 + iq]));
    strcpy(unit_str, ctl->qnt_unit[iq]);
    NC(nc_put_att_text
       (ncid, varid[4 + iq], "units", strlen(unit_str), unit_str));
  }

  /* Finish definition... */
  nc_enddef(ncid);

  /* Write data... */
  start[0] = (size_t) 0;
  count[0] = (size_t) atm->np;
  NC(nc_put_vara_double(ncid, varid[0], start, count, atm->time));
  NC(nc_put_vara_double(ncid, varid[1], start, count, atm->p));
  NC(nc_put_vara_double(ncid, varid[2], start, count, atm->lon));
  NC(nc_put_vara_double(ncid, varid[3], start, count, atm->lat));
  for (iq = 0; iq < ctl->nq; iq++)
    NC(nc_put_vara_double(ncid, varid[4 + iq], start, count, atm->q[iq]));

  /* Close file... */
  NC(nc_close(ncid));
}

/*****************************************************************************/

void write_csi(
  const char *dirname,
  const char *filename,
  atm_t * atm,
  ctl_t * ctl,
  double t,
  double dt,
  int write) {

  static FILE *in, *out;

  static char line[LEN];

  static double modmean[GX][GY][GZ], obsmean[GX][GY][GZ],
    rt, rz, rlon, rlat, t0, t1, area, dlon, dlat, lat;

  static double mod_val[SPEARMAN_NP], obs_val[SPEARMAN_NP],
    work[2 * SPEARMAN_NP];

  static int init = 0, obscounter[GX][GY][GZ],
    gw[GX][GY][GZ], gx[GX][GY][GZ], gy[GX][GY][GZ], gz[GX][GY][GZ], cw, cx,
    cy, cz, ip, ix, iy, iz, robs;

  static size_t mod_n;

  char file[LEN];

  /* Set filename... */
  if (dirname != NULL)
    sprintf(file, "%s/%s", dirname, filename);
  else
    sprintf(file, "%s", filename);

  /* Write info... */
  if (write)
    printf("Write CSI data: %s\n", file);
  else
    printf("Update CSI data: %s\n", file);

  /* Check quantity index for mass... */
  if (ctl->qnt_mass < 0)
    ERRMSG("Need quantity mass to analyze CSI!");

  /* Set time interval... */
  t0 = t - 0.5 * dt;
  t1 = t + 0.5 * dt;

  /* Initialize grid cells... */
  for (ix = 0; ix < ctl->csi_nx; ix++)
    for (iy = 0; iy < ctl->csi_ny; iy++)
      for (iz = 0; iz < ctl->csi_nz; iz++) {
	modmean[ix][iy][iz] = 0;
	obsmean[ix][iy][iz] = 0;
	obscounter[ix][iy][iz] = 0;
      }

  /* Open observation data file... */
  if (!init) {
    init = 1;
    if (!(in = fopen(ctl->csi_obsfile, "r")))
      ERRMSG("Cannot open observation data file!");
  }

  /* Read data... */
  while (fgets(line, LEN, in)) {

    /* Read data... */
    if (sscanf(line, "%lg %lg %lg %lg %d", &rt, &rz, &rlon, &rlat, &robs) !=
	5)
      continue;

    /* Check time... */
    if (rt < t0)
      continue;
    if (rt > t1)
      break;

    /* Calculate indices */
    ix = (int) ((rlon - ctl->csi_lon0)
		/ (ctl->csi_lon1 - ctl->csi_lon0) * ctl->csi_nx);
    iy = (int) ((rlat - ctl->csi_lat0)
		/ (ctl->csi_lat1 - ctl->csi_lat0) * ctl->csi_ny);
    iz = (int) ((rz - ctl->csi_z0)
		/ (ctl->csi_z1 - ctl->csi_z0) * ctl->csi_nz);

    /* Check indices... */
    if (ix < 0 || ix >= ctl->csi_nx ||
	iy < 0 || iy >= ctl->csi_ny || iz < 0 || iz >= ctl->csi_nz)
      continue;

    /* Get mean observation index... */
    obsmean[ix][iy][iz] += robs;
    obscounter[ix][iy][iz]++;
  }

  /* Analyze model data... */
  for (ip = 0; ip < atm->np; ip++) {

    /* Check time... */
    if (atm->time[ip] < t0 || atm->time[ip] > t1)
      continue;

    /* Get indices... */
    ix = (int) ((atm->lon[ip] - ctl->csi_lon0)
		/ (ctl->csi_lon1 - ctl->csi_lon0) * ctl->csi_nx);
    iy = (int) ((atm->lat[ip] - ctl->csi_lat0)
		/ (ctl->csi_lat1 - ctl->csi_lat0) * ctl->csi_ny);
    iz = (int) ((Z(atm->p[ip]) - ctl->csi_z0)
		/ (ctl->csi_z1 - ctl->csi_z0) * ctl->csi_nz);

    /* Check indices... */
    if (ix < 0 || ix >= ctl->csi_nx ||
	iy < 0 || iy >= ctl->csi_ny || iz < 0 || iz >= ctl->csi_nz)
      continue;

    /* Get total mass in grid cell... */
    modmean[ix][iy][iz] += atm->q[ctl->qnt_mass][ip];
  }

  /* Analyze all grid cells... */
  for (ix = 0; ix < ctl->csi_nx; ix++)
    for (iy = 0; iy < ctl->csi_ny; iy++)
      for (iz = 0; iz < ctl->csi_nz; iz++) {

	/* Calculate mean observation index... */
	if (obscounter[ix][iy][iz] > 0)
	  obsmean[ix][iy][iz] /= obscounter[ix][iy][iz];

	/* Calculate column density... */
	if (modmean[ix][iy][iz] > 0) {
	  dlon = (ctl->csi_lon1 - ctl->csi_lon0) / ctl->csi_nx;
	  dlat = (ctl->csi_lat1 - ctl->csi_lat0) / ctl->csi_ny;
	  lat = ctl->csi_lat0 + dlat * (iy + 0.5);
	  area = dlat * M_PI * RE / 180. * dlon * M_PI * RE / 180.
	    * cos(lat * M_PI / 180.);
	  modmean[ix][iy][iz] /= (1e6 * area);
	}

	/* Calculate CSI... */
	if (obscounter[ix][iy][iz] > 0) {
	  if (obsmean[ix][iy][iz] >= ctl->csi_obsmin &&
	      modmean[ix][iy][iz] >= ctl->csi_modmin) {
	    cx++;
	    gx[ix][iy][iz]++;
	  } else if (obsmean[ix][iy][iz] >= ctl->csi_obsmin &&
		     modmean[ix][iy][iz] < ctl->csi_modmin) {
	    cy++;
	    gy[ix][iy][iz]++;
	  } else if (obsmean[ix][iy][iz] < ctl->csi_obsmin &&
		     modmean[ix][iy][iz] >= ctl->csi_modmin) {
	    cz++;
	    gz[ix][iy][iz]++;
	  } else {
	    cw++;
	    gw[ix][iy][iz]++;
	  }
	}

	/* Fill Array for Spearman */
	if (mod_n < SPEARMAN_NP && obscounter[ix][iy][iz] > 0
	    && modmean[ix][iy][iz] >= ctl->spearman_modmin
	    && obsmean[ix][iy][iz] >= ctl->spearman_obsmin) {
	  mod_val[mod_n] = modmean[ix][iy][iz];
	  obs_val[mod_n] = obsmean[ix][iy][iz];
	  mod_n = mod_n + 1;
	}
      }

  /* Check if output is requested... */
  if (write) {

    /* Create CSI data file... */
    if (!(out = fopen(file, "w")))
      ERRMSG("Cannot create CSI data file!");

    /* Write header... */
    fprintf(out,
	    "# $1 = time [s]\n"
	    "# $2 = altitude [km]\n"
	    "# $3 = longitude [deg]\n"
	    "# $4 = latitude [deg]\n"
	    "# $5 = W\n" "# $6 = X\n" "# $7 = Y\n" "# $8 = Z\n");

    /* Write data... */
    for (iz = 0; iz < ctl->csi_nz; iz++)
      for (ix = 0; ix < ctl->csi_nx; ix++) {
	fprintf(out, "\n");
	for (iy = 0; iy < ctl->csi_ny; iy++)
	  fprintf(out, "%.2f %g %g %g %d %d %d %d\n",
		  t,
		  ctl->csi_z0 + (ctl->csi_z1 -
				 ctl->csi_z0) / ctl->csi_nz * (iz + 0.5),
		  ctl->csi_lon0 + (ctl->csi_lon1 -
				   ctl->csi_lon0) / ctl->csi_nx * (ix + 0.5),
		  ctl->csi_lat0 + (ctl->csi_lat1 -
				   ctl->csi_lat0) / ctl->csi_ny * (iy + 0.5),
		  gw[ix][iy][iz], gx[ix][iy][iz], gy[ix][iy][iz],
		  gz[ix][iy][iz]);
      }

    /* Print CSI data... */
    fprintf(out, "\n");
    fprintf(out, "# t= %.2f / O= %d\n", t, cx + cy);
    fprintf(out, "# t= %.2f / F= %d\n", t, cx + cz);
    fprintf(out, "# t= %.2f / X= %d\n", t, cx);
    fprintf(out, "# t= %.2f / Y= %d\n", t, cy);
    fprintf(out, "# t= %.2f / Z= %d\n", t, cz);
    fprintf(out, "# t= %.2f / B= %g\n", t, 100. * (cx + cz) / (cx + cy));
    fprintf(out, "# t= %.2f / POD= %g\n", t, (100. * cx) / (cx + cy));
    fprintf(out, "# t= %.2f / FAR= %g\n", t, (100. * cz) / (cx + cz));
    fprintf(out, "# t= %.2f / CSI= %g\n", t, (100. * cx) / (cx + cy + cz));
    fprintf(out, "# t= %.2f / GS= %g\n", t,
	    (100. * (cx * cw - cy * cz)) / ((cy + cz) * (cx + cy + cw + cz) *
					    (cx * cw - cy * cz)));
    fprintf(out, "# t= %.2f / SPEARMAN= %.8f\n", t,
	    mod_n > 0 ? gsl_stats_spearman(mod_val, 1, obs_val, 1, mod_n,
					   work) : 0);
    fprintf(out, "# t= %.2f / SPEARMAN_N= %i\n", t, (int) mod_n);

    /* Close data file... */
    fclose(out);

    /* Reset counters... */
    mod_n = 0;
    cw = cx = cy = cz = 0;
    for (ix = 0; ix < ctl->csi_nx; ix++)
      for (iy = 0; iy < ctl->csi_ny; iy++)
	for (iz = 0; iz < ctl->csi_nz; iz++)
	  gw[ix][iy][iz] = gx[ix][iy][iz] = gy[ix][iy][iz] = gz[ix][iy][iz] =
	    0;
  }
}

/*****************************************************************************/

void write_grid(
  const char *dirname,
  const char *filename,
  atm_t * atm,
  ctl_t * ctl,
  double t,
  double dt) {

  static double grid_q[NQ][GX][GY][GZ], z, dz, lon, dlon, lat, dlat,
    area, rho_air, cd, mmr, t0, t1;

  static int grid_n[GX][GY][GZ], ip, ix, iy, iz, iq;

  FILE *out;

  char file[LEN];

  /* Set filename... */
  if (dirname != NULL)
    sprintf(file, "%s/%s", dirname, filename);
  else
    sprintf(file, "%s", filename);

  /* Check dimensions... */
  if (ctl->grid_nx > GX || ctl->grid_ny > GY || ctl->grid_nz > GZ)
    ERRMSG("Grid dimensions too large!");

  /* Set time interval for output... */
  t0 = t - 0.5 * dt;
  t1 = t + 0.5 * dt;

  /* Set grid box sizes... */
  dz = (ctl->grid_z1 - ctl->grid_z0) / ctl->grid_nz;
  dlon = (ctl->grid_lon1 - ctl->grid_lon0) / ctl->grid_nx;
  dlat = (ctl->grid_lat1 - ctl->grid_lat0) / ctl->grid_ny;

  /* Initialize grid... */
  for (ix = 0; ix < ctl->grid_nx; ix++)
    for (iy = 0; iy < ctl->grid_ny; iy++)
      for (iz = 0; iz < ctl->grid_nz; iz++) {
	grid_n[ix][iy][iz] = 0;
	for (iq = 0; iq < ctl->nq; iq++)
	  grid_q[iq][ix][iy][iz] = 0;
      }

  /* Average data... */
  for (ip = 0; ip < atm->np; ip++)
    if (atm->time[ip] >= t0 && atm->time[ip] <= t1) {
      ix = (int) ((atm->lon[ip] - ctl->grid_lon0) / dlon);
      iy = (int) ((atm->lat[ip] - ctl->grid_lat0) / dlat);
      iz = (int) ((Z(atm->p[ip]) - ctl->grid_z0) / dz);

      /* Check indices... */
      if (ix < 0 || ix >= ctl->grid_nx ||
	  iy < 0 || iy >= ctl->grid_ny || iz < 0 || iz >= ctl->grid_nz)
	continue;

      /* Add data... */
      grid_n[ix][iy][iz]++;
      for (iq = 0; iq < ctl->nq; iq++)
	grid_q[iq][ix][iy][iz] += atm->q[iq][ip];
    }

  /* Calculate mean values... */
  for (iq = 0; iq < ctl->nq; iq++)
    for (ix = 0; ix < ctl->grid_nx; ix++)
      for (iy = 0; iy < ctl->grid_ny; iy++)
	for (iz = 0; iz < ctl->grid_nz; iz++)
	  if (grid_n[ix][iy][iz] > 0)
	    grid_q[iq][ix][iy][iz] /= grid_n[ix][iy][iz];

  /* Write info... */
  printf("Write gridded data: %s\n", file);

  /* Create new file... */
  if (!(out = fopen(file, "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = altitude [km]\n"
	  "# $3 = longitude [deg]\n"
	  "# $4 = latitude [deg]\n"
	  "# $5 = surface area [km^2]\n"
	  "# $6 = layer width [km]\n"
	  "# $7 = number of particles\n"
	  "# $8 = column density [kg/m^2]\n"
	  "# $9 = mass mixing ratio [kg/kg]\n");
  for (iq = 0; iq < ctl->nq; iq++)
    fprintf(out, "# $%i = %s (mean) [%s]\n", (10 + iq),
	    ctl->qnt_name[iq], ctl->qnt_unit[iq]);

  /* Write data... */
  for (iz = 0; iz < ctl->grid_nz; iz++)
    for (ix = 0; ix < ctl->grid_nx; ix++) {
      fprintf(out, "\n");
      for (iy = 0; iy < ctl->grid_ny; iy++) {

	/* Set coordinates... */
	z = ctl->grid_z0 + dz * (iz + 0.5);
	lon = ctl->grid_lon0 + dlon * (ix + 0.5);
	lat = ctl->grid_lat0 + dlat * (iy + 0.5);

	/* Calculate surface area... */
	area = dlat * M_PI * RE / 180. * dlon * M_PI * RE / 180.
	  * cos(lat * M_PI / 180.);

	/* Calculate column density... */
	if (ctl->qnt_mass >= 0)
	  cd =
	    grid_q[ctl->qnt_mass][ix][iy][iz] * grid_n[ix][iy][iz] / (1e6 *
								      area);

	/* Calculate mass mixing ratio... */
	if (ctl->qnt_mass >= 0 && ctl->qnt_temp >= 0) {
	  rho_air =
	    100. * P(z) / (287.058 * grid_q[ctl->qnt_temp][ix][iy][iz]);
	  mmr = grid_q[ctl->qnt_mass][ix][iy][iz] * grid_n[ix][iy][iz]
	    / (rho_air * 1e6 * area * 1e3 * dz);
	}

	/* Write to file... */
	fprintf(out, "%.2f %g %g %g %g %g %d %g %g",
		t, z, lon, lat, area, dz, grid_n[ix][iy][iz], cd, mmr);
	for (iq = 0; iq < ctl->nq; iq++)
	  fprintf(out, " %g", grid_q[iq][ix][iy][iz]);
	fprintf(out, "\n");
      }
    }

  /* Close file... */
  fclose(out);
}

/*****************************************************************************/

void write_station(
  const char *dirname,
  const char *filename,
  atm_t * atm,
  ctl_t * ctl,
  double t,
  double dt) {

  static FILE *out;

  double rmax2, x0[3], x1[3], t0, t1;

  int ip, iq;

  char file[LEN];

  /* Set filename... */
  if (dirname != NULL)
    sprintf(file, "%s/%s", dirname, filename);
  else
    sprintf(file, "%s", filename);

  /* Init... */
  rmax2 = gsl_pow_2(ctl->stat_r);
  geo2cart(0, ctl->stat_lon, ctl->stat_lat, x0);

  /* Set time interval for output... */
  t0 = t - 0.5 * dt;
  t1 = t + 0.5 * dt;

  /* Write info... */
  printf("Write station data: %s\n", file);

  /* Create new file... */
  if (!(out = fopen(file, "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = altitude [km]\n"
	  "# $3 = longitude [deg]\n" "# $4 = latitude [deg]\n");
  for (iq = 0; iq < ctl->nq; iq++)
    fprintf(out, "# $%i = %s [%s]\n", (iq + 5),
	    ctl->qnt_name[iq], ctl->qnt_unit[iq]);
  fprintf(out, "\n");

  /* Loop over air parcels... */
  for (ip = 0; ip < atm->np; ip++)
    if (atm->time[ip] >= t0 && atm->time[ip] <= t1) {

      /* Check station flag... */
      if (ctl->qnt_station >= 0)
	if (atm->q[ctl->qnt_station][ip])
	  continue;

      /* Get Cartesian coordinates... */
      geo2cart(0, atm->lon[ip], atm->lat[ip], x1);

      /* Check horizontal distance... */
      if (DIST2(x0, x1) <= rmax2) {

	/* Set station flag... */
	if (ctl->qnt_station >= 0)
	  atm->q[ctl->qnt_station][ip] = 1;

	/* Write data... */
	fprintf(out, "%.2f %g %g %g",
		atm->time[ip], Z(atm->p[ip]), atm->lon[ip], atm->lat[ip]);
	for (iq = 0; iq < ctl->nq; iq++) {
	  fprintf(out, " ");
	  fprintf(out, ctl->qnt_format[iq], atm->q[iq][ip]);
	}
	fprintf(out, "\n");
      }
    }

  /* Close file... */
  fclose(out);
}
