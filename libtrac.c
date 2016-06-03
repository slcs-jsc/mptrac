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
  ctl_t * ctl,
  char *metbase,
  double t,
  met_t * met0,
  met_t * met1) {

  char filename[LEN];

  static int init;

  /* Initialize at first call... */
  if (!init) {
    init = 1;

    get_met_help(t, -1, metbase, ctl->dt_met, filename);
    read_met(filename, met0);

    get_met_help(t + 1.0 * ctl->direction, 1, metbase, ctl->dt_met, filename);
    read_met(filename, met1);
  }

  /* Read new data for forward trajectories... */
  if (t > met1->time && ctl->direction == 1) {
    memcpy(met0, met1, sizeof(met_t));
    get_met_help(t, 1, metbase, ctl->dt_met, filename);
    read_met(filename, met1);
  }

  /* Read new data for backward trajectories... */
  if (t < met0->time && ctl->direction == -1) {
    memcpy(met1, met0, sizeof(met_t));
    get_met_help(t, -1, metbase, ctl->dt_met, filename);
    read_met(filename, met0);
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

void intpol_met_2d(
  double array[EX][EY],
  int ix,
  int iy,
  double wx,
  double wy,
  double *var) {

  double aux00, aux01, aux10, aux11;

  /* Set variables... */
  aux00 = array[ix][iy];
  aux01 = array[ix][iy + 1];
  aux10 = array[ix + 1][iy];
  aux11 = array[ix + 1][iy + 1];

  /* Interpolate horizontally... */
  aux00 = wy * (aux00 - aux01) + aux01;
  aux11 = wy * (aux10 - aux11) + aux11;
  *var = wx * (aux00 - aux11) + aux11;
}

/*****************************************************************************/

void intpol_met_3d(
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
  double *ps,
  double *t,
  double *u,
  double *v,
  double *w,
  double *h2o,
  double *o3) {

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
  if (ps != NULL)
    intpol_met_2d(met->ps, ix, iy, wx, wy, ps);
  if (t != NULL)
    intpol_met_3d(met->t, ip, ix, iy, wp, wx, wy, t);
  if (u != NULL)
    intpol_met_3d(met->u, ip, ix, iy, wp, wx, wy, u);
  if (v != NULL)
    intpol_met_3d(met->v, ip, ix, iy, wp, wx, wy, v);
  if (w != NULL)
    intpol_met_3d(met->w, ip, ix, iy, wp, wx, wy, w);
  if (h2o != NULL)
    intpol_met_3d(met->h2o, ip, ix, iy, wp, wx, wy, h2o);
  if (o3 != NULL)
    intpol_met_3d(met->o3, ip, ix, iy, wp, wx, wy, o3);
}

/*****************************************************************************/

void intpol_met_time(
  met_t * met0,
  met_t * met1,
  double ts,
  double p,
  double lon,
  double lat,
  double *ps,
  double *t,
  double *u,
  double *v,
  double *w,
  double *h2o,
  double *o3) {

  double h2o0, h2o1, o30, o31, ps0, ps1, t0, t1, u0, u1, v0, v1, w0, w1, wt;

  /* Spatial interpolation... */
  intpol_met_space(met0, p, lon, lat,
		   ps == NULL ? NULL : &ps0,
		   t == NULL ? NULL : &t0,
		   u == NULL ? NULL : &u0,
		   v == NULL ? NULL : &v0,
		   w == NULL ? NULL : &w0,
		   h2o == NULL ? NULL : &h2o0, o3 == NULL ? NULL : &o30);
  intpol_met_space(met1, p, lon, lat,
		   ps == NULL ? NULL : &ps1,
		   t == NULL ? NULL : &t1,
		   u == NULL ? NULL : &u1,
		   v == NULL ? NULL : &v1,
		   w == NULL ? NULL : &w1,
		   h2o == NULL ? NULL : &h2o1, o3 == NULL ? NULL : &o31);

  /* Get weighting factor... */
  wt = (met1->time - ts) / (met1->time - met0->time);

  /* Interpolate... */
  if (ps != NULL)
    *ps = wt * (ps0 - ps1) + ps1;
  if (t != NULL)
    *t = wt * (t0 - t1) + t1;
  if (u != NULL)
    *u = wt * (u0 - u1) + u1;
  if (v != NULL)
    *v = wt * (v0 - v1) + v1;
  if (w != NULL)
    *w = wt * (w0 - w1) + w1;
  if (h2o != NULL)
    *h2o = wt * (h2o0 - h2o1) + h2o1;
  if (o3 != NULL)
    *o3 = wt * (o30 - o31) + o31;
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
  const char *filename,
  atm_t * atm,
  ctl_t * ctl) {

  FILE *in;

  char line[LEN], *tok;

  int iq;

  /* Initialize... */
  atm->np = 0;

  /* Write info... */
  printf("Read atmospheric data: %s\n", filename);

  /* Open file... */
  if (!(in = fopen(filename, "r")))
    ERRMSG("Cannot open file!");

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

void read_ctl(
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
  ctl->qnt_m = -1;
  ctl->qnt_r = -1;
  ctl->qnt_rho = -1;
  ctl->qnt_ps = -1;
  ctl->qnt_t = -1;
  ctl->qnt_u = -1;
  ctl->qnt_v = -1;
  ctl->qnt_w = -1;
  ctl->qnt_h2o = -1;
  ctl->qnt_o3 = -1;
  ctl->qnt_theta = -1;
  ctl->qnt_stat = -1;

  /* Read quantities... */
  ctl->nq = (int) scan_ctl(filename, argc, argv, "NQ", -1, "0", NULL);
  for (iq = 0; iq < ctl->nq; iq++) {

    /* Read quantity names, units, and formats... */
    scan_ctl(filename, argc, argv, "QNT_NAME", iq, "", ctl->qnt_name[iq]);
    scan_ctl(filename, argc, argv, "QNT_UNIT", iq, "", ctl->qnt_unit[iq]);
    scan_ctl(filename, argc, argv, "QNT_FORMAT", iq, "%g",
	     ctl->qnt_format[iq]);

    /* Try to identify quantities... */
    if (strcmp(ctl->qnt_name[iq], "m") == 0)
      ctl->qnt_m = iq;
    else if (strcmp(ctl->qnt_name[iq], "r") == 0)
      ctl->qnt_r = iq;
    else if (strcmp(ctl->qnt_name[iq], "rho") == 0)
      ctl->qnt_rho = iq;
    else if (strcmp(ctl->qnt_name[iq], "ps") == 0)
      ctl->qnt_ps = iq;
    else if (strcmp(ctl->qnt_name[iq], "t") == 0)
      ctl->qnt_t = iq;
    else if (strcmp(ctl->qnt_name[iq], "u") == 0)
      ctl->qnt_u = iq;
    else if (strcmp(ctl->qnt_name[iq], "v") == 0)
      ctl->qnt_v = iq;
    else if (strcmp(ctl->qnt_name[iq], "w") == 0)
      ctl->qnt_w = iq;
    else if (strcmp(ctl->qnt_name[iq], "h2o") == 0)
      ctl->qnt_h2o = iq;
    else if (strcmp(ctl->qnt_name[iq], "o3") == 0)
      ctl->qnt_o3 = iq;
    else if (strcmp(ctl->qnt_name[iq], "theta") == 0)
      ctl->qnt_theta = iq;
    else if (strcmp(ctl->qnt_name[iq], "stat") == 0)
      ctl->qnt_stat = iq;
  }

  /* Time steps of simulation... */
  ctl->direction =
    (int) scan_ctl(filename, argc, argv, "DIRECTION", -1, "1", NULL);
  if (ctl->direction != -1 && ctl->direction != 1)
    ERRMSG("Set DIRECTION to -1 or 1!");
  ctl->t_start =
    scan_ctl(filename, argc, argv, "T_START", -1, "-1e100", NULL);
  ctl->t_stop = scan_ctl(filename, argc, argv, "T_STOP", -1, "-1e100", NULL);
  ctl->dt_mod = scan_ctl(filename, argc, argv, "DT_MOD", -1, "600", NULL);

  /* Meteorological data... */
  ctl->dt_met = scan_ctl(filename, argc, argv, "DT_MET", -1, "21600", NULL);

  /* Isosurface parameter... */
  ctl->isosurf
    = (int) scan_ctl(filename, argc, argv, "ISOSURF", -1, "0", NULL);

  /* Diffusion parameters... */
  ctl->turb_dx_trop
    = scan_ctl(filename, argc, argv, "TURB_DX_TROP", -1, "50.0", NULL);
  ctl->turb_dx_strat
    = scan_ctl(filename, argc, argv, "TURB_DX_STRAT", -1, "0.0", NULL);
  ctl->turb_dz_trop
    = scan_ctl(filename, argc, argv, "TURB_DZ_TROP", -1, "0.0", NULL);
  ctl->turb_dz_strat
    = scan_ctl(filename, argc, argv, "TURB_DZ_STRAT", -1, "0.1", NULL);
  ctl->turb_meso =
    scan_ctl(filename, argc, argv, "TURB_MESO", -1, "0.16", NULL);

  /* Life time of particles... */
  ctl->tdec_trop = scan_ctl(filename, argc, argv, "TDEC_TROP", -1, "0", NULL);
  ctl->tdec_strat =
    scan_ctl(filename, argc, argv, "TDEC_STRAT", -1, "0", NULL);

  /* Output of atmospheric data... */
  scan_ctl(filename, argc, argv, "ATM_BASENAME", -1, "-", ctl->atm_basename);
  scan_ctl(filename, argc, argv, "ATM_GPFILE", -1, "-", ctl->atm_gpfile);
  ctl->atm_dt_out =
    scan_ctl(filename, argc, argv, "ATM_DT_OUT", -1, "86400", NULL);

  /* Output of CSI data... */
  scan_ctl(filename, argc, argv, "CSI_BASENAME", -1, "-", ctl->csi_basename);
  ctl->csi_dt_out =
    scan_ctl(filename, argc, argv, "CSI_DT_OUT", -1, "86400", NULL);
  scan_ctl(filename, argc, argv, "CSI_OBSFILE", -1, "obs.tab",
	   ctl->csi_obsfile);
  ctl->csi_obsmin =
    scan_ctl(filename, argc, argv, "CSI_OBSMIN", -1, "0", NULL);
  ctl->csi_modmin =
    scan_ctl(filename, argc, argv, "CSI_MODMIN", -1, "0", NULL);
  ctl->csi_z0 = scan_ctl(filename, argc, argv, "CSI_Z0", -1, "0", NULL);
  ctl->csi_z1 = scan_ctl(filename, argc, argv, "CSI_Z1", -1, "100", NULL);
  ctl->csi_nz = (int) scan_ctl(filename, argc, argv, "CSI_NZ", -1, "1", NULL);
  ctl->csi_lon0 =
    scan_ctl(filename, argc, argv, "CSI_LON0", -1, "-180", NULL);
  ctl->csi_lon1 = scan_ctl(filename, argc, argv, "CSI_LON1", -1, "180", NULL);
  ctl->csi_nx =
    (int) scan_ctl(filename, argc, argv, "CSI_NX", -1, "360", NULL);
  ctl->csi_lat0 = scan_ctl(filename, argc, argv, "CSI_LAT0", -1, "-90", NULL);
  ctl->csi_lat1 = scan_ctl(filename, argc, argv, "CSI_LAT1", -1, "90", NULL);
  ctl->csi_ny =
    (int) scan_ctl(filename, argc, argv, "CSI_NY", -1, "180", NULL);

  /* Output of gridded data... */
  scan_ctl(filename, argc, argv, "GRID_BASENAME", -1, "-",
	   ctl->grid_basename);
  scan_ctl(filename, argc, argv, "GRID_GPFILE", -1, "-", ctl->grid_gpfile);
  ctl->grid_dt_out =
    scan_ctl(filename, argc, argv, "GRID_DT_OUT", -1, "86400", NULL);
  ctl->grid_z0 = scan_ctl(filename, argc, argv, "GRID_Z0", -1, "0", NULL);
  ctl->grid_z1 = scan_ctl(filename, argc, argv, "GRID_Z1", -1, "100", NULL);
  ctl->grid_nz =
    (int) scan_ctl(filename, argc, argv, "GRID_NZ", -1, "1", NULL);
  ctl->grid_lon0 =
    scan_ctl(filename, argc, argv, "GRID_LON0", -1, "-180", NULL);
  ctl->grid_lon1 =
    scan_ctl(filename, argc, argv, "GRID_LON1", -1, "180", NULL);
  ctl->grid_nx =
    (int) scan_ctl(filename, argc, argv, "GRID_NX", -1, "360", NULL);
  ctl->grid_lat0 =
    scan_ctl(filename, argc, argv, "GRID_LAT0", -1, "-90", NULL);
  ctl->grid_lat1 =
    scan_ctl(filename, argc, argv, "GRID_LAT1", -1, "90", NULL);
  ctl->grid_ny =
    (int) scan_ctl(filename, argc, argv, "GRID_NY", -1, "180", NULL);

  /* Output of sample data... */
  scan_ctl(filename, argc, argv, "SAMPLE_BASENAME", -1, "-",
	   ctl->sample_basename);
  scan_ctl(filename, argc, argv, "SAMPLE_OBSFILE", -1, "-",
	   ctl->sample_obsfile);

  /* Output of station data... */
  scan_ctl(filename, argc, argv, "STAT_BASENAME", -1, "-",
	   ctl->stat_basename);
  ctl->stat_lon = scan_ctl(filename, argc, argv, "STAT_LON", -1, "0", NULL);
  ctl->stat_lat = scan_ctl(filename, argc, argv, "STAT_LAT", -1, "0", NULL);
  ctl->stat_r = scan_ctl(filename, argc, argv, "STAT_R", -1, "50", NULL);
}

/*****************************************************************************/

void read_met(
  char *filename,
  met_t * met) {

  static float help[EX * EY];

  int ix, iy, ip, dimid, ncid, varid, year, mon, day, hour;

  size_t np, nx, ny;

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

  /* Check and convert pressure levels... */
  for (ip = 0; ip < met->np; ip++) {
    if (ip > 0 && met->p[ip - 1] > met->p[ip])
      ERRMSG("Pressure levels must be in descending order!");
    met->p[ip] /= 100.;
  }

  /* Read surface pressure... */
  if (nc_inq_varid(ncid, "LNSP", &varid) == NC_NOERR) {
    NC(nc_get_var_float(ncid, varid, help));
    for (iy = 0; iy < met->ny; iy++)
      for (ix = 0; ix < met->nx; ix++)
	met->ps[ix][iy] = exp(help[iy * met->nx + ix]) / 100.;
  } else {
    for (ix = 0; ix < met->nx; ix++)
      for (iy = 0; iy < met->ny; iy++)
	met->ps[ix][iy] = met->p[0];
  }

  /* Read meteorological data... */
  read_met_help(ncid, "T", met, met->np, met->t, 1.0);
  read_met_help(ncid, "U", met, met->np, met->u, 1.0);
  read_met_help(ncid, "V", met, met->np, met->v, 1.0);
  read_met_help(ncid, "W", met, met->np, met->w, 0.01f);
  read_met_help(ncid, "Q", met, met->np, met->h2o, 1.608f);
  read_met_help(ncid, "O3", met, met->np, met->o3, 0.602f);

  /* Extrapolate data for lower boundary... */
  read_met_extrapolate(met);

  /* Copy data to obtain periodic boundary conditions... */
  read_met_periodic(met);

  /* Close file... */
  NC(nc_close(ncid));
}

/*****************************************************************************/

void read_met_extrapolate(
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
	met->h2o[ix][iy][ip]
	  = (float) LIN(met->p[ip + 1], met->h2o[ix][iy][ip + 1],
			met->p[ip + 2], met->h2o[ix][iy][ip + 2], met->p[ip]);
	met->o3[ix][iy][ip]
	  = (float) LIN(met->p[ip + 1], met->o3[ix][iy][ip + 1],
			met->p[ip + 2], met->o3[ix][iy][ip + 2], met->p[ip]);
      }
    }
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

  int ip, ix, iy, n = 0, varid;

  /* Check if variable exists... */
  if (nc_inq_varid(ncid, varname, &varid) != NC_NOERR)
    return;

  /* Read data... */
  NC(nc_get_var_float(ncid, varid, help));

  /* Copy and check data... */
  for (ip = 0; ip < np; ip++)
    for (iy = 0; iy < met->ny; iy++)
      for (ix = 0; ix < met->nx; ix++) {
	dest[ix][iy][ip] = scl * help[n++];
	if (dest[ix][iy][ip] < -1e10 || dest[ix][iy][ip] > 1e10)
	  dest[ix][iy][ip] = GSL_NAN;
      }
}

/*****************************************************************************/

void read_met_periodic(
  met_t * met) {

  int ip, iy;

  /* Check longitudes... */
  if (fabs(met->lon[met->nx - 1] - met->lon[0] - 360) < 0.01)
    return;

  /* Increase longitude counter... */
  if ((++met->nx) > EX)
    ERRMSG("Cannot create periodic boundary conditions!");

  /* Set longitude... */
  met->lon[met->nx - 1] = met->lon[met->nx - 2] + met->lon[1] - met->lon[0];

  /* Loop over latitudes and pressure levels... */
  for (iy = 0; iy < met->ny; iy++)
    for (ip = 0; ip < met->np; ip++) {
      met->ps[met->nx - 1][iy] = met->ps[0][iy];
      met->u[met->nx - 1][iy][ip] = met->u[0][iy][ip];
      met->v[met->nx - 1][iy][ip] = met->v[0][iy][ip];
      met->w[met->nx - 1][iy][ip] = met->w[0][iy][ip];
      met->t[met->nx - 1][iy][ip] = met->t[0][iy][ip];
      met->h2o[met->nx - 1][iy][ip] = met->h2o[0][iy][ip];
      met->o3[met->nx - 1][iy][ip] = met->o3[0][iy][ip];
    }
}

/*****************************************************************************/

double scan_ctl(
  const char *filename,
  int argc,
  char *argv[],
  const char *varname,
  int arridx,
  const char *defvalue,
  char *value) {

  FILE *in = NULL;

  char dummy[LEN], fullname1[LEN], fullname2[LEN], line[LEN],
    msg[LEN], rvarname[LEN], rval[LEN];

  int contain = 0, i;

  /* Open file... */
  if (filename[strlen(filename) - 1] != '-')
    if (!(in = fopen(filename, "r")))
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

void timer(
  const char *name,
  int id,
  int mode) {

  static double starttime[NTIMER], runtime[NTIMER];

  /* Check id... */
  if (id < 0 || id >= NTIMER)
    ERRMSG("Too many timers!");

  /* Start timer... */
  if (mode == 1) {
    if (starttime[id] <= 0)
      starttime[id] = omp_get_wtime();
    else
      ERRMSG("Timer already started!");
  }

  /* Stop timer... */
  else if (mode == 2) {
    if (starttime[id] > 0) {
      runtime[id] = runtime[id] + omp_get_wtime() - starttime[id];
      starttime[id] = -1;
    } else
      ERRMSG("Timer not started!");
  }

  /* Print timer... */
  else if (mode == 3)
    printf("%s = %g s\n", name, runtime[id]);
}

/*****************************************************************************/

double tropopause(
  double t,
  double lat) {

  static double doys[12]
  = { 1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 };

  static double lats[73]
    = { -90, -87.5, -85, -82.5, -80, -77.5, -75, -72.5, -70, -67.5,
    -65, -62.5, -60, -57.5, -55, -52.5, -50, -47.5, -45, -42.5,
    -40, -37.5, -35, -32.5, -30, -27.5, -25, -22.5, -20, -17.5,
    -15, -12.5, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 12.5,
    15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5,
    45, 47.5, 50, 52.5, 55, 57.5, 60, 62.5, 65, 67.5, 70, 72.5,
    75, 77.5, 80, 82.5, 85, 87.5, 90
  };

  static double tps[12][73]
    = { {324.1, 325.6, 325, 324.3, 322.5, 319.7, 314, 307.2, 301.8, 299.6,
	 297.1, 292.2, 285.6, 276.1, 264, 248.9, 231.9, 213.5, 194.4,
	 175.3, 157, 140.4, 126.7, 116.3, 109.5, 105.4, 103, 101.4, 100.4,
	 99.69, 99.19, 98.84, 98.56, 98.39, 98.39, 98.42, 98.44, 98.54,
	 98.68, 98.81, 98.89, 98.96, 99.12, 99.65, 101.4, 105.4, 113.5, 128,
	 152.1, 184.7, 214, 234.1, 247.3, 255.8, 262.6, 267.7, 271.7, 275,
	 277.2, 279, 280.1, 280.4, 280.6, 280.1, 279.3, 278.3, 276.8, 275.8,
	 275.3, 275.6, 275.4, 274.1, 273.5},
  {337.3, 338.7, 337.8, 336.4, 333, 328.8, 321.1, 312.6, 306.6, 303.7,
   300.2, 293.8, 285.4, 273.8, 259.6, 242.7, 224.4, 205.2, 186, 167.5,
   150.3, 135, 122.8, 113.9, 108.2, 104.7, 102.5, 101.1, 100.2, 99.42,
   98.88, 98.52, 98.25, 98.09, 98.07, 98.1, 98.12, 98.2, 98.25, 98.27,
   98.26, 98.27, 98.36, 98.79, 100.2, 104.2, 113.7, 131.2, 159.5, 193,
   220.4, 238.1, 250.2, 258.1, 264.7, 269.7, 273.7, 277.3, 280.2, 282.8,
   284.9, 286.5, 288.1, 288.8, 289, 288.5, 287.2, 286.3, 286.1, 287.2,
   287.5, 286.2, 285.8},
  {335, 336, 335.7, 335.1, 332.3, 328.1, 320.6, 311.8, 305.1, 301.9,
   297.6, 290, 280.4, 268.3, 254.6, 239.6, 223.9, 207.9, 192.2, 176.9,
   161.7, 146.4, 132.2, 120.6, 112.3, 107.2, 104.3, 102.4, 101.3,
   100.4, 99.86, 99.47, 99.16, 98.97, 98.94, 98.97, 99, 99.09, 99.2,
   99.31, 99.35, 99.41, 99.51, 99.86, 101.1, 104.9, 114.3, 131, 156.8,
   186.3, 209.3, 224.6, 236.8, 246.3, 254.9, 262.3, 268.8, 274.8,
   279.9, 284.6, 288.6, 291.6, 294.9, 297.5, 299.8, 301.8, 303.1,
   304.3, 304.9, 306, 306.6, 306.2, 306},
  {306.2, 306.7, 305.7, 307.1, 307.3, 306.4, 301.8, 296.2, 292.4,
   290.3, 287.1, 280.9, 273.4, 264.3, 254.1, 242.8, 231, 219, 207.2,
   195.5, 183.3, 169.7, 154.7, 138.7, 124.1, 113.6, 107.8, 104.7,
   102.8, 101.7, 100.9, 100.4, 100, 99.79, 99.7, 99.66, 99.68, 99.79,
   99.94, 100.2, 100.5, 100.9, 101.4, 102.1, 103.4, 107, 115.2, 129.1,
   148.7, 171, 190.8, 205.6, 218.4, 229.4, 239.6, 248.6, 256.5,
   263.7, 270.3, 276.6, 282.6, 288.1, 294.5, 300.4, 306.3, 311.4,
   315.1, 318.3, 320.3, 322.2, 322.8, 321.5, 321.1},
  {266.5, 264.9, 260.8, 261, 262, 263, 261.3, 259.7, 259.2, 259.8,
   260.1, 258.6, 256.7, 253.6, 249.5, 243.9, 237.4, 230, 222.1, 213.9,
   205, 194.4, 180.4, 161.8, 140.7, 122.9, 112.1, 106.7, 104.1, 102.7,
   101.8, 101.4, 101.1, 101, 101, 101, 101.1, 101.2, 101.5, 101.9,
   102.4, 103, 103.8, 104.9, 106.8, 110.1, 115.6, 124, 135.2, 148.9,
   165.2, 181.3, 198, 211.8, 223.5, 233.8, 242.9, 251.5, 259, 266.2,
   273.1, 279.2, 286.2, 292.8, 299.6, 306, 311.1, 315.5, 318.8, 322.6,
   325.3, 325.8, 325.8},
  {220.1, 218.1, 210.8, 207.2, 207.6, 210.5, 211.4, 213.5, 217.3,
   222.4, 227.9, 232.8, 237.4, 240.8, 242.8, 243, 241.5, 238.6, 234.2,
   228.5, 221, 210.7, 195.1, 172.9, 147.8, 127.6, 115.6, 109.9, 107.1,
   105.7, 105, 104.8, 104.8, 104.9, 105, 105.1, 105.3, 105.5, 105.8,
   106.4, 107, 107.6, 108.1, 108.8, 110, 111.8, 114.2, 117.4, 121.6,
   127.9, 137.3, 151.2, 169.5, 189, 205.8, 218.9, 229.1, 237.8, 245,
   251.5, 257.1, 262.3, 268.2, 274, 280.4, 286.7, 292.4, 297.9, 302.9,
   308.5, 312.2, 313.1, 313.3},
  {187.4, 184.5, 173.3, 166.1, 165.4, 167.8, 169.6, 173.6, 179.6,
   187.9, 198.9, 210, 220.5, 229.2, 235.7, 239.9, 241.8, 241.6, 239.6,
   235.8, 229.4, 218.6, 200.9, 175.9, 149.4, 129.4, 118.3, 113.1,
   110.8, 109.7, 109.3, 109.4, 109.7, 110, 110.2, 110.4, 110.5, 110.7,
   111, 111.4, 111.8, 112.1, 112.3, 112.7, 113.2, 113.9, 115, 116.4,
   117.9, 120.4, 124.1, 130.9, 142.2, 159.6, 179.6, 198.5, 212.9,
   224.2, 232.7, 239.1, 243.8, 247.7, 252.4, 257.3, 263.2, 269.5,
   275.4, 281.1, 286.3, 292, 296.3, 298.2, 298.8},
  {166, 166.4, 155.7, 148.3, 147.1, 149, 152.1, 157, 163.6, 172.4,
   185.3, 199.2, 212.6, 224, 233.2, 239.6, 243.3, 244.6, 243.6, 240.3,
   233.9, 222.6, 203.7, 177, 149.5, 129.7, 119, 114, 111.7, 110.7,
   110.3, 110.3, 110.6, 110.9, 111.1, 111.3, 111.5, 111.6, 111.9,
   112.2, 112.5, 112.6, 112.8, 113, 113.4, 114, 115.1, 116.5, 118.3,
   120.9, 124.4, 130.2, 139.4, 154.6, 173.8, 193.1, 208.1, 220.4,
   230.1, 238.2, 244.7, 249.5, 254.5, 259.3, 264.5, 269.4, 273.7,
   278.2, 282.6, 287.4, 290.9, 292.5, 293},
  {171.9, 172.8, 166.2, 162.3, 161.4, 162.5, 165.2, 169.6, 175.3,
   183.1, 193.8, 205.9, 218.3, 229.6, 238.5, 244.3, 246.9, 246.7,
   243.8, 238.4, 230.2, 217.9, 199.6, 174.9, 148.9, 129.8, 119.5,
   114.8, 112.3, 110.9, 110.3, 110.1, 110.2, 110.3, 110.4, 110.5,
   110.6, 110.8, 111, 111.4, 111.8, 112, 112.2, 112.4, 112.9, 113.6,
   114.7, 116.3, 118.4, 121.9, 127.1, 136.1, 149.8, 168.4, 186.9,
   203.3, 217, 229.1, 238.7, 247, 254, 259.3, 264.3, 268.3, 272.5,
   276.6, 280.4, 284.4, 288.4, 293.3, 297.2, 298.7, 299.1},
  {191.6, 192.2, 189, 188.1, 190.2, 193.7, 197.8, 202.9, 208.5,
   215.6, 224.2, 233.1, 241.2, 247.3, 250.8, 251.3, 248.9, 244.2,
   237.3, 228.4, 217.2, 202.9, 184.5, 162.5, 140.7, 124.8, 116.2,
   111.8, 109.4, 107.9, 107, 106.7, 106.6, 106.6, 106.7, 106.7,
   106.8, 107, 107.4, 108, 108.7, 109.3, 109.8, 110.4, 111.2,
   112.4, 114.2, 116.9, 121.1, 127.9, 139.3, 155.2, 173.6, 190.7,
   206.1, 220.1, 232.3, 243, 251.8, 259.2, 265.7, 270.6, 275.3,
   279.3, 283.3, 286.9, 289.7, 292.8, 296.1, 300.5, 303.9, 304.8,
   305.1},
  {241.5, 239.6, 236.8, 237.4, 239.4, 242.3, 244.2, 246.4, 249.2,
   253.6, 258.6, 262.7, 264.8, 264.2, 260.6, 254.1, 245.5, 235.3,
   223.9, 211.7, 198.3, 183.1, 165.6, 147.1, 130.5, 118.7, 111.9,
   108.1, 105.8, 104.3, 103.4, 102.8, 102.5, 102.4, 102.5, 102.5,
   102.5, 102.7, 103.1, 103.8, 104.6, 105.4, 106.1, 107, 108.2,
   109.9, 112.8, 117.5, 126, 140.4, 161, 181.9, 201.2, 216.8, 230.4,
   241.8, 251.4, 259.9, 266.9, 272.8, 277.4, 280.4, 282.9, 284.6,
   286.1, 287.4, 288.3, 289.5, 290.9, 294.2, 296.9, 297.5, 297.6},
  {301.2, 300.3, 296.6, 295.4, 295, 294.3, 291.2, 287.4, 284.9, 284.7,
   284.1, 281.5, 277.1, 270.4, 261.7, 250.6, 237.6, 223.1, 207.9, 192,
   175.8, 158.8, 142.1, 127.6, 116.8, 109.9, 106, 103.6, 102.1, 101.1,
   100.4, 99.96, 99.6, 99.37, 99.32, 99.32, 99.31, 99.46, 99.77, 100.2,
   100.7, 101.3, 101.8, 102.7, 104.1, 106.8, 111.9, 121, 136.7, 160,
   186.9, 209.9, 228.1, 241.2, 251.5, 259.5, 265.7, 270.9, 274.8, 278,
   280.3, 281.8, 283, 283.3, 283.7, 283.8, 283, 282.2, 281.2, 281.4,
   281.7, 281.1, 281.2}
  };

  double doy, p0, p1, pt;

  int imon, ilat;

  /* Get day of year... */
  doy = fmod(t / 86400., 365.25);
  while (doy < 0)
    doy += 365.25;

  /* Get indices... */
  imon = locate(doys, 12, doy);
  ilat = locate(lats, 73, lat);

  /* Get tropopause pressure... */
  p0 = LIN(lats[ilat], tps[imon][ilat],
	   lats[ilat + 1], tps[imon][ilat + 1], lat);
  p1 = LIN(lats[ilat], tps[imon + 1][ilat],
	   lats[ilat + 1], tps[imon + 1][ilat + 1], lat);
  pt = LIN(doys[imon], p0, doys[imon + 1], p1, doy);

  /* Return tropopause pressure... */
  return pt;
}

/*****************************************************************************/

void write_atm(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl,
  double t) {

  FILE *in, *out;

  char line[LEN];

  double r;

  int ip, iq, year, mon, day, hour, min, sec;

  /* Check if gnuplot output is requested... */
  if (ctl->atm_gpfile[0] != '-') {

    /* Write info... */
    printf("Plot atmospheric data: %s.png\n", filename);

    /* Create gnuplot pipe... */
    if (!(out = popen("gnuplot", "w")))
      ERRMSG("Cannot create pipe to gnuplot!");

    /* Set plot filename... */
    fprintf(out, "set out \"%s.png\"\n", filename);

    /* Set time string... */
    jsec2time(t, &year, &mon, &day, &hour, &min, &sec, &r);
    fprintf(out, "timestr=\"%d-%02d-%02d, %02d:%02d UTC\"\n",
	    year, mon, day, hour, min);

    /* Dump gnuplot file to pipe... */
    if (!(in = fopen(ctl->atm_gpfile, "r")))
      ERRMSG("Cannot open file!");
    while (fgets(line, LEN, in))
      fprintf(out, "%s", line);
    fclose(in);
  }

  else {

    /* Write info... */
    printf("Write atmospheric data: %s\n", filename);

    /* Create file... */
    if (!(out = fopen(filename, "w")))
      ERRMSG("Cannot create file!");
  }

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

void write_csi(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl,
  double t) {

  static FILE *in, *out;

  static char line[LEN];

  static double modmean[GX][GY][GZ], obsmean[GX][GY][GZ],
    rt, rz, rlon, rlat, t0, t1, area, dlon, dlat, lat;

  static int init, obscount[GX][GY][GZ], cx, cy, cz, ip, ix, iy, iz, robs;

  /* Init... */
  if (!init) {
    init = 1;

    /* Check quantity index for mass... */
    if (ctl->qnt_m < 0)
      ERRMSG("Need quantity mass to analyze CSI!");

    /* Open observation data file... */
    printf("Read CSI observation data: %s\n", ctl->csi_obsfile);
    if (!(in = fopen(ctl->csi_obsfile, "r")))
      ERRMSG("Cannot open observation data file!");

    /* Create new file... */
    printf("Write CSI data: %s\n", filename);
    if (!(out = fopen(filename, "w")))
      ERRMSG("Cannot create model data file!");

    /* Write header... */
    fprintf(out,
	    "# $1 = time [s]\n"
	    "# $2 = number of hits (cx)\n"
	    "# $3 = number of misses (cy)\n"
	    "# $4 = number of false alarms (cz)\n"
	    "# $5 = number of observations (cx + cy)\n"
	    "# $6 = number of forecasts (cx + cz)\n"
	    "# $7 = bias (forecasts/observations) [%%]\n"
	    "# $8 = probability of detection (POD) [%%]\n"
	    "# $9 = false alarm rate (FAR) [%%]\n"
	    "# $10 = critical success index (CSI) [%%]\n\n");
  }

  /* Set time interval... */
  t0 = t - 0.5 * ctl->dt_mod;
  t1 = t + 0.5 * ctl->dt_mod;

  /* Initialize grid cells... */
  for (ix = 0; ix < ctl->csi_nx; ix++)
    for (iy = 0; iy < ctl->csi_ny; iy++)
      for (iz = 0; iz < ctl->csi_nz; iz++)
	modmean[ix][iy][iz] = obsmean[ix][iy][iz] = obscount[ix][iy][iz] = 0;

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

    /* Calculate indices... */
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
    obscount[ix][iy][iz]++;
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
    modmean[ix][iy][iz] += atm->q[ctl->qnt_m][ip];
  }

  /* Analyze all grid cells... */
  for (ix = 0; ix < ctl->csi_nx; ix++)
    for (iy = 0; iy < ctl->csi_ny; iy++)
      for (iz = 0; iz < ctl->csi_nz; iz++) {

	/* Calculate mean observation index... */
	if (obscount[ix][iy][iz] > 0)
	  obsmean[ix][iy][iz] /= obscount[ix][iy][iz];

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
	if (obscount[ix][iy][iz] > 0) {
	  if (obsmean[ix][iy][iz] >= ctl->csi_obsmin &&
	      modmean[ix][iy][iz] >= ctl->csi_modmin)
	    cx++;
	  else if (obsmean[ix][iy][iz] >= ctl->csi_obsmin &&
		   modmean[ix][iy][iz] < ctl->csi_modmin)
	    cy++;
	  else if (obsmean[ix][iy][iz] < ctl->csi_obsmin &&
		   modmean[ix][iy][iz] >= ctl->csi_modmin)
	    cz++;
	}
      }

  /* Write output... */
  if (fmod(t, ctl->csi_dt_out) == 0) {

    /* Write... */
    fprintf(out, "%.2f %d %d %d %d %d %g %g %g %g\n",
	    t, cx, cy, cz, cx + cy, cx + cz,
	    (cx + cy > 0) ? 100. * (cx + cz) / (cx + cy) : GSL_NAN,
	    (cx + cy > 0) ? (100. * cx) / (cx + cy) : GSL_NAN,
	    (cx + cz > 0) ? (100. * cz) / (cx + cz) : GSL_NAN,
	    (cx + cy + cz > 0) ? (100. * cx) / (cx + cy + cz) : GSL_NAN);

    /* Set counters to zero... */
    cx = cy = cz = 0;
  }

  /* Close file... */
  if (t == ctl->t_stop)
    fclose(out);
}

/*****************************************************************************/

void write_grid(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl,
  double t) {

  FILE *in, *out;

  char line[LEN];

  static double grid_q[NQ][GX][GY][GZ], z, dz, lon, dlon, lat, dlat,
    area, rho_air, cd, mmr, t0, t1, r;

  static int grid_n[GX][GY][GZ], ip, ix, iy, iz, iq,
    year, mon, day, hour, min, sec;

  /* Check dimensions... */
  if (ctl->grid_nx > GX || ctl->grid_ny > GY || ctl->grid_nz > GZ)
    ERRMSG("Grid dimensions too large!");

  /* Set time interval for output... */
  t0 = t - 0.5 * ctl->dt_mod;
  t1 = t + 0.5 * ctl->dt_mod;

  /* Set grid box size... */
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

  /* Check if gnuplot output is requested... */
  if (ctl->grid_gpfile[0] != '-') {

    /* Write info... */
    printf("Plot gridded data: %s.png\n", filename);

    /* Create gnuplot pipe... */
    if (!(out = popen("gnuplot", "w")))
      ERRMSG("Cannot create pipe to gnuplot!");

    /* Set plot filename... */
    fprintf(out, "set out \"%s.png\"\n", filename);

    /* Set time string... */
    jsec2time(t, &year, &mon, &day, &hour, &min, &sec, &r);
    fprintf(out, "timestr=\"%d-%02d-%02d, %02d:%02d UTC\"\n",
	    year, mon, day, hour, min);

    /* Dump gnuplot file to pipe... */
    if (!(in = fopen(ctl->grid_gpfile, "r")))
      ERRMSG("Cannot open file!");
    while (fgets(line, LEN, in))
      fprintf(out, "%s", line);
    fclose(in);
  }

  else {

    /* Write info... */
    printf("Write gridded data: %s\n", filename);

    /* Create file... */
    if (!(out = fopen(filename, "w")))
      ERRMSG("Cannot create file!");
  }

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
	if (ctl->qnt_m >= 0)
	  cd =
	    grid_q[ctl->qnt_m][ix][iy][iz] * grid_n[ix][iy][iz] / (1e6 *
								   area);

	/* Calculate mass mixing ratio... */
	if (ctl->qnt_m >= 0 && ctl->qnt_t >= 0) {
	  rho_air = 100. * P(z) / (287.058 * grid_q[ctl->qnt_t][ix][iy][iz]);
	  mmr = grid_q[ctl->qnt_m][ix][iy][iz] * grid_n[ix][iy][iz]
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

void write_sample(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl,
  double t) {

  static FILE *in, *out;

  static char line[LEN];

  static double area, cd, mmr, rho_air, rt, rz, rlon, rlat, rdz, rdx,
    dlat, dx2, p0, p1, t0, t1, x0[3], x1[3], q_mean[NQ];

  static int init, ip, iq, np, ridx, ridx_old = -999;

  /* Open files... */
  if (!init) {
    init = 1;

    /* Open observation data file... */
    printf("Read sample observation data: %s\n", ctl->sample_obsfile);
    if (!(in = fopen(ctl->sample_obsfile, "r")))
      ERRMSG("Cannot open observation data file!");

    /* Create new file... */
    printf("Write sample model data: %s\n", filename);
    if (!(out = fopen(filename, "w")))
      ERRMSG("Cannot create model data file!");

    /* Write header... */
    fprintf(out,
	    "# $1 = time [s]\n"
	    "# $2 = altitude [km]\n"
	    "# $3 = longitude [deg]\n"
	    "# $4 = latitude [deg]\n"
	    "# $5 = vertical layer width [km]\n"
	    "# $6 = search radius [km]\n"
	    "# $7 = data set index\n"
	    "# $8 = number of particles\n"
	    "# $9 = column density [kg/m^2]\n"
	    "# $10 = mass mixing ratio [kg/kg]\n");
    for (iq = 0; iq < ctl->nq; iq++)
      fprintf(out, "# $%i = %s (mean) [%s]\n", iq + 11, ctl->qnt_name[iq],
	      ctl->qnt_unit[iq]);
  }

  /* Set time interval... */
  t0 = t - 0.5 * ctl->dt_mod;
  t1 = t + 0.5 * ctl->dt_mod;

  /* Read data... */
  while (fgets(line, LEN, in)) {

    /* Read data... */
    if (sscanf(line, "%lg %lg %lg %lg %lg %lg %d",
	       &rt, &rz, &rlon, &rlat, &rdz, &rdx, &ridx) != 7)
      continue;

    /* Check time... */
    if (rt < t0)
      continue;
    if (rt > t1)
      break;

    /* Set search ranges... */
    p0 = P(rz - 0.5 * rdz);
    p1 = P(rz + 0.5 * rdz);
    dlat = dy2deg(rdx);
    dx2 = gsl_pow_2(rdx);

    /* Get geolocation... */
    geo2cart(0, rlon, rlat, x0);

    /* Init... */
    np = 0;
    for (iq = 0; iq < ctl->nq; iq++)
      q_mean[iq] = 0;

    /* Loop over air parcles... */
    for (ip = 0; ip < atm->np; ip++) {

      /* Check latitudinal distance... */
      if (fabs(atm->lat[ip] - rlat) > dlat)
	continue;

      /* Check vertical distance... */
      if (atm->p[ip] > p0 || atm->p[ip] < p1)
	continue;

      /* Check horizontal distance... */
      geo2cart(0, atm->lon[ip], atm->lat[ip], x1);
      if (DIST2(x0, x1) > dx2)
	continue;

      /* Calculate mean values... */
      for (iq = 0; iq < ctl->nq; iq++)
	q_mean[iq] += atm->q[iq][ip];
      np++;
    }

    /* Calculate mean values... */
    if (np > 0)
      for (iq = 0; iq < ctl->nq; iq++)
	q_mean[iq] /= np;

    /* Calculate surface area... */
    area = M_PI * dx2;

    /* Calculate column density... */
    if (ctl->qnt_m >= 0)
      cd = q_mean[ctl->qnt_m] * np / (1e6 * area);

    /* Calculate mass mixing ratio... */
    if (ctl->qnt_m >= 0 && ctl->qnt_t >= 0) {
      rho_air = 100. * P(rz) / (287.058 * q_mean[ctl->qnt_t]);
      mmr = q_mean[ctl->qnt_m] * np / (rho_air * 1e6 * area * 1e3 * rdz);
    }

    /* Write output... */
    if (ridx != ridx_old)
      fprintf(out, "\n");
    ridx_old = ridx;
    fprintf(out, "%.2f %g %g %g %g %g %d %d %g %g",
	    rt, rz, rlon, rlat, rdz, rdx, ridx, np, cd, mmr);
    for (iq = 0; iq < ctl->nq; iq++) {
      fprintf(out, " ");
      fprintf(out, ctl->qnt_format[iq], q_mean[iq]);
    }
    fprintf(out, "\n");
  }

  /* Close files... */
  if (t == ctl->t_stop) {
    fclose(out);
    fclose(in);
  }
}

/*****************************************************************************/

void write_station(
  const char *filename,
  atm_t * atm,
  ctl_t * ctl,
  double t) {

  static FILE *out;

  static double rmax2, t0, t1, x0[3], x1[3];

  static int init, ip, iq;

  /* Init... */
  if (!init) {
    init = 1;

    /* Write info... */
    printf("Write station data: %s\n", filename);

    /* Create new file... */
    if (!(out = fopen(filename, "w")))
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

    /* Set geolocation and search radius... */
    geo2cart(0, ctl->stat_lon, ctl->stat_lat, x0);
    rmax2 = gsl_pow_2(ctl->stat_r);
  }

  /* Set time interval for output... */
  t0 = t - 0.5 * ctl->dt_mod;
  t1 = t + 0.5 * ctl->dt_mod;

  /* Loop over air parcels... */
  for (ip = 0; ip < atm->np; ip++) {

    /* Check time... */
    if (atm->time[ip] < t0 || atm->time[ip] > t1)
      continue;

    /* Check station flag... */
    if (ctl->qnt_stat >= 0)
      if (atm->q[ctl->qnt_stat][ip])
	continue;

    /* Get Cartesian coordinates... */
    geo2cart(0, atm->lon[ip], atm->lat[ip], x1);

    /* Check horizontal distance... */
    if (DIST2(x0, x1) > rmax2)
      continue;

    /* Set station flag... */
    if (ctl->qnt_stat >= 0)
      atm->q[ctl->qnt_stat][ip] = 1;

    /* Write data... */
    fprintf(out, "%.2f %g %g %g",
	    atm->time[ip], Z(atm->p[ip]), atm->lon[ip], atm->lat[ip]);
    for (iq = 0; iq < ctl->nq; iq++) {
      fprintf(out, " ");
      fprintf(out, ctl->qnt_format[iq], atm->q[iq][ip]);
    }
    fprintf(out, "\n");
  }

  /* Close file... */
  if (t == ctl->t_stop)
    fclose(out);
}
