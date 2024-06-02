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
  
  Copyright (C) 2013-2023 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Create tropopause data set from meteorological data.
*/

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  clim_t *clim;

  met_t *met;

  static double ps[EX * EY], pt[EX * EY], qt[EX * EY], o3t[EX * EY],
    zs[EX * EY], zt[EX * EY], tt[EX * EY], lon, lon0, lon1, lons[EX], dlon,
    lat, lat0, lat1, lats[EY], dlat;

  static int init, i, nx, ny, nt, ncid, varid, dims[3], h2o, o3;

  static size_t count[10], start[10];

  /* Allocate... */
  ALLOC(clim, clim_t, 1);
  ALLOC(met, met_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <tropo.nc> <met0> [ <met1> ... ]");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  lon0 = scan_ctl(argv[1], argc, argv, "TROPO_LON0", -1, "-180", NULL);
  lon1 = scan_ctl(argv[1], argc, argv, "TROPO_LON1", -1, "180", NULL);
  dlon = scan_ctl(argv[1], argc, argv, "TROPO_DLON", -1, "-999", NULL);
  lat0 = scan_ctl(argv[1], argc, argv, "TROPO_LAT0", -1, "-90", NULL);
  lat1 = scan_ctl(argv[1], argc, argv, "TROPO_LAT1", -1, "90", NULL);
  dlat = scan_ctl(argv[1], argc, argv, "TROPO_DLAT", -1, "-999", NULL);
  h2o = (int) scan_ctl(argv[1], argc, argv, "TROPO_H2O", -1, "1", NULL);
  o3 = (int) scan_ctl(argv[1], argc, argv, "TROPO_O3", -1, "1", NULL);

  /* Read climatological data... */
  read_clim(&ctl, clim);

  /* Loop over files... */
  for (i = 3; i < argc; i++) {

    /* Set control parameters... */
    ctl.met_tropo = 0;

    /* Read meteorological data... */
    if (!read_met(argv[i], &ctl, clim, met))
      continue;

    /* Set horizontal grid... */
    if (!init) {
      init = 1;

      /* Get grid... */
      if (dlon <= 0)
	dlon = fabs(met->lon[1] - met->lon[0]);
      if (dlat <= 0)
	dlat = fabs(met->lat[1] - met->lat[0]);
      if (lon0 < -360 && lon1 > 360) {
	lon0 = gsl_stats_min(met->lon, 1, (size_t) met->nx);
	lon1 = gsl_stats_max(met->lon, 1, (size_t) met->nx);
      }
      nx = ny = 0;
      for (lon = lon0; lon <= lon1; lon += dlon) {
	lons[nx] = lon;
	if ((++nx) > EX)
	  ERRMSG("Too many longitudes!");
      }
      if (lat0 < -90 && lat1 > 90) {
	lat0 = gsl_stats_min(met->lat, 1, (size_t) met->ny);
	lat1 = gsl_stats_max(met->lat, 1, (size_t) met->ny);
      }
      for (lat = lat0; lat <= lat1; lat += dlat) {
	lats[ny] = lat;
	if ((++ny) > EY)
	  ERRMSG("Too many latitudes!");
      }

      /* Create netCDF file... */
      LOG(1, "Write tropopause data file: %s", argv[2]);
      NC(nc_create(argv[2], NC_CLOBBER, &ncid));

      /* Create dimensions... */
      NC(nc_def_dim(ncid, "time", (size_t) NC_UNLIMITED, &dims[0]));
      NC(nc_def_dim(ncid, "lat", (size_t) ny, &dims[1]));
      NC(nc_def_dim(ncid, "lon", (size_t) nx, &dims[2]));

      /* Create variables... */
      NC_DEF_VAR("time", NC_DOUBLE, 1, &dims[0], "time",
		 "seconds since 2000-01-01 00:00:00 UTC");
      NC_DEF_VAR("lat", NC_DOUBLE, 1, &dims[1], "latitude", "degrees_north");
      NC_DEF_VAR("lon", NC_DOUBLE, 1, &dims[2], "longitude", "degrees_east");

      NC_DEF_VAR("clp_z", NC_FLOAT, 3, &dims[0], "cold point height", "km");
      NC_DEF_VAR("clp_p", NC_FLOAT, 3, &dims[0], "cold point pressure",
		 "hPa");
      NC_DEF_VAR("clp_t", NC_FLOAT, 3, &dims[0], "cold point temperature",
		 "K");
      if (h2o)
	NC_DEF_VAR("clp_q", NC_FLOAT, 3, &dims[0], "cold point water vapor",
		   "ppv");
      if (o3)
	NC_DEF_VAR("clp_o3", NC_FLOAT, 3, &dims[0], "cold point ozone",
		   "ppv");

      NC_DEF_VAR("dyn_z", NC_FLOAT, 3, &dims[0],
		 "dynamical tropopause height", "km");
      NC_DEF_VAR("dyn_p", NC_FLOAT, 3, &dims[0],
		 "dynamical tropopause pressure", "hPa");
      NC_DEF_VAR("dyn_t", NC_FLOAT, 3, &dims[0],
		 "dynamical tropopause temperature", "K");
      if (h2o)
	NC_DEF_VAR("dyn_q", NC_FLOAT, 3, &dims[0],
		   "dynamical tropopause water vapor", "ppv");
      if (o3)
	NC_DEF_VAR("dyn_o3", NC_FLOAT, 3, &dims[0],
		   "dynamical tropopause ozone", "ppv");

      NC_DEF_VAR("wmo_1st_z", NC_FLOAT, 3, &dims[0],
		 "WMO 1st tropopause height", "km");
      NC_DEF_VAR("wmo_1st_p", NC_FLOAT, 3, &dims[0],
		 "WMO 1st tropopause pressure", "hPa");
      NC_DEF_VAR("wmo_1st_t", NC_FLOAT, 3, &dims[0],
		 "WMO 1st tropopause temperature", "K");
      if (h2o)
	NC_DEF_VAR("wmo_1st_q", NC_FLOAT, 3, &dims[0],
		   "WMO 1st tropopause water vapor", "ppv");
      if (o3)
	NC_DEF_VAR("wmo_1st_o3", NC_FLOAT, 3, &dims[0],
		   "WMO 1st tropopause ozone", "ppv");

      NC_DEF_VAR("wmo_2nd_z", NC_FLOAT, 3, &dims[0],
		 "WMO 2nd tropopause height", "km");
      NC_DEF_VAR("wmo_2nd_p", NC_FLOAT, 3, &dims[0],
		 "WMO 2nd tropopause pressure", "hPa");
      NC_DEF_VAR("wmo_2nd_t", NC_FLOAT, 3, &dims[0],
		 "WMO 2nd tropopause temperature", "K");
      if (h2o)
	NC_DEF_VAR("wmo_2nd_q", NC_FLOAT, 3, &dims[0],
		   "WMO 2nd tropopause water vapor", "ppv");
      if (o3)
	NC_DEF_VAR("wmo_2nd_o3", NC_FLOAT, 3, &dims[0],
		   "WMO 2nd tropopause ozone", "ppv");

      NC_DEF_VAR("ps", NC_FLOAT, 3, &dims[0], "surface pressure", "hPa");
      NC_DEF_VAR("zs", NC_FLOAT, 3, &dims[0], "surface height", "km");

      /* End definition... */
      NC(nc_enddef(ncid));

      /* Write longitude and latitude... */
      NC_PUT_DOUBLE("lat", lats, 0);
      NC_PUT_DOUBLE("lon", lons, 0);
    }

    /* Write time... */
    start[0] = (size_t) nt;
    count[0] = 1;
    start[1] = 0;
    count[1] = (size_t) ny;
    start[2] = 0;
    count[2] = (size_t) nx;
    NC_PUT_DOUBLE("time", &met->time, 1);

    /* Get cold point... */
    get_tropo(2, &ctl, clim, met, lons, nx, lats, ny, pt, zt, tt, qt, o3t, ps,
	      zs);
    NC_PUT_DOUBLE("clp_z", zt, 1);
    NC_PUT_DOUBLE("clp_p", pt, 1);
    NC_PUT_DOUBLE("clp_t", tt, 1);
    if (h2o)
      NC_PUT_DOUBLE("clp_q", qt, 1);
    if (o3)
      NC_PUT_DOUBLE("clp_o3", o3t, 1);

    /* Get dynamical tropopause... */
    get_tropo(5, &ctl, clim, met, lons, nx, lats, ny, pt, zt, tt, qt, o3t, ps,
	      zs);
    NC_PUT_DOUBLE("dyn_z", zt, 1);
    NC_PUT_DOUBLE("dyn_p", pt, 1);
    NC_PUT_DOUBLE("dyn_t", tt, 1);
    if (h2o)
      NC_PUT_DOUBLE("dyn_q", qt, 1);
    if (o3)
      NC_PUT_DOUBLE("dyn_o3", o3t, 1);

    /* Get WMO 1st tropopause... */
    get_tropo(3, &ctl, clim, met, lons, nx, lats, ny, pt, zt, tt, qt, o3t, ps,
	      zs);
    NC_PUT_DOUBLE("wmo_1st_z", zt, 1);
    NC_PUT_DOUBLE("wmo_1st_p", pt, 1);
    NC_PUT_DOUBLE("wmo_1st_t", tt, 1);
    if (h2o)
      NC_PUT_DOUBLE("wmo_1st_q", qt, 1);
    if (o3)
      NC_PUT_DOUBLE("wmo_1st_o3", o3t, 1);

    /* Get WMO 2nd tropopause... */
    get_tropo(4, &ctl, clim, met, lons, nx, lats, ny, pt, zt, tt, qt, o3t, ps,
	      zs);
    NC_PUT_DOUBLE("wmo_2nd_z", zt, 1);
    NC_PUT_DOUBLE("wmo_2nd_p", pt, 1);
    NC_PUT_DOUBLE("wmo_2nd_t", tt, 1);
    if (h2o)
      NC_PUT_DOUBLE("wmo_2nd_q", qt, 1);
    if (o3)
      NC_PUT_DOUBLE("wmo_2nd_o3", o3t, 1);

    /* Write surface data... */
    NC_PUT_DOUBLE("ps", ps, 1);
    NC_PUT_DOUBLE("zs", zs, 1);

    /* Increment time step counter... */
    nt++;
  }

  /* Close file... */
  NC(nc_close(ncid));

  /* Free... */
  free(clim);
  free(met);

  return EXIT_SUCCESS;
}
