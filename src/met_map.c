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
  
  Copyright (C) 2013-2019 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Extract global map from meteorological data.
*/

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  met_t *met;

  FILE *out;

  static double timem[EX][EY], p, ps, psm[EX][EY], pt, ptm[EX][EY], t,
    tm[EX][EY], u, um[EX][EY], v, vm[EX][EY], w, wm[EX][EY], h2o,
    h2om[EX][EY], o3, o3m[EX][EY], z, zm[EX][EY], pv, pvm[EX][EY], zt,
    ztm[EX][EY], tt, ttm[EX][EY], lon, lon0, lon1, lons[EX], dlon, lat, lat0,
    lat1, lats[EY], dlat;

  static int i, ix, iy, np[EX][EY], nx, ny;

  /* Allocate... */
  ALLOC(met, met_t, 1);

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <map.tab> <z> <met0> [ <met1> ... ]");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  lon0 = scan_ctl(argv[1], argc, argv, "LON0", -1, "-180", NULL);
  lon1 = scan_ctl(argv[1], argc, argv, "LON1", -1, "180", NULL);
  dlon = scan_ctl(argv[1], argc, argv, "DLON", -1, "-999", NULL);
  lat0 = scan_ctl(argv[1], argc, argv, "LAT0", -1, "-90", NULL);
  lat1 = scan_ctl(argv[1], argc, argv, "LAT1", -1, "90", NULL);
  dlat = scan_ctl(argv[1], argc, argv, "DLAT", -1, "-999", NULL);
  p = P(atof(argv[3]));

  /* Loop over files... */
  for (i = 4; i < argc; i++) {

    /* Read meteorological data... */
    read_met(&ctl, argv[i], met);

    /* Set horizontal grid... */
    if (dlon <= 0)
      dlon = fabs(met->lon[1] - met->lon[0]);
    if (dlat <= 0)
      dlat = fabs(met->lat[1] - met->lat[0]);
    if (lon0 < -360 && lon1 > 360) {
      lon0 = gsl_stats_min(met->lon, 1, (size_t) met->nx);
      lon1 = gsl_stats_max(met->lon, 1, (size_t) met->nx);
    }
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

    /* Average data... */
    for (ix = 0; ix < nx; ix++)
      for (iy = 0; iy < ny; iy++) {

	/* Interpolate to given log-pressure height... */
	intpol_met_space(met, p, lons[ix], lats[iy], &ps, &pt,
			 &z, &t, &u, &v, &w, &pv, &h2o, &o3);

	/* Get tropopause data... */
	intpol_met_space(met, pt, lons[ix], lats[iy], NULL, NULL,
			 &zt, &tt, NULL, NULL, NULL, NULL, NULL, NULL);

	/* Sum up data... */
	timem[ix][iy] += met->time;
	zm[ix][iy] += z;
	tm[ix][iy] += t;
	um[ix][iy] += u;
	vm[ix][iy] += v;
	wm[ix][iy] += w;
	pvm[ix][iy] += pv;
	h2om[ix][iy] += h2o;
	o3m[ix][iy] += o3;
	psm[ix][iy] += ps;
	ptm[ix][iy] += pt;
	ztm[ix][iy] += zt;
	ttm[ix][iy] += tt;
	np[ix][iy]++;
      }
  }

  /* Create output file... */
  printf("Write meteorological data file: %s\n", argv[2]);
  if (!(out = fopen(argv[2], "w")))
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
	  "# $9  = vertical wind [hPa/s]\n"
	  "# $10 = H2O volume mixing ratio [1]\n");
  fprintf(out,
	  "# $11 = O3 volume mixing ratio [1]\n"
	  "# $12 = geopotential height [km]\n"
	  "# $13 = potential vorticity [PVU]\n"
	  "# $14 = surface pressure [hPa]\n"
	  "# $15 = tropopause pressure [hPa]\n"
	  "# $16 = tropopause geopotential height [km]\n"
	  "# $17 = tropopause temperature [K]\n");

  /* Write data... */
  for (iy = 0; iy < ny; iy++) {
    fprintf(out, "\n");
    for (ix = 0; ix < nx; ix++)
      fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	      timem[ix][iy] / np[ix][iy], Z(p), lons[ix], lats[iy], p,
	      tm[ix][iy] / np[ix][iy], um[ix][iy] / np[ix][iy],
	      vm[ix][iy] / np[ix][iy], wm[ix][iy] / np[ix][iy],
	      h2om[ix][iy] / np[ix][iy], o3m[ix][iy] / np[ix][iy],
	      zm[ix][iy] / np[ix][iy], pvm[ix][iy] / np[ix][iy],
	      psm[ix][iy] / np[ix][iy], ptm[ix][iy] / np[ix][iy],
	      ztm[ix][iy] / np[ix][iy], ttm[ix][iy] / np[ix][iy]);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(met);

  return EXIT_SUCCESS;
}
