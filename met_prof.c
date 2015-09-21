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
  Extract vertical profile from meteorological data.
*/

#include "libtrac.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum number of altitudes. */
#define NZ 1000

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  met_t *met;

  FILE *in, *out;

  static double timem[NZ], z, z0, z1, dz, lon, lon0, lon1, dlon, lonm[NZ],
    lat, lat0, lat1, dlat, latm[NZ], t, tm[NZ],
    u, um[NZ], v, vm[NZ], w, wm[NZ];

  static int epol, i, iz, np[NZ], redx, redy, redp;

  /* Allocate... */
  ALLOC(met, met_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <prof.tab> <met0> [ <met1> ... ]");

  /* Read control parameters... */
  read_ctl(NULL, NULL, argc, argv, &ctl);
  z0 = scan_ctl(NULL, NULL, argc, argv, "Z0", -1, "0", NULL);
  z1 = scan_ctl(NULL, NULL, argc, argv, "Z1", -1, "60", NULL);
  dz = scan_ctl(NULL, NULL, argc, argv, "DZ", -1, "1", NULL);
  lon0 = scan_ctl(NULL, NULL, argc, argv, "LON0", -1, "0", NULL);
  lon1 = scan_ctl(NULL, NULL, argc, argv, "LON1", -1, "0", NULL);
  dlon = scan_ctl(NULL, NULL, argc, argv, "DLON", -1, "1", NULL);
  lat0 = scan_ctl(NULL, NULL, argc, argv, "LAT0", -1, "0", NULL);
  lat1 = scan_ctl(NULL, NULL, argc, argv, "LAT1", -1, "0", NULL);
  dlat = scan_ctl(NULL, NULL, argc, argv, "DLAT", -1, "1", NULL);
  epol = (int) scan_ctl(NULL, NULL, argc, argv, "EPOL", -1, "0", NULL);
  redx = (int) scan_ctl(NULL, NULL, argc, argv, "REDX", -1, "1", NULL);
  redy = (int) scan_ctl(NULL, NULL, argc, argv, "REDY", -1, "1", NULL);
  redp = (int) scan_ctl(NULL, NULL, argc, argv, "REDP", -1, "1", NULL);

  /* Loop over input files... */
  for (i = 3; i < argc; i++) {

    /* Read meteorological data... */
    if (!(in = fopen(argv[i], "r")))
      continue;
    else
      fclose(in);
    read_met(argv[i], met);

    /* Extrapolate... */
    if (epol)
      extrapolate_met(met);

    /* Reduce resolution... */
    reduce_met(met, redx, redy, redp);

    /* Average... */
    for (z = z0; z <= z1; z += dz) {
      iz = (int) ((z - z0) / dz);
      if (iz < 0 || iz > NZ)
	ERRMSG("Too many altitudes!");
      for (lon = lon0; lon <= lon1; lon += dlon)
	for (lat = lat0; lat <= lat1; lat += dlat) {
	  intpol_met_space(met, P(z), lon, lat, &t, &u, &v, &w);
	  if (gsl_finite(t) && gsl_finite(u)
	      && gsl_finite(v) && gsl_finite(w)) {
	    timem[iz] += met->time;
	    lonm[iz] += lon;
	    latm[iz] += lat;
	    tm[iz] += t;
	    um[iz] += u;
	    vm[iz] += v;
	    wm[iz] += w;
	    np[iz]++;
	  }
	}
    }
  }

  /* Normalize... */
  for (z = z0; z <= z1; z += dz) {
    iz = (int) ((z - z0) / dz);
    if (np[iz] > 0) {
      timem[iz] /= np[iz];
      lonm[iz] /= np[iz];
      latm[iz] /= np[iz];
      tm[iz] /= np[iz];
      um[iz] /= np[iz];
      vm[iz] /= np[iz];
      wm[iz] /= np[iz];
    } else {
      timem[iz] = GSL_NAN;
      lonm[iz] = GSL_NAN;
      latm[iz] = GSL_NAN;
      tm[iz] = GSL_NAN;
      um[iz] = GSL_NAN;
      vm[iz] = GSL_NAN;
      wm[iz] = GSL_NAN;
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
	  "# $9  = vertical wind [hPa/s]\n\n");

  /* Write data... */
  for (z = z0; z <= z1; z += dz) {
    iz = (int) ((z - z0) / dz);
    fprintf(out, "%.2f %g %g %g %g %g %g %g %g\n",
	    timem[iz], z, lonm[iz], latm[iz], P(z),
	    tm[iz], um[iz], vm[iz], wm[iz]);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(met);

  return EXIT_SUCCESS;
}
