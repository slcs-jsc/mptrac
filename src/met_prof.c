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
  
  Copyright (C) 2013-2018 Forschungszentrum Juelich GmbH
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
  int const argc,
  char const *argv[]) {

  ctl_t ctl;

  met_t *met;

  FILE *out;

  static double timem[NZ], z, z0, z1, dz, lon, lon0, lon1, dlon, lonm[NZ],
    lat, lat0, lat1, dlat, latm[NZ], t, tm[NZ], uvw[3], um[NZ], vm[NZ],
    wm[NZ], h2o, h2om[NZ], o3, o3m[NZ], ps, psm[NZ], pt, ptm[NZ], tt, ttm[NZ],
    zg, zgm[NZ], zt, ztm[NZ], pv, pvm[NZ];

  static int i, iz, np[NZ];

  /* Allocate... */
  ALLOC(met, met_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <prof.tab> <met0> [ <met1> ... ]");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  z0 = scan_ctl(argv[1], argc, argv, "Z0", -1, "0", NULL);
  z1 = scan_ctl(argv[1], argc, argv, "Z1", -1, "60", NULL);
  dz = scan_ctl(argv[1], argc, argv, "DZ", -1, "1", NULL);
  lon0 = scan_ctl(argv[1], argc, argv, "LON0", -1, "0", NULL);
  lon1 = scan_ctl(argv[1], argc, argv, "LON1", -1, "0", NULL);
  dlon = scan_ctl(argv[1], argc, argv, "DLON", -1, "1", NULL);
  lat0 = scan_ctl(argv[1], argc, argv, "LAT0", -1, "0", NULL);
  lat1 = scan_ctl(argv[1], argc, argv, "LAT1", -1, "0", NULL);
  dlat = scan_ctl(argv[1], argc, argv, "DLAT", -1, "1", NULL);

  /* Loop over input files... */
  for (i = 3; i < argc; i++) {

    /* Read meteorological data... */
    read_met(&ctl, argv[i], met);

    /* Average... */
    for (z = z0; z <= z1; z += dz) {
      iz = (int) ((z - z0) / dz);
      if (iz < 0 || iz > NZ)
	ERRMSG("Too many altitudes!");
      for (lon = lon0; lon <= lon1; lon += dlon)
	for (lat = lat0; lat <= lat1; lat += dlat) {
	  intpol_met_space(met, P(z), lon, lat, 
          &ps, &pt, &zg, &t, &pv, &h2o, &o3);
      intpol_winds_space(met, P(z), lon, lat, uvw);
	  intpol_met_space(met, pt, lon, lat, 
          NULL, NULL, &zt, &tt, NULL, NULL, NULL);
	  if (gsl_finite(t) && gsl_finite(uvw[0]) && gsl_finite(uvw[1]) && gsl_finite(uvw[2])) {
	    timem[iz] += met->time;
	    lonm[iz] += lon;
	    latm[iz] += lat;
	    zgm[iz] += zg;
	    tm[iz] += t;
	    um[iz] += uvw[0];
	    vm[iz] += uvw[1];
	    wm[iz] += uvw[2];
	    pvm[iz] += pv;
	    h2om[iz] += h2o;
	    o3m[iz] += o3;
	    psm[iz] += ps;
	    ptm[iz] += pt;
	    ztm[iz] += zt;
	    ttm[iz] += tt;
	    np[iz]++;
	  }
	}
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
	  "# $9  = vertical wind [hPa/s]\n");
  fprintf(out,
	  "# $10 = H2O volume mixing ratio [1]\n"
	  "# $11 = O3 volume mixing ratio [1]\n"
	  "# $12 = geopotential height [km]\n"
	  "# $13 = potential vorticity [PVU]\n"
	  "# $14 = surface pressure [hPa]\n"
	  "# $15 = tropopause pressure [hPa]\n"
	  "# $16 = tropopause geopotential height [km]\n"
	  "# $17 = tropopause temperature [K]\n\n");

  /* Write data... */
  for (z = z0; z <= z1; z += dz) {
    iz = (int) ((z - z0) / dz);
    fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	    timem[iz] / np[iz], z, lonm[iz] / np[iz], latm[iz] / np[iz], P(z),
	    tm[iz] / np[iz], um[iz] / np[iz], vm[iz] / np[iz],
	    wm[iz] / np[iz], h2om[iz] / np[iz], o3m[iz] / np[iz],
	    zgm[iz] / np[iz], pvm[iz] / np[iz], psm[iz] / np[iz],
	    ptm[iz] / np[iz], ztm[iz] / np[iz], ttm[iz] / np[iz]);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(met);

  return EXIT_SUCCESS;
}
