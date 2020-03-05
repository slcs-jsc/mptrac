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
  Extract zonal mean from meteorological data.
*/

#include "libtrac.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Maximum number of altitudes. */
#define NZ 1000

/*! Maximum number of latitudes. */
#define NY 721

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  met_t *met;

  FILE *out;

  static double timem[NZ][NY], psm[NZ][NY], ptm[NZ][NY], ttm[NZ][NY],
    ztm[NZ][NY], tm[NZ][NY], um[NZ][NY], vm[NZ][NY], wm[NZ][NY], h2om[NZ][NY],
    h2otm[NZ][NY], pvm[NZ][NY], o3m[NZ][NY], lwcm[NZ][NY], iwcm[NZ][NY],
    zm[NZ][NY], z, z0, z1, dz, zt, tt, plev[NZ], ps, pt, t, u, v, w, pv, h2o,
    h2ot, o3, lwc, iwc, lat, lat0, lat1, dlat, lats[NY], cw[3];

  static int i, ix, iy, iz, np[NZ][NY], npt[NZ][NY], ny, nz, ci[3];

  /* Allocate... */
  ALLOC(met, met_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <zm.tab> <met0> [ <met1> ... ]");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  z0 = scan_ctl(argv[1], argc, argv, "ZM_Z0", -1, "-999", NULL);
  z1 = scan_ctl(argv[1], argc, argv, "ZM_Z1", -1, "-999", NULL);
  dz = scan_ctl(argv[1], argc, argv, "ZM_DZ", -1, "-999", NULL);
  lat0 = scan_ctl(argv[1], argc, argv, "ZM_LAT0", -1, "-90", NULL);
  lat1 = scan_ctl(argv[1], argc, argv, "ZM_LAT1", -1, "90", NULL);
  dlat = scan_ctl(argv[1], argc, argv, "ZM_DLAT", -1, "-999", NULL);

  /* Loop over files... */
  for (i = 3; i < argc; i++) {

    /* Read meteorological data... */
    if (!read_met(&ctl, argv[i], met))
      continue;

    /* Set vertical grid... */
    if (z0 < 0)
      z0 = Z(met->p[0]);
    if (z1 < 0)
      z1 = Z(met->p[met->np - 1]);
    nz = 0;
    if (dz < 0) {
      for (iz = 0; iz < met->np; iz++)
	if (Z(met->p[iz]) >= z0 && Z(met->p[iz]) <= z1) {
	  plev[nz] = met->p[iz];
	  if ((++nz) > NZ)
	    ERRMSG("Too many pressure levels!");
	}
    } else
      for (z = z0; z <= z1; z += dz) {
	plev[nz] = P(z);
	if ((++nz) > NZ)
	  ERRMSG("Too many pressure levels!");
      }

    /* Set horizontal grid... */
    if (dlat <= 0)
      dlat = fabs(met->lat[1] - met->lat[0]);
    ny = 0;
    if (lat0 < -90 && lat1 > 90) {
      lat0 = gsl_stats_min(met->lat, 1, (size_t) met->ny);
      lat1 = gsl_stats_max(met->lat, 1, (size_t) met->ny);
    }
    for (lat = lat0; lat <= lat1; lat += dlat) {
      lats[ny] = lat;
      if ((++ny) > NY)
	ERRMSG("Too many latitudes!");
    }

    /* Average... */
    for (ix = 0; ix < met->nx; ix++)
      for (iy = 0; iy < ny; iy++)
	for (iz = 0; iz < nz; iz++) {

	  /* Interpolate meteo data... */
	  intpol_met_space_3d(met, met->z, plev[iz], met->lon[ix],
			      met->lat[iy], &z, ci, cw, 1);
	  intpol_met_space_3d(met, met->t, plev[iz], met->lon[ix],
			      met->lat[iy], &t, ci, cw, 0);
	  intpol_met_space_3d(met, met->u, plev[iz], met->lon[ix],
			      met->lat[iy], &u, ci, cw, 0);
	  intpol_met_space_3d(met, met->v, plev[iz], met->lon[ix],
			      met->lat[iy], &v, ci, cw, 0);
	  intpol_met_space_3d(met, met->w, plev[iz], met->lon[ix],
			      met->lat[iy], &w, ci, cw, 0);
	  intpol_met_space_3d(met, met->pv, plev[iz], met->lon[ix],
			      met->lat[iy], &pv, ci, cw, 0);
	  intpol_met_space_3d(met, met->h2o, plev[iz], met->lon[ix],
			      met->lat[iy], &h2o, ci, cw, 0);
	  intpol_met_space_3d(met, met->o3, plev[iz], met->lon[ix],
			      met->lat[iy], &o3, ci, cw, 0);
	  intpol_met_space_3d(met, met->lwc, plev[iz], met->lon[ix],
			      met->lat[iy], &lwc, ci, cw, 0);
	  intpol_met_space_3d(met, met->iwc, plev[iz], met->lon[ix],
			      met->lat[iy], &iwc, ci, cw, 0);
	  intpol_met_space_2d(met, met->ps, met->lon[ix], met->lat[iy], &ps,
			      ci, cw, 0);
	  intpol_met_space_2d(met, met->pt, met->lon[ix], met->lat[iy], &pt,
			      ci, cw, 0);

	  /* Interpolate tropopause data... */
	  intpol_met_space_3d(met, met->z, pt, met->lon[ix], met->lat[iy],
			      &zt, ci, cw, 1);
	  intpol_met_space_3d(met, met->t, pt, met->lon[ix], met->lat[iy],
			      &tt, ci, cw, 0);
	  intpol_met_space_3d(met, met->h2o, pt, met->lon[ix], met->lat[iy],
			      &h2ot, ci, cw, 0);

	  /* Averaging... */
	  timem[iz][iy] += met->time;
	  zm[iz][iy] += z;
	  tm[iz][iy] += t;
	  um[iz][iy] += u;
	  vm[iz][iy] += v;
	  wm[iz][iy] += w;
	  pvm[iz][iy] += pv;
	  h2om[iz][iy] += h2o;
	  o3m[iz][iy] += o3;
	  lwcm[iz][iy] += lwc;
	  iwcm[iz][iy] += iwc;
	  psm[iz][iy] += ps;
	  if (gsl_finite(pt)) {
	    ptm[iz][iy] += pt;
	    ztm[iz][iy] += zt;
	    ttm[iz][iy] += tt;
	    h2otm[iz][iy] += h2ot;
	    npt[iz][iy]++;
	  }
	  np[iz][iy]++;
	}
  }

  /* Create output file... */
  printf("Write meteorological data file: %s\n", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = altitude [km]\n"
	  "# $3 = longitude [deg]\n"
	  "# $4 = latitude [deg]\n"
	  "# $5 = pressure [hPa]\n"
	  "# $6 = temperature [K]\n"
	  "# $7 = zonal wind [m/s]\n"
	  "# $8 = meridional wind [m/s]\n" "# $9 = vertical wind [hPa/s]\n");
  fprintf(out,
	  "# $10 = H2O volume mixing ratio [1]\n"
	  "# $11 = O3 volume mixing ratio [1]\n"
	  "# $12 = geopotential height [km]\n"
	  "# $13 = potential vorticity [PVU]\n"
	  "# $14 = surface pressure [hPa]\n"
	  "# $15 = tropopause pressure [hPa]\n"
	  "# $16 = tropopause geopotential height [km]\n"
	  "# $17 = tropopause temperature [K]\n"
	  "# $18 = tropopause water vapor [ppv]\n"
	  "# $19 = cloud liquid water content [kg/kg]\n"
	  "# $20 = cloud ice water content [kg/kg]\n");

  /* Write data... */
  for (iz = 0; iz < nz; iz++) {
    fprintf(out, "\n");
    for (iy = 0; iy < ny; iy++)
      fprintf(out,
	      "%.2f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	      timem[iz][iy] / np[iz][iy], Z(plev[iz]), 0.0, lats[iy],
	      plev[iz], tm[iz][iy] / np[iz][iy], um[iz][iy] / np[iz][iy],
	      vm[iz][iy] / np[iz][iy], wm[iz][iy] / np[iz][iy],
	      h2om[iz][iy] / np[iz][iy], o3m[iz][iy] / np[iz][iy],
	      zm[iz][iy] / np[iz][iy], pvm[iz][iy] / np[iz][iy],
	      psm[iz][iy] / np[iz][iy], ptm[iz][iy] / npt[iz][iy],
	      ztm[iz][iy] / npt[iz][iy], ttm[iz][iy] / npt[iz][iy],
	      h2otm[iz][iy] / npt[iz][iy], lwcm[iz][iy] / np[iz][iy],
	      iwcm[iz][iy] / np[iz][iy]);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(met);

  return EXIT_SUCCESS;
}
