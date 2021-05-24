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
  Extract map from meteorological data.
*/

#include "libtrac.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Maximum number of longitudes. */
#define NX 1441

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

  static double timem[NX][NY], p0, ps, psm[NX][NY], ts, tsm[NX][NY], zs,
    zsm[NX][NY], us, usm[NX][NY], vs, vsm[NX][NY], pbl, pblm[NX][NY], pt,
    ptm[NX][NY], t, pm[NX][NY], tm[NX][NY], u, um[NX][NY], v, vm[NX][NY],
    w, wm[NX][NY], h2o, h2om[NX][NY], h2ot, h2otm[NX][NY], o3, o3m[NX][NY],
    lwc, lwcm[NX][NY], iwc, iwcm[NX][NY], z, zm[NX][NY], pv, pvm[NX][NY],
    zt, ztm[NX][NY], tt, ttm[NX][NY], pc, pcm[NX][NY], cl, clm[NX][NY],
    plcl, plclm[NX][NY], plfc, plfcm[NX][NY], pel, pelm[NX][NY],
    cape, capem[NX][NY], theta, ptop, pbot, t0, lon, lon0, lon1,
    lons[NX], dlon, lat, lat0, lat1, lats[NY], dlat, cw[3];

  static int i, ix, iy, np[NX][NY], nx, ny, ci[3];

  /* Allocate... */
  ALLOC(met, met_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <map.tab> <met0> [ <met1> ... ]");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  p0 = P(scan_ctl(argv[1], argc, argv, "MAP_Z0", -1, "10", NULL));
  lon0 = scan_ctl(argv[1], argc, argv, "MAP_LON0", -1, "-180", NULL);
  lon1 = scan_ctl(argv[1], argc, argv, "MAP_LON1", -1, "180", NULL);
  dlon = scan_ctl(argv[1], argc, argv, "MAP_DLON", -1, "-999", NULL);
  lat0 = scan_ctl(argv[1], argc, argv, "MAP_LAT0", -1, "-90", NULL);
  lat1 = scan_ctl(argv[1], argc, argv, "MAP_LAT1", -1, "90", NULL);
  dlat = scan_ctl(argv[1], argc, argv, "MAP_DLAT", -1, "-999", NULL);
  theta = scan_ctl(argv[1], argc, argv, "MAP_THETA", -1, "-999", NULL);

  /* Loop over files... */
  for (i = 3; i < argc; i++) {

    /* Read meteorological data... */
    if (!read_met(&ctl, argv[i], met))
      continue;

    /* Set horizontal grid... */
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
      if ((++nx) > NX)
	ERRMSG("Too many longitudes!");
    }
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
    for (ix = 0; ix < nx; ix++)
      for (iy = 0; iy < ny; iy++) {

	/* Find pressure level for given theta level... */
	if (theta > 0) {
	  ptop = met->p[met->np - 1];
	  pbot = met->p[0];
	  do {
	    p0 = 0.5 * (ptop + pbot);
	    intpol_met_space_3d(met, met->t, p0, lons[ix], lats[iy],
				&t0, ci, cw, 1);
	    if (THETA(p0, t0) > theta)
	      ptop = p0;
	    else
	      pbot = p0;
	  } while (fabs(ptop - pbot) > 1e-5);
	}

	/* Interpolate meteo data... */
	INTPOL_SPACE_ALL(p0, lons[ix], lats[iy]);

	/* Averaging... */
	timem[ix][iy] += met->time;
	zm[ix][iy] += z;
	pm[ix][iy] += p0;
	tm[ix][iy] += t;
	um[ix][iy] += u;
	vm[ix][iy] += v;
	wm[ix][iy] += w;
	pvm[ix][iy] += pv;
	h2om[ix][iy] += h2o;
	o3m[ix][iy] += o3;
	lwcm[ix][iy] += lwc;
	iwcm[ix][iy] += iwc;
	psm[ix][iy] += ps;
	tsm[ix][iy] += ts;
	zsm[ix][iy] += zs;
	usm[ix][iy] += us;
	vsm[ix][iy] += vs;
	pblm[ix][iy] += pbl;
	ptm[ix][iy] += pt;
	pcm[ix][iy] += pc;
	clm[ix][iy] += cl;
	plclm[ix][iy] += plcl;
	plfcm[ix][iy] += plfc;
	pelm[ix][iy] += pel;
	capem[ix][iy] += cape;
	ztm[ix][iy] += zt;
	ttm[ix][iy] += tt;
	h2otm[ix][iy] += h2ot;
	np[ix][iy]++;
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
	  "# $8 = meridional wind [m/s]\n"
	  "# $9 = vertical velocity [hPa/s]\n"
	  "# $10 = H2O volume mixing ratio [ppv]\n");
  fprintf(out,
	  "# $11 = O3 volume mixing ratio [ppv]\n"
	  "# $12 = geopotential height [km]\n"
	  "# $13 = potential vorticity [PVU]\n"
	  "# $14 = surface pressure [hPa]\n"
	  "# $15 = surface temperature [K]\n"
	  "# $16 = surface geopotential height [km]\n"
	  "# $17 = surface zonal wind [m/s]\n"
	  "# $18 = surface meridional wind [m/s]\n"
	  "# $19 = tropopause pressure [hPa]\n"
	  "# $20 = tropopause geopotential height [km]\n");
  fprintf(out,
	  "# $21 = tropopause temperature [K]\n"
	  "# $22 = tropopause water vapor [ppv]\n"
	  "# $23 = cloud liquid water content [kg/kg]\n"
	  "# $24 = cloud ice water content [kg/kg]\n"
	  "# $25 = total column cloud water [kg/m^2]\n"
	  "# $26 = cloud top pressure [hPa]\n"
	  "# $27 = pressure at lifted condensation level (LCL) [hPa]\n"
	  "# $28 = pressure at level of free convection (LFC) [hPa]\n"
	  "# $29 = pressure at equilibrium level (EL) [hPa]\n"
	  "# $30 = convective available potential energy (CAPE) [J/kg]\n");
  fprintf(out,
	  "# $31 = relative humidity over water [%%]\n"
	  "# $32 = relative humidity over ice [%%]\n"
	  "# $33 = dew point temperature [K]\n"
	  "# $34 = frost point temperature [K]\n"
	  "# $35 = boundary layer pressure [hPa]\n");

  /* Write data... */
  for (iy = 0; iy < ny; iy++) {
    fprintf(out, "\n");
    for (ix = 0; ix < nx; ix++)
      fprintf(out,
	      "%.2f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g"
	      " %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	      timem[ix][iy] / np[ix][iy], Z(pm[ix][iy] / np[ix][iy]),
	      lons[ix], lats[iy], pm[ix][iy] / np[ix][iy],
	      tm[ix][iy] / np[ix][iy], um[ix][iy] / np[ix][iy],
	      vm[ix][iy] / np[ix][iy], wm[ix][iy] / np[ix][iy],
	      h2om[ix][iy] / np[ix][iy], o3m[ix][iy] / np[ix][iy],
	      zm[ix][iy] / np[ix][iy], pvm[ix][iy] / np[ix][iy],
	      psm[ix][iy] / np[ix][iy], tsm[ix][iy] / np[ix][iy],
	      zsm[ix][iy] / np[ix][iy], usm[ix][iy] / np[ix][iy],
	      vsm[ix][iy] / np[ix][iy], ptm[ix][iy] / np[ix][iy],
	      ztm[ix][iy] / np[ix][iy], ttm[ix][iy] / np[ix][iy],
	      h2otm[ix][iy] / np[ix][iy], lwcm[ix][iy] / np[ix][iy],
	      iwcm[ix][iy] / np[ix][iy], clm[ix][iy] / np[ix][iy],
	      pcm[ix][iy] / np[ix][iy], plclm[ix][iy] / np[ix][iy],
	      plfcm[ix][iy] / np[ix][iy], pelm[ix][iy] / np[ix][iy],
	      capem[ix][iy] / np[ix][iy],
	      RH(pm[ix][iy] / np[ix][iy], tm[ix][iy] / np[ix][iy],
		 h2om[ix][iy] / np[ix][iy]),
	      RHICE(pm[ix][iy] / np[ix][iy], tm[ix][iy] / np[ix][iy],
		    h2om[ix][iy] / np[ix][iy]),
	      TDEW(pm[ix][iy] / np[ix][iy], h2om[ix][iy] / np[ix][iy]),
	      TICE(pm[ix][iy] / np[ix][iy], h2om[ix][iy] / np[ix][iy]),
	      pblm[ix][iy] / np[ix][iy]);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(met);

  return EXIT_SUCCESS;
}
