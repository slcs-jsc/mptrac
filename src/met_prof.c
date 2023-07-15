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
  Extract vertical profile from meteorological data.
*/

#include "libtrac.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Maximum number of altitudes. */
#define NZ 1000

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  clim_t *clim;

  met_t *met;

  FILE *out;

  static double timem[NZ], z, z0, z1, dz, lon, lon0, lon1, dlon, lonm[NZ],
    lat, lat0, lat1, dlat, latm[NZ], t, tm[NZ], u, um[NZ], v, vm[NZ], w,
    wm[NZ], h2o, h2om[NZ], h2ot, h2otm[NZ], o3, o3m[NZ], lwc, lwcm[NZ],
    iwc, iwcm[NZ], cc, ccm[NZ], ps, psm[NZ], ts, tsm[NZ], zs, zsm[NZ],
    us, usm[NZ], vs, vsm[NZ], lsm, lsmm[NZ], sst, sstm[NZ], pbl, pblm[NZ],
    pt, ptm[NZ], pct, pctm[NZ], pcb, pcbm[NZ], cl, clm[NZ], plcl, plclm[NZ],
    plfc, plfcm[NZ], pel, pelm[NZ], cape, capem[NZ], cin, cinm[NZ], tt,
    ttm[NZ], zm[NZ], zt, ztm[NZ], pv, pvm[NZ], plev[NZ], rhm[NZ], rhicem[NZ],
    tdewm[NZ], ticem[NZ], tnatm[NZ], hno3m[NZ], ohm[NZ], h2o2m[NZ], cw[3];

  static int i, iz, np[NZ], npc[NZ], npt[NZ], nz, ci[3];

  /* Allocate... */
  ALLOC(clim, clim_t, 1);
  ALLOC(met, met_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <prof.tab> <met0> [ <met1> ... ]");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  z0 = scan_ctl(argv[1], argc, argv, "PROF_Z0", -1, "-999", NULL);
  z1 = scan_ctl(argv[1], argc, argv, "PROF_Z1", -1, "-999", NULL);
  dz = scan_ctl(argv[1], argc, argv, "PROF_DZ", -1, "-999", NULL);
  lon0 = scan_ctl(argv[1], argc, argv, "PROF_LON0", -1, "0", NULL);
  lon1 = scan_ctl(argv[1], argc, argv, "PROF_LON1", -1, "0", NULL);
  dlon = scan_ctl(argv[1], argc, argv, "PROF_DLON", -1, "-999", NULL);
  lat0 = scan_ctl(argv[1], argc, argv, "PROF_LAT0", -1, "0", NULL);
  lat1 = scan_ctl(argv[1], argc, argv, "PROF_LAT1", -1, "0", NULL);
  dlat = scan_ctl(argv[1], argc, argv, "PROF_DLAT", -1, "-999", NULL);

  /* Read climatological data... */
  read_clim(&ctl, clim);

  /* Loop over input files... */
  for (i = 3; i < argc; i++) {

    /* Read meteorological data... */
    if (!read_met(argv[i], &ctl, clim, met))
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
    if (dlon <= 0)
      dlon = fabs(met->lon[1] - met->lon[0]);
    if (dlat <= 0)
      dlat = fabs(met->lat[1] - met->lat[0]);

    /* Average... */
    for (iz = 0; iz < nz; iz++)
      for (lon = lon0; lon <= lon1; lon += dlon)
	for (lat = lat0; lat <= lat1; lat += dlat) {

	  /* Interpolate meteo data... */
	  INTPOL_SPACE_ALL(plev[iz], lon, lat);

	  /* Averaging... */
	  if (gsl_finite(t) && gsl_finite(u)
	      && gsl_finite(v) && gsl_finite(w)) {
	    timem[iz] += met->time;
	    lonm[iz] += lon;
	    latm[iz] += lat;
	    zm[iz] += z;
	    tm[iz] += t;
	    um[iz] += u;
	    vm[iz] += v;
	    wm[iz] += w;
	    pvm[iz] += pv;
	    h2om[iz] += h2o;
	    o3m[iz] += o3;
	    lwcm[iz] += lwc;
	    iwcm[iz] += iwc;
	    ccm[iz] += cc;
	    psm[iz] += ps;
	    tsm[iz] += ts;
	    zsm[iz] += zs;
	    usm[iz] += us;
	    vsm[iz] += vs;
	    lsmm[iz] += lsm;
	    sstm[iz] += sst;
	    pblm[iz] += pbl;
	    pctm[iz] += pct;
	    pcbm[iz] += pcb;
	    clm[iz] += cl;
	    if (gsl_finite(plfc) && gsl_finite(pel) && cape >= ctl.conv_cape
		&& (ctl.conv_cin <= 0 || cin < ctl.conv_cin)) {
	      plclm[iz] += plcl;
	      plfcm[iz] += plfc;
	      pelm[iz] += pel;
	      capem[iz] += cape;
	      cinm[iz] += cin;
	      npc[iz]++;
	    }
	    if (gsl_finite(pt)) {
	      ptm[iz] += pt;
	      ztm[iz] += zt;
	      ttm[iz] += tt;
	      h2otm[iz] += h2ot;
	      npt[iz]++;
	    }
	    rhm[iz] += RH(plev[iz], t, h2o);
	    rhicem[iz] += RHICE(plev[iz], t, h2o);
	    tdewm[iz] += TDEW(plev[iz], h2o);
	    ticem[iz] += TICE(plev[iz], h2o);
	    hno3m[iz] += clim_hno3(clim, met->time, lat, plev[iz]);
	    tnatm[iz] +=
	      nat_temperature(plev[iz], h2o,
			      clim_hno3(clim, met->time, lat, plev[iz]));
	    ohm[iz] += clim_oh(&ctl, clim, met->time, lon, lat, plev[iz]);
	    h2o2m[iz] += clim_h2o2(clim, met->time, lat, plev[iz]);
	    np[iz]++;
	  }
	}
  }

  /* Create output file... */
  LOG(1, "Write meteorological data file: %s", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  MET_HEADER;

  /* Write data... */
  fprintf(out, "\n");
  for (iz = 0; iz < nz; iz++)
    fprintf(out,
	    "%.2f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g"
	    " %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g"
	    " %g %g %g %g %g %g %g %g %g %g %d %d %d\n",
	    timem[iz] / np[iz], Z(plev[iz]), lonm[iz] / np[iz],
	    latm[iz] / np[iz], plev[iz], tm[iz] / np[iz], um[iz] / np[iz],
	    vm[iz] / np[iz], wm[iz] / np[iz], h2om[iz] / np[iz],
	    o3m[iz] / np[iz], zm[iz] / np[iz], pvm[iz] / np[iz],
	    psm[iz] / np[iz], tsm[iz] / np[iz], zsm[iz] / np[iz],
	    usm[iz] / np[iz], vsm[iz] / np[iz], lsmm[iz] / np[iz],
	    sstm[iz] / np[iz], ptm[iz] / npt[iz], ztm[iz] / npt[iz],
	    ttm[iz] / npt[iz], h2otm[iz] / npt[iz],
	    lwcm[iz] / np[iz], iwcm[iz] / np[iz], ccm[iz] / np[iz],
	    clm[iz] / np[iz], pctm[iz] / np[iz], pcbm[iz] / np[iz],
	    plclm[iz] / npc[iz], plfcm[iz] / npc[iz], pelm[iz] / npc[iz],
	    capem[iz] / npc[iz], cinm[iz] / npc[iz], rhm[iz] / np[iz],
	    rhicem[iz] / np[iz], tdewm[iz] / np[iz], ticem[iz] / np[iz],
	    tnatm[iz] / np[iz], hno3m[iz] / np[iz], ohm[iz] / np[iz],
	    h2o2m[iz] / np[iz], pblm[iz] / np[iz], np[iz], npt[iz], npc[iz]);

  /* Close file... */
  fclose(out);

  /* Free... */
  free(clim);
  free(met);

  return EXIT_SUCCESS;
}
