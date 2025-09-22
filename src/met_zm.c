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
  Extract zonal mean from meteorological data.
*/

#include "mptrac.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Maximum number of altitudes. */
#define NZ 1000

/*! Maximum number of latitudes. */
#define NY EY

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  clim_t *clim;

  met_t *met;
  
  dd_t *dd;

  FILE *out;

  static double timem[NZ][NY], psm[NZ][NY], tsm[NZ][NY], zsm[NZ][NY],
    usm[NZ][NY], vsm[NZ][NY], essm[NZ][NY], nssm[NZ][NY], shfm[NZ][NY],
    lsmm[NZ][NY], sstm[NZ][NY], pblm[NZ][NY],
    ptm[NZ][NY], pctm[NZ][NY], pcbm[NZ][NY], clm[NZ][NY], plclm[NZ][NY],
    plfcm[NZ][NY], pelm[NZ][NY], capem[NZ][NY], cinm[NZ][NY], o3cm[NZ][NY],
    ttm[NZ][NY], ztm[NZ][NY], tm[NZ][NY], um[NZ][NY], vm[NZ][NY], wm[NZ][NY],
    h2om[NZ][NY], h2otm[NZ][NY], pvm[NZ][NY], o3m[NZ][NY], lwcm[NZ][NY],
    rwcm[NZ][NY], iwcm[NZ][NY], swcm[NZ][NY], ccm[NZ][NY], zm[NZ][NY],
    rhm[NZ][NY], rhicem[NZ][NY], tdewm[NZ][NY], ticem[NZ][NY], tnatm[NZ][NY],
    hno3m[NZ][NY], ohm[NZ][NY], h2o2m[NZ][NY], ho2m[NZ][NY], o1dm[NZ][NY], z,
    zt, tt, plev[NZ], ps, ts, zs, us, vs, ess, nss, shf, lsm, sst, pbl,
    pt, pct, pcb, plcl, plfc, pel, cape, cin, o3c, cl, t, u, v, w, pv,
    h2o, h2ot, o3, lwc, rwc, iwc, swc, cc, lat, lats[NY], lonm[NZ][NY], cw[3];

  static int np[NZ][NY], npc[NZ][NY], npt[NZ][NY], ny, nz, ci[3];

  /* Allocate... */
  ALLOC(clim, clim_t, 1);
  ALLOC(met, met_t, 1);
  ALLOC(dd, dd_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <zm.tab> <met0> [ <met1> ... ]");

  /* Read control parameters... */
  mptrac_read_ctl(argv[1], argc, argv, &ctl);
  double z0 = scan_ctl(argv[1], argc, argv, "ZM_Z0", -1, "-999", NULL);
  double z1 = scan_ctl(argv[1], argc, argv, "ZM_Z1", -1, "-999", NULL);
  double dz = scan_ctl(argv[1], argc, argv, "ZM_DZ", -1, "-999", NULL);
  double lon0 = scan_ctl(argv[1], argc, argv, "ZM_LON0", -1, "-360", NULL);
  double lon1 = scan_ctl(argv[1], argc, argv, "ZM_LON1", -1, "360", NULL);
  double lat0 = scan_ctl(argv[1], argc, argv, "ZM_LAT0", -1, "-90", NULL);
  double lat1 = scan_ctl(argv[1], argc, argv, "ZM_LAT1", -1, "90", NULL);
  double dlat = scan_ctl(argv[1], argc, argv, "ZM_DLAT", -1, "-999", NULL);

  /* Read climatological data... */
  mptrac_read_clim(&ctl, clim);

  /* Loop over files... */
  for (int i = 3; i < argc; i++) {

    /* Read meteorological data... */
    if (!mptrac_read_met(argv[i], &ctl, clim, met, dd))
      continue;

    /* Set vertical grid... */
    if (z0 < 0)
      z0 = Z(met->p[0]);
    if (z1 < 0)
      z1 = Z(met->p[met->np - 1]);
    nz = 0;
    if (dz < 0) {
      for (int iz = 0; iz < met->np; iz++)
	if (Z(met->p[iz]) >= z0 && Z(met->p[iz]) <= z1) {
	  plev[nz] = met->p[iz];
	  if ((++nz) >= NZ)
	    ERRMSG("Too many pressure levels!");
	}
    } else
      for (z = z0; z <= z1; z += dz) {
	plev[nz] = P(z);
	if ((++nz) >= NZ)
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
    for (lat = lat0; lat <= lat1 + 0.001; lat += dlat) {
      lats[ny] = round(lat * 1e3) / 1e3;
      if ((++ny) >= NY)
	ERRMSG("Too many latitudes!");
    }

    /* Average... */
    for (int ix = 0; ix < met->nx; ix++)
      if (met->lon[ix] >= lon0 && met->lon[ix] <= lon1)
	for (int iy = 0; iy < ny; iy++)
	  for (int iz = 0; iz < nz; iz++) {

	    /* Interpolate meteo data... */
	    INTPOL_SPACE_ALL(plev[iz], met->lon[ix], lats[iy]);

	    /* Averaging... */
	    timem[iz][iy] += met->time;
	    lonm[iz][iy] += met->lon[ix];
	    zm[iz][iy] += z;
	    tm[iz][iy] += t;
	    um[iz][iy] += u;
	    vm[iz][iy] += v;
	    wm[iz][iy] += w;
	    pvm[iz][iy] += pv;
	    h2om[iz][iy] += h2o;
	    o3m[iz][iy] += o3;
	    lwcm[iz][iy] += lwc;
	    rwcm[iz][iy] += rwc;
	    iwcm[iz][iy] += iwc;
	    swcm[iz][iy] += swc;
	    ccm[iz][iy] += cc;
	    psm[iz][iy] += ps;
	    tsm[iz][iy] += ts;
	    zsm[iz][iy] += zs;
	    usm[iz][iy] += us;
	    vsm[iz][iy] += vs;
	    essm[iz][iy] += ess;
	    nssm[iz][iy] += nss;
	    shfm[iz][iy] += shf;
	    lsmm[iz][iy] += lsm;
	    sstm[iz][iy] += sst;
	    pblm[iz][iy] += pbl;
	    pctm[iz][iy] += pct;
	    pcbm[iz][iy] += pcb;
	    clm[iz][iy] += cl;
	    if (isfinite(plfc) && isfinite(pel) && cape >= ctl.conv_cape
		&& (ctl.conv_cin <= 0 || cin < ctl.conv_cin)) {
	      plclm[iz][iy] += plcl;
	      plfcm[iz][iy] += plfc;
	      pelm[iz][iy] += pel;
	      capem[iz][iy] += cape;
	      cinm[iz][iy] += cin;
	      npc[iz][iy]++;
	    }
	    if (isfinite(pt)) {
	      ptm[iz][iy] += pt;
	      ztm[iz][iy] += zt;
	      ttm[iz][iy] += tt;
	      h2otm[iz][iy] += h2ot;
	      npt[iz][iy]++;
	    }
	    o3cm[iz][iy] += o3c;
	    rhm[iz][iy] += RH(plev[iz], t, h2o);
	    rhicem[iz][iy] += RHICE(plev[iz], t, h2o);
	    tdewm[iz][iy] += TDEW(plev[iz], h2o);
	    ticem[iz][iy] += TICE(plev[iz], h2o);
	    hno3m[iz][iy] +=
	      clim_zm(&clim->hno3, met->time, lats[iy], plev[iz]);
	    tnatm[iz][iy] +=
	      nat_temperature(plev[iz], h2o,
			      clim_zm(&clim->hno3, met->time, lats[iy],
				      plev[iz]));
	    ohm[iz][iy] +=
	      clim_oh(&ctl, clim, met->time, met->lon[ix], lats[iy],
		      plev[iz]);
	    h2o2m[iz][iy]
	      += clim_zm(&clim->h2o2, met->time, lats[iy], plev[iz]);
	    ho2m[iz][iy]
	      += clim_zm(&clim->ho2, met->time, lats[iy], plev[iz]);
	    o1dm[iz][iy]
	      += clim_zm(&clim->o1d, met->time, lats[iy], plev[iz]);
	    np[iz][iy]++;
	  }
  }

  /* Create output file... */
  LOG(1, "Write meteorological data file: %s", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  MET_HEADER;

  /* Write data... */
  for (int iz = 0; iz < nz; iz++) {
    fprintf(out, "\n");
    for (int iy = 0; iy < ny; iy++)
      fprintf(out,
	      "%.2f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g"
	      " %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g"
	      " %g %g %g %g %g %g %g %g %g %g %g %g %g %d %d %d\n",
	      timem[iz][iy] / np[iz][iy], Z(plev[iz]),
	      lonm[iz][iy] / np[iz][iy], lats[iy],
	      plev[iz], tm[iz][iy] / np[iz][iy], um[iz][iy] / np[iz][iy],
	      vm[iz][iy] / np[iz][iy], wm[iz][iy] / np[iz][iy],
	      h2om[iz][iy] / np[iz][iy], o3m[iz][iy] / np[iz][iy],
	      zm[iz][iy] / np[iz][iy], pvm[iz][iy] / np[iz][iy],
	      psm[iz][iy] / np[iz][iy], tsm[iz][iy] / np[iz][iy],
	      zsm[iz][iy] / np[iz][iy], usm[iz][iy] / np[iz][iy],
	      vsm[iz][iy] / np[iz][iy], essm[iz][iy] / np[iz][iy],
	      nssm[iz][iy] / np[iz][iy], shfm[iz][iy] / np[iz][iy],
	      lsmm[iz][iy] / np[iz][iy],
	      sstm[iz][iy] / np[iz][iy], ptm[iz][iy] / npt[iz][iy],
	      ztm[iz][iy] / npt[iz][iy], ttm[iz][iy] / npt[iz][iy],
	      h2otm[iz][iy] / npt[iz][iy], lwcm[iz][iy] / np[iz][iy],
	      rwcm[iz][iy] / np[iz][iy], iwcm[iz][iy] / np[iz][iy],
	      swcm[iz][iy] / np[iz][iy], ccm[iz][iy] / np[iz][iy],
	      clm[iz][iy] / np[iz][iy], pctm[iz][iy] / np[iz][iy],
	      pcbm[iz][iy] / np[iz][iy], plclm[iz][iy] / npc[iz][iy],
	      plfcm[iz][iy] / npc[iz][iy], pelm[iz][iy] / npc[iz][iy],
	      capem[iz][iy] / npc[iz][iy], cinm[iz][iy] / npc[iz][iy],
	      rhm[iz][iy] / np[iz][iy], rhicem[iz][iy] / np[iz][iy],
	      tdewm[iz][iy] / np[iz][iy], ticem[iz][iy] / np[iz][iy],
	      tnatm[iz][iy] / np[iz][iy], hno3m[iz][iy] / np[iz][iy],
	      ohm[iz][iy] / np[iz][iy], h2o2m[iz][iy] / np[iz][iy],
	      ho2m[iz][iy] / np[iz][iy], o1dm[iz][iy] / np[iz][iy],
	      pblm[iz][iy] / np[iz][iy], o3cm[iz][iy] / np[iz][iy],
	      np[iz][iy], npt[iz][iy], npc[iz][iy]);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(clim);
  free(met);

  return EXIT_SUCCESS;
}
