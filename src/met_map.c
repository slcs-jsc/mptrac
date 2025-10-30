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
  Extract map from meteorological data.
*/

#include "mptrac.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Maximum number of longitudes. */
#define NX EX

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

  static double *timem, ps, *psm, ts, *tsm, zs, *zsm, us, *usm, vs, *vsm,
    ess, *essm, nss, *nssm, shf, *shfm, lsm, *lsmm, sst, *sstm, pbl, *pblm,
    pt, *ptm, t, *pm, *tm, u, *um, v, *vm, w, *wm, h2o, *h2om, h2ot, *h2otm,
    o3, *o3m, *hno3m, *ohm, *h2o2m, *ho2m, *o1dm, *tdewm, *ticem, *tnatm,
    lwc, *lwcm, rwc, *rwcm, iwc, *iwcm, swc, *swcm, cc, *ccm, z, *zm,
    pv, *pvm, zt, *ztm, tt, *ttm, pct, *pctm, pcb, *pcbm, cl, *clm,
    plcl, *plclm, plfc, *plfcm, pel, *pelm, cape, *capem, cin, *cinm,
    o3c, *o3cm, *rhm, *rhicem, ptop, pbot, t0, lon, lons[NX], lat, lats[NY];

  static int *np, *npc, *npt, nx, ny;

  /* Allocate... */
  ALLOC(clim, clim_t, 1);
  ALLOC(met, met_t, 1);
  ALLOC(dd, dd_t, 1);
  ALLOC(timem, double,
	NX * NY);
  ALLOC(psm, double,
	NX * NY);
  ALLOC(tsm, double,
	NX * NY);
  ALLOC(zsm, double,
	NX * NY);
  ALLOC(usm, double,
	NX * NY);
  ALLOC(vsm, double,
	NX * NY);
  ALLOC(essm, double,
	NX * NY);
  ALLOC(nssm, double,
	NX * NY);
  ALLOC(shfm, double,
	NX * NY);
  ALLOC(lsmm, double,
	NX * NY);
  ALLOC(sstm, double,
	NX * NY);
  ALLOC(pblm, double,
	NX * NY);
  ALLOC(ptm, double,
	NX * NY);
  ALLOC(pm, double,
	NX * NY);
  ALLOC(tm, double,
	NX * NY);
  ALLOC(um, double,
	NX * NY);
  ALLOC(vm, double,
	NX * NY);
  ALLOC(wm, double,
	NX * NY);
  ALLOC(h2om, double,
	NX * NY);
  ALLOC(h2otm, double,
	NX * NY);
  ALLOC(o3m, double,
	NX * NY);
  ALLOC(hno3m, double,
	NX * NY);
  ALLOC(ohm, double,
	NX * NY);
  ALLOC(h2o2m, double,
	NX * NY);
  ALLOC(ho2m, double,
	NX * NY);
  ALLOC(o1dm, double,
	NX * NY);
  ALLOC(tdewm, double,
	NX * NY);
  ALLOC(ticem, double,
	NX * NY);
  ALLOC(tnatm, double,
	NX * NY);
  ALLOC(lwcm, double,
	NX * NY);
  ALLOC(rwcm, double,
	NX * NY);
  ALLOC(iwcm, double,
	NX * NY);
  ALLOC(swcm, double,
	NX * NY);
  ALLOC(ccm, double,
	NX * NY);
  ALLOC(zm, double,
	NX * NY);
  ALLOC(pvm, double,
	NX * NY);
  ALLOC(ztm, double,
	NX * NY);
  ALLOC(ttm, double,
	NX * NY);
  ALLOC(pctm, double,
	NX * NY);
  ALLOC(pcbm, double,
	NX * NY);
  ALLOC(clm, double,
	NX * NY);
  ALLOC(plclm, double,
	NX * NY);
  ALLOC(plfcm, double,
	NX * NY);
  ALLOC(pelm, double,
	NX * NY);
  ALLOC(capem, double,
	NX * NY);
  ALLOC(cinm, double,
	NX * NY);
  ALLOC(o3cm, double,
	NX * NY);
  ALLOC(rhm, double,
	NX * NY);
  ALLOC(rhicem, double,
	NX * NY);
  ALLOC(np, int,
	NX * NY);
  ALLOC(npc, int,
	NX * NY);
  ALLOC(npt, int,
	NX * NY);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <map.tab> <met0> [ <met1> ... ]");

  /* Read control parameters... */
  mptrac_read_ctl(argv[1], argc, argv, &ctl);
  double p0 = P(scan_ctl(argv[1], argc, argv, "MAP_Z0", -1, "10", NULL));
  double lon0 = scan_ctl(argv[1], argc, argv, "MAP_LON0", -1, "-180", NULL);
  double lon1 = scan_ctl(argv[1], argc, argv, "MAP_LON1", -1, "180", NULL);
  double dlon = scan_ctl(argv[1], argc, argv, "MAP_DLON", -1, "-999", NULL);
  double lat0 = scan_ctl(argv[1], argc, argv, "MAP_LAT0", -1, "-90", NULL);
  double lat1 = scan_ctl(argv[1], argc, argv, "MAP_LAT1", -1, "90", NULL);
  double dlat = scan_ctl(argv[1], argc, argv, "MAP_DLAT", -1, "-999", NULL);
  double theta = scan_ctl(argv[1], argc, argv, "MAP_THETA", -1, "-999", NULL);

  /* Read climatological data... */
  mptrac_read_clim(&ctl, clim);

  /* Loop over files... */
  for (int i = 3; i < argc; i++) {

    /* Read meteorological data... */
    if (!mptrac_read_met(argv[i], &ctl, clim, met, dd))
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
    for (lon = lon0; lon <= lon1 + 0.001; lon += dlon) {
      lons[nx] = round(lon * 1e3) / 1e3;
      if ((++nx) >= NX)
	ERRMSG("Too many longitudes!");
    }
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
    for (int ix = 0; ix < nx; ix++)
      for (int iy = 0; iy < ny; iy++) {

	/* Find pressure level for given theta level... */
	INTPOL_INIT;
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
	timem[iy * nx + ix] += met->time;
	zm[iy * nx + ix] += z;
	pm[iy * nx + ix] += p0;
	tm[iy * nx + ix] += t;
	um[iy * nx + ix] += u;
	vm[iy * nx + ix] += v;
	wm[iy * nx + ix] += w;
	pvm[iy * nx + ix] += pv;
	h2om[iy * nx + ix] += h2o;
	o3m[iy * nx + ix] += o3;
	lwcm[iy * nx + ix] += lwc;
	rwcm[iy * nx + ix] += rwc;
	iwcm[iy * nx + ix] += iwc;
	swcm[iy * nx + ix] += swc;
	ccm[iy * nx + ix] += cc;
	psm[iy * nx + ix] += ps;
	tsm[iy * nx + ix] += ts;
	zsm[iy * nx + ix] += zs;
	usm[iy * nx + ix] += us;
	vsm[iy * nx + ix] += vs;
	essm[iy * nx + ix] += ess;
	nssm[iy * nx + ix] += nss;
	shfm[iy * nx + ix] += shf;
	lsmm[iy * nx + ix] += lsm;
	sstm[iy * nx + ix] += sst;
	pblm[iy * nx + ix] += pbl;
	pctm[iy * nx + ix] += pct;
	pcbm[iy * nx + ix] += pcb;
	clm[iy * nx + ix] += cl;
	if (isfinite(plfc) && isfinite(pel) && cape >= ctl.conv_cape
	    && (ctl.conv_cin <= 0 || cin < ctl.conv_cin)) {
	  plclm[iy * nx + ix] += plcl;
	  plfcm[iy * nx + ix] += plfc;
	  pelm[iy * nx + ix] += pel;
	  capem[iy * nx + ix] += cape;
	  cinm[iy * nx + ix] += cin;
	  npc[iy * nx + ix]++;
	}
	if (isfinite(pt)) {
	  ptm[iy * nx + ix] += pt;
	  ztm[iy * nx + ix] += zt;
	  ttm[iy * nx + ix] += tt;
	  h2otm[iy * nx + ix] += h2ot;
	  npt[iy * nx + ix]++;
	}
	o3cm[iy * nx + ix] += o3c;
	hno3m[iy * nx + ix] += clim_zm(&clim->hno3, met->time, lats[iy], p0);
	tnatm[iy * nx + ix] +=
	  nat_temperature(p0, h2o,
			  clim_zm(&clim->hno3, met->time, lats[iy], p0));
	ohm[iy * nx + ix] +=
	  clim_oh(&ctl, clim, met->time, lons[ix], lats[iy], p0);
	h2o2m[iy * nx + ix] += clim_zm(&clim->h2o2, met->time, lats[iy], p0);
	ho2m[iy * nx + ix] += clim_zm(&clim->ho2, met->time, lats[iy], p0);
	o1dm[iy * nx + ix] += clim_zm(&clim->o1d, met->time, lats[iy], p0);
	rhm[iy * nx + ix] += RH(p0, t, h2o);
	rhicem[iy * nx + ix] += RHICE(p0, t, h2o);
	tdewm[iy * nx + ix] += TDEW(p0, h2o);
	ticem[iy * nx + ix] += TICE(p0, h2o);
	np[iy * nx + ix]++;
      }
  }

  /* Create output file... */
  LOG(1, "Write meteorological data file: %s", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  MET_HEADER;

  /* Write data... */
  for (int iy = 0; iy < ny; iy++) {
    fprintf(out, "\n");
    for (int ix = 0; ix < nx; ix++)
      fprintf(out,
	      "%.2f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g"
	      " %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g"
	      " %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d %d %d\n",
	      timem[iy * nx + ix] / np[iy * nx + ix],
	      Z(pm[iy * nx + ix] / np[iy * nx + ix]), lons[ix], lats[iy],
	      pm[iy * nx + ix] / np[iy * nx + ix],
	      tm[iy * nx + ix] / np[iy * nx + ix],
	      um[iy * nx + ix] / np[iy * nx + ix],
	      vm[iy * nx + ix] / np[iy * nx + ix],
	      wm[iy * nx + ix] / np[iy * nx + ix],
	      h2om[iy * nx + ix] / np[iy * nx + ix],
	      o3m[iy * nx + ix] / np[iy * nx + ix],
	      zm[iy * nx + ix] / np[iy * nx + ix],
	      pvm[iy * nx + ix] / np[iy * nx + ix],
	      psm[iy * nx + ix] / np[iy * nx + ix],
	      tsm[iy * nx + ix] / np[iy * nx + ix],
	      zsm[iy * nx + ix] / np[iy * nx + ix],
	      usm[iy * nx + ix] / np[iy * nx + ix],
	      vsm[iy * nx + ix] / np[iy * nx + ix],
	      essm[iy * nx + ix] / np[iy * nx + ix],
	      nssm[iy * nx + ix] / np[iy * nx + ix],
	      shfm[iy * nx + ix] / np[iy * nx + ix],
	      lsmm[iy * nx + ix] / np[iy * nx + ix],
	      sstm[iy * nx + ix] / np[iy * nx + ix],
	      ptm[iy * nx + ix] / npt[iy * nx + ix],
	      ztm[iy * nx + ix] / npt[iy * nx + ix],
	      ttm[iy * nx + ix] / npt[iy * nx + ix],
	      h2otm[iy * nx + ix] / npt[iy * nx + ix],
	      lwcm[iy * nx + ix] / np[iy * nx + ix],
	      rwcm[iy * nx + ix] / np[iy * nx + ix],
	      iwcm[iy * nx + ix] / np[iy * nx + ix],
	      swcm[iy * nx + ix] / np[iy * nx + ix],
	      ccm[iy * nx + ix] / np[iy * nx + ix],
	      clm[iy * nx + ix] / np[iy * nx + ix],
	      pctm[iy * nx + ix] / np[iy * nx + ix],
	      pcbm[iy * nx + ix] / np[iy * nx + ix],
	      plclm[iy * nx + ix] / npc[iy * nx + ix],
	      plfcm[iy * nx + ix] / npc[iy * nx + ix],
	      pelm[iy * nx + ix] / npc[iy * nx + ix],
	      capem[iy * nx + ix] / npc[iy * nx + ix],
	      cinm[iy * nx + ix] / npc[iy * nx + ix],
	      rhm[iy * nx + ix] / np[iy * nx + ix],
	      rhicem[iy * nx + ix] / np[iy * nx + ix],
	      tdewm[iy * nx + ix] / np[iy * nx + ix],
	      ticem[iy * nx + ix] / np[iy * nx + ix],
	      tnatm[iy * nx + ix] / np[iy * nx + ix],
	      hno3m[iy * nx + ix] / np[iy * nx + ix],
	      ohm[iy * nx + ix] / np[iy * nx + ix],
	      h2o2m[iy * nx + ix] / np[iy * nx + ix],
	      ho2m[iy * nx + ix] / np[iy * nx + ix],
	      o1dm[iy * nx + ix] / np[iy * nx + ix],
	      pblm[iy * nx + ix] / np[iy * nx + ix],
	      o3cm[iy * nx + ix] / np[iy * nx + ix], np[iy * nx + ix],
	      npt[iy * nx + ix], npc[iy * nx + ix]);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(clim);
  free(met);
  free(dd);
  free(timem);
  free(psm);
  free(tsm);
  free(zsm);
  free(usm);
  free(vsm);
  free(essm);
  free(nssm);
  free(shfm);
  free(lsmm);
  free(sstm);
  free(pblm);
  free(ptm);
  free(pm);
  free(tm);
  free(um);
  free(vm);
  free(wm);
  free(h2om);
  free(h2otm);
  free(o3m);
  free(hno3m);
  free(ohm);
  free(h2o2m);
  free(ho2m);
  free(o1dm);
  free(tdewm);
  free(ticem);
  free(tnatm);
  free(lwcm);
  free(rwcm);
  free(iwcm);
  free(swcm);
  free(ccm);
  free(zm);
  free(pvm);
  free(ztm);
  free(ttm);
  free(pctm);
  free(pcbm);
  free(clm);
  free(plclm);
  free(plfcm);
  free(pelm);
  free(capem);
  free(cinm);
  free(o3cm);
  free(rhm);
  free(rhicem);
  free(np);
  free(npc);
  free(npt);

  return EXIT_SUCCESS;
}
