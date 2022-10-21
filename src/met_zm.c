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

  clim_t *clim;

  FILE *out;

  static double timem[NZ][NY], psm[NZ][NY], tsm[NZ][NY], zsm[NZ][NY],
    usm[NZ][NY], vsm[NZ][NY], pblm[NZ][NY], ptm[NZ][NY], pctm[NZ][NY],
    pcbm[NZ][NY], clm[NZ][NY], plclm[NZ][NY], plfcm[NZ][NY], pelm[NZ][NY],
    capem[NZ][NY], cinm[NZ][NY], ttm[NZ][NY], ztm[NZ][NY], tm[NZ][NY],
    um[NZ][NY], vm[NZ][NY], wm[NZ][NY], h2om[NZ][NY], h2otm[NZ][NY],
    pvm[NZ][NY], o3m[NZ][NY], lwcm[NZ][NY], iwcm[NZ][NY], zm[NZ][NY],
    rhm[NZ][NY], rhicem[NZ][NY], tdewm[NZ][NY], ticem[NZ][NY], tnatm[NZ][NY],
    hno3m[NZ][NY], ohm[NZ][NY], z, z0, z1, dz, zt, tt, plev[NZ],
    ps, ts, zs, us, vs, pbl, pt, pct, pcb, plcl, plfc, pel,
    cape, cin, cl, t, u, v, w, pv, h2o, h2ot, o3, lwc, iwc,
    lat, lat0, lat1, dlat, lats[NY], lon0, lon1, lonm[NZ][NY], cw[3];

  static int i, ix, iy, iz, np[NZ][NY], npc[NZ][NY], npt[NZ][NY], ny, nz,
    ci[3];

  /* Allocate... */
  ALLOC(met, met_t, 1);
  ALLOC(clim, clim_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <zm.tab> <met0> [ <met1> ... ]");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  z0 = scan_ctl(argv[1], argc, argv, "ZM_Z0", -1, "-999", NULL);
  z1 = scan_ctl(argv[1], argc, argv, "ZM_Z1", -1, "-999", NULL);
  dz = scan_ctl(argv[1], argc, argv, "ZM_DZ", -1, "-999", NULL);
  lon0 = scan_ctl(argv[1], argc, argv, "ZM_LON0", -1, "-360", NULL);
  lon1 = scan_ctl(argv[1], argc, argv, "ZM_LON1", -1, "360", NULL);
  lat0 = scan_ctl(argv[1], argc, argv, "ZM_LAT0", -1, "-90", NULL);
  lat1 = scan_ctl(argv[1], argc, argv, "ZM_LAT1", -1, "90", NULL);
  dlat = scan_ctl(argv[1], argc, argv, "ZM_DLAT", -1, "-999", NULL);

  /* Read climatological data... */
  read_clim(&ctl, clim);

  /* Loop over files... */
  for (i = 3; i < argc; i++) {

    /* Read meteorological data... */
    if (!read_met(argv[i], &ctl, met, clim))
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
      if (met->lon[ix] >= lon0 && met->lon[ix] <= lon1)
	for (iy = 0; iy < ny; iy++)
	  for (iz = 0; iz < nz; iz++) {

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
	    iwcm[iz][iy] += iwc;
	    psm[iz][iy] += ps;
	    tsm[iz][iy] += ts;
	    zsm[iz][iy] += zs;
	    usm[iz][iy] += us;
	    vsm[iz][iy] += vs;
	    pblm[iz][iy] += pbl;
	    pctm[iz][iy] += pct;
	    pcbm[iz][iy] += pcb;
	    clm[iz][iy] += cl;
	    if (gsl_finite(plfc) && gsl_finite(pel) && cape >= ctl.conv_cape
		&& (ctl.conv_cin <= 0 || cin < ctl.conv_cin)) {
	      plclm[iz][iy] += plcl;
	      plfcm[iz][iy] += plfc;
	      pelm[iz][iy] += pel;
	      capem[iz][iy] += cape;
	      cinm[iz][iy] += cin;
	      npc[iz][iy]++;
	    }
	    if (gsl_finite(pt)) {
	      ptm[iz][iy] += pt;
	      ztm[iz][iy] += zt;
	      ttm[iz][iy] += tt;
	      h2otm[iz][iy] += h2ot;
	      npt[iz][iy]++;
	    }
	    rhm[iz][iy] += RH(plev[iz], t, h2o);
	    rhicem[iz][iy] += RHICE(plev[iz], t, h2o);
	    tdewm[iz][iy] += TDEW(plev[iz], h2o);
	    ticem[iz][iy] += TICE(plev[iz], h2o);
	    hno3m[iz][iy] += clim_hno3(met->time, lats[iy], plev[iz], clim);
	    tnatm[iz][iy] +=
	      nat_temperature(plev[iz], h2o,
			      clim_hno3(met->time, lats[iy], plev[iz], clim));
	    ohm[iz][iy] +=
	      clim_oh_diurnal(&ctl, clim, met->time, plev[iz], met->lon[ix],
			      lats[iy]);
	    np[iz][iy]++;
	  }
  }

  /* Create output file... */
  LOG(1, "Write meteorological data file: %s", argv[2]);
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
	  "# $27 = cloud bottom pressure [hPa]\n"
	  "# $28 = pressure at lifted condensation level (LCL) [hPa]\n"
	  "# $29 = pressure at level of free convection (LFC) [hPa]\n"
	  "# $30 = pressure at equilibrium level (EL) [hPa]\n");
  fprintf(out,
	  "# $31 = convective available potential energy (CAPE) [J/kg]\n"
	  "# $32 = convective inhibition (CIN) [J/kg]\n"
	  "# $33 = relative humidity over water [%%]\n"
	  "# $34 = relative humidity over ice [%%]\n"
	  "# $35 = dew point temperature [K]\n"
	  "# $36 = frost point temperature [K]\n"
	  "# $37 = NAT temperature [K]\n"
	  "# $38 = HNO3 volume mixing ratio [ppv]\n"
	  "# $39 = OH concentration [molec/cm^3]\n"
	  "# $40 = boundary layer pressure [hPa]\n");
  fprintf(out,
	  "# $41 = number of data points\n"
	  "# $42 = number of tropopause data points\n"
	  "# $43 = number of CAPE data points\n");

  /* Write data... */
  for (iz = 0; iz < nz; iz++) {
    fprintf(out, "\n");
    for (iy = 0; iy < ny; iy++)
      fprintf(out,
	      "%.2f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g"
	      "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g"
	      " %g %g %g %g %g %d %d %d\n",
	      timem[iz][iy] / np[iz][iy], Z(plev[iz]),
	      lonm[iz][iy] / np[iz][iy], lats[iy],
	      plev[iz], tm[iz][iy] / np[iz][iy], um[iz][iy] / np[iz][iy],
	      vm[iz][iy] / np[iz][iy], wm[iz][iy] / np[iz][iy],
	      h2om[iz][iy] / np[iz][iy], o3m[iz][iy] / np[iz][iy],
	      zm[iz][iy] / np[iz][iy], pvm[iz][iy] / np[iz][iy],
	      psm[iz][iy] / np[iz][iy], tsm[iz][iy] / np[iz][iy],
	      zsm[iz][iy] / np[iz][iy], usm[iz][iy] / np[iz][iy],
	      vsm[iz][iy] / np[iz][iy], ptm[iz][iy] / npt[iz][iy],
	      ztm[iz][iy] / npt[iz][iy], ttm[iz][iy] / npt[iz][iy],
	      h2otm[iz][iy] / npt[iz][iy], lwcm[iz][iy] / np[iz][iy],
	      iwcm[iz][iy] / np[iz][iy], clm[iz][iy] / np[iz][iy],
	      pctm[iz][iy] / np[iz][iy], pcbm[iz][iy] / np[iz][iy],
	      plclm[iz][iy] / npc[iz][iy], plfcm[iz][iy] / npc[iz][iy],
	      pelm[iz][iy] / npc[iz][iy], capem[iz][iy] / npc[iz][iy],
	      cinm[iz][iy] / npc[iz][iy], rhm[iz][iy] / np[iz][iy],
	      rhicem[iz][iy] / np[iz][iy], tdewm[iz][iy] / np[iz][iy],
	      ticem[iz][iy] / np[iz][iy], tnatm[iz][iy] / np[iz][iy],
	      hno3m[iz][iy] / np[iz][iy], ohm[iz][iy] / np[iz][iy],
	      pblm[iz][iy] / np[iz][iy], np[iz][iy],
	      npt[iz][iy], npc[iz][iy]);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(met);
  free(clim);

  return EXIT_SUCCESS;
}
