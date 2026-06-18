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
  
  Copyright (C) 2013-2026 Forschungszentrum Juelich GmbH
*/

/*!
  \file
  Create meteorological data files with synthetic wind fields.
*/

#include "mptrac.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Print command-line help. */
void usage(
  void);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;
  met_t *met;

  static char filename[LEN];

  static int year, mon, day, hour, min, sec;
  static double r;

  /* Print usage information... */
  USAGE;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: wind <ctl> <metbase> [KEY VALUE ...]\n\n"
	   "Use -h for full help.");

  /* Read control parameters... */
  mptrac_read_ctl(argv[1], argc, argv, &ctl);
  const double t0 = scan_ctl(argv[1], argc, argv, "WIND_T0", -1, "0", NULL);
  const int nx =
    (int) scan_ctl(argv[1], argc, argv, "WIND_NX", -1, "360", NULL);
  const int ny =
    (int) scan_ctl(argv[1], argc, argv, "WIND_NY", -1, "181", NULL);
  const int nz =
    (int) scan_ctl(argv[1], argc, argv, "WIND_NZ", -1, "61", NULL);
  const double z0 = scan_ctl(argv[1], argc, argv, "WIND_Z0", -1, "0", NULL);
  const double z1 = scan_ctl(argv[1], argc, argv, "WIND_Z1", -1, "60", NULL);
  const double u0 =
    scan_ctl(argv[1], argc, argv, "WIND_U0", -1, "38.587660177302", NULL);
  const double u1 =
    scan_ctl(argv[1], argc, argv, "WIND_U1", -1, "38.587660177302", NULL);
  const double w0 = scan_ctl(argv[1], argc, argv, "WIND_W0", -1, "0", NULL);
  const double alpha =
    scan_ctl(argv[1], argc, argv, "WIND_ALPHA", -1, "0.0", NULL);
  const int lat_reverse =
    (int) scan_ctl(argv[1], argc, argv, "WIND_LAT_REVERSE", -1, "0", NULL);
  const double temp0 =
    scan_ctl(argv[1], argc, argv, "WIND_TEMP0", -1, "280", NULL);
  const double temp1 =
    scan_ctl(argv[1], argc, argv, "WIND_TEMP1", -1, "280", NULL);
  const double ps =
    scan_ctl(argv[1], argc, argv, "WIND_PS", -1, "1013.25", NULL);
  const double zs = scan_ctl(argv[1], argc, argv, "WIND_ZS", -1, "0", NULL);
  const double t2m =
    scan_ctl(argv[1], argc, argv, "WIND_T2M", -1, "280", NULL);
  const double iews =
    scan_ctl(argv[1], argc, argv, "WIND_IEWS", -1, "0", NULL);
  const double inss =
    scan_ctl(argv[1], argc, argv, "WIND_INSS", -1, "0", NULL);
  const double ishf =
    scan_ctl(argv[1], argc, argv, "WIND_ISHF", -1, "0", NULL);
  const double lsm = scan_ctl(argv[1], argc, argv, "WIND_LSM", -1, "1", NULL);
  const double sst =
    scan_ctl(argv[1], argc, argv, "WIND_SST", -1, "280", NULL);
  const double blh =
    scan_ctl(argv[1], argc, argv, "WIND_BLH", -1, "1.0", NULL);
  const double q = scan_ctl(argv[1], argc, argv, "WIND_Q", -1, "0", NULL);
  const double o3 = scan_ctl(argv[1], argc, argv, "WIND_O3", -1, "0", NULL);

  /* Check dimensions... */
  if (nx < 2 || nx > EX)
    ERRMSG("Set 2 <= NX <= MAX!");
  if (ny < 2 || ny > EY)
    ERRMSG("Set 2 <= NY <= MAX!");
  if (nz < 2 || nz > EP)
    ERRMSG("Set 2 <= NZ <= MAX!");

  /* Get time and output filename... */
  jsec2time(t0, &year, &mon, &day, &hour, &min, &sec, &r);
  sprintf(filename, "%s_%d_%02d_%02d_%02d.nc", argv[2], year, mon, day, hour);

  /* Initialize synthetic meteorological data. */
  ALLOC(met, met_t, 1);
  met->time = t0;
  met->coord_type = 0;
  met->nx = nx;
  met->ny = ny;
  met->np = nz;

  /* Set grid... */
  for (int ix = 0; ix < nx; ix++)
    met->lon[ix] = 360.0 / nx * (double) ix;
  for (int iy = 0; iy < ny; iy++)
    met->lat[iy] = (lat_reverse ? -(180.0 / (ny - 1) * (double) iy - 90.0)
		    : (180.0 / (ny - 1) * (double) iy - 90.0));
  for (int iz = 0; iz < nz; iz++)
    met->p[iz] = P(LIN(0.0, z0, nz - 1.0, z1, iz));

  /* Set meteo data... */
  for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++) {
      const double usfc = u0
	* (cos(DEG2RAD(met->lat[iy])) * cos(DEG2RAD(alpha))
	   + sin(DEG2RAD(met->lat[iy])) * cos(DEG2RAD(met->lon[ix]))
	   * sin(DEG2RAD(alpha)));
      const double vsfc =
	-u0 * sin(DEG2RAD(met->lon[ix])) * sin(DEG2RAD(alpha));

      met->ps[ix][iy] = (float) ps;
      met->zs[ix][iy] = (float) zs;
      met->ts[ix][iy] = (float) t2m;
      met->us[ix][iy] = (float) usfc;
      met->vs[ix][iy] = (float) vsfc;
      met->ess[ix][iy] = (float) iews;
      met->nss[ix][iy] = (float) inss;
      met->shf[ix][iy] = (float) ishf;
      met->lsm[ix][iy] = (float) lsm;
      met->sst[ix][iy] = (float) sst;
      met->pbl[ix][iy] = (float) P(zs + blh);

      for (int iz = 0; iz < nz; iz++) {
	const double u = LIN(0.0, u0, nz - 1.0, u1, iz)
	  * (cos(DEG2RAD(met->lat[iy])) * cos(DEG2RAD(alpha))
	     + sin(DEG2RAD(met->lat[iy])) * cos(DEG2RAD(met->lon[ix]))
	     * sin(DEG2RAD(alpha)));
	const double v = -LIN(0.0, u0, nz - 1.0, u1, iz)
	  * sin(DEG2RAD(met->lon[ix])) * sin(DEG2RAD(alpha));

	met->t[ix][iy][iz] = (float) LIN(0.0, temp0, nz - 1.0, temp1, iz);
	met->u[ix][iy][iz] = (float) u;
	met->v[ix][iy][iz] = (float) v;
	met->w[ix][iy][iz] = (float) DZ2DP(1e-3 * w0, met->p[iz]);
	met->h2o[ix][iy][iz] = (float) (q * MA / MH2O);
	met->o3[ix][iy][iz] = (float) (o3 * MA / MO3);
      }
    }

  /* Write synthetic meteorological data. */
  write_met_nc(filename, &ctl, met);

  /* Free... */
  free(met);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

/*! Print command-line help. */
void usage(
  void) {

  printf("\nMPTRAC wind tool.\n\n");
  printf("Create meteorological data with synthetic wind fields.\n");
  printf("\n");
  printf("Usage:\n");
  printf("  wind <ctl> <metbase> [KEY VALUE ...]\n");
  printf("\n");
  printf("Arguments:\n");
  printf("  <ctl>      Control file.\n");
  printf("  <metbase>  Basename of the output meteorological data.\n");
  printf("  [KEY VALUE]  Optional control parameters.\n");
  printf("\nFurther information:\n");
  printf("  Manual: https://slcs-jsc.github.io/mptrac/\n");
}
