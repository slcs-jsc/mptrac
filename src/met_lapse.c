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
  Calculate lapse rate statistics.
*/

#include "mptrac.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Lapse rate minimum [K/km. */
#define LAPSEMIN -20.0

/*! Lapse rate bin size [K/km]. */
#define DLAPSE 0.1

/*! Maximum number of histogram bins. */
#define IDXMAX 400

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

  static double p2[1000], t[1000], t2[1000], z[1000], z2[1000], lat_mean,
    z_mean;

  static int hist_max[1000], hist_min[1000], hist_mean[1000], hist_sig[1000],
    nhist_max, nhist_min, nhist_mean, nhist_sig, np;

  /* Allocate... */
  ALLOC(clim, clim_t, 1);
  ALLOC(met, met_t, 1);
  ALLOC(dd, dd_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <lapse.tab> <met0> [ <met1> ... ]");

  /* Read control parameters... */
  mptrac_read_ctl(argv[1], argc, argv, &ctl);
  const int dz =
    (int) scan_ctl(argv[1], argc, argv, "LAPSE_DZ", -1, "20", NULL);
  const double lat0 =
    (int) scan_ctl(argv[1], argc, argv, "LAPSE_LAT0", -1, "-90", NULL);
  const double lat1 =
    (int) scan_ctl(argv[1], argc, argv, "LAPSE_LAT1", -1, "90", NULL);
  const double z0 =
    (int) scan_ctl(argv[1], argc, argv, "LAPSE_Z0", -1, "0", NULL);
  const double z1 =
    (int) scan_ctl(argv[1], argc, argv, "LAPSE_Z1", -1, "100", NULL);
  const int intpol =
    (int) scan_ctl(argv[1], argc, argv, "LAPSE_INTPOL", -1, "1", NULL);

  /* Read climatological data... */
  mptrac_read_clim(&ctl, clim);

  /* Loop over files... */
  for (int i = 3; i < argc; i++) {

    /* Read meteorological data... */
    if (!mptrac_read_met(argv[i], &ctl, clim, met, dd))
      continue;

    /* Get altitude and pressure profiles... */
    for (int iz = 0; iz < met->np; iz++)
      z[iz] = Z(met->p[iz]);
    for (int iz = 0; iz <= 250; iz++) {
      z2[iz] = 0.0 + 0.1 * iz;
      p2[iz] = P(z2[iz]);
    }

    /* Loop over grid points... */
    for (int ix = 0; ix < met->nx; ix++)
      for (int iy = 0; iy < met->ny; iy++) {

	/* Check latitude range... */
	if (met->lat[iy] < lat0 || met->lat[iy] > lat1)
	  continue;

	/* Interpolate temperature profile... */
	for (int iz = 0; iz < met->np; iz++)
	  t[iz] = met->t[ix][iy][iz];
	if (intpol == 1)
	  spline(z, t, met->np, z2, t2, 251, ctl.met_tropo_spline);
	else
	  for (int iz = 0; iz <= 250; iz++) {
	    int idx = locate_irr(z, met->np, z2[iz]);
	    t2[iz] = LIN(z[idx], t[idx], z[idx + 1], t[idx + 1], z2[iz]);
	  }

	/* Loop over vertical levels... */
	for (int iz = 0; iz <= 250; iz++) {

	  /* Check height range... */
	  if (z2[iz] < z0 || z2[iz] > z1)
	    continue;

	  /* Check surface pressure... */
	  if (p2[iz] > met->ps[ix][iy])
	    continue;

	  /* Get mean latitude and height... */
	  lat_mean += met->lat[iy];
	  z_mean += z2[iz];
	  np++;

	  /* Get lapse rates within a vertical layer... */
	  int nlapse = 0;
	  double lapse_max = -1e99, lapse_min = 1e99, lapse_mean =
	    0, lapse_sig = 0;
	  for (int iz2 = iz + 1; iz2 <= iz + dz; iz2++) {
	    lapse_max =
	      MAX(LAPSE(p2[iz], t2[iz], p2[iz2], t2[iz2]), lapse_max);
	    lapse_min =
	      MIN(LAPSE(p2[iz], t2[iz], p2[iz2], t2[iz2]), lapse_min);
	    lapse_mean += LAPSE(p2[iz], t2[iz], p2[iz2], t2[iz2]);
	    lapse_sig += SQR(LAPSE(p2[iz], t2[iz], p2[iz2], t2[iz2]));
	    nlapse++;
	  }
	  lapse_mean /= nlapse;
	  lapse_sig = sqrt(MAX(lapse_sig / nlapse - SQR(lapse_mean), 0));

	  /* Get histograms... */
	  int idx = (int) ((lapse_max - LAPSEMIN) / DLAPSE);
	  if (idx >= 0 && idx < IDXMAX) {
	    hist_max[idx]++;
	    nhist_max++;
	  }

	  idx = (int) ((lapse_min - LAPSEMIN) / DLAPSE);
	  if (idx >= 0 && idx < IDXMAX) {
	    hist_min[idx]++;
	    nhist_min++;
	  }

	  idx = (int) ((lapse_mean - LAPSEMIN) / DLAPSE);
	  if (idx >= 0 && idx < IDXMAX) {
	    hist_mean[idx]++;
	    nhist_mean++;
	  }

	  idx = (int) ((lapse_sig - LAPSEMIN) / DLAPSE);
	  if (idx >= 0 && idx < IDXMAX) {
	    hist_sig[idx]++;
	    nhist_sig++;
	  }
	}
      }
  }

  /* Create output file... */
  LOG(1, "Write lapse rate data: %s", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = mean altitude [km]\n"
	  "# $2 = mean latitude [deg]\n"
	  "# $3 = lapse rate [K/km]\n"
	  "# $4 = counts of maxima per bin\n"
	  "# $5 = total number of maxima\n"
	  "# $6 = normalized frequency of maxima\n"
	  "# $7 = counts of minima per bin\n"
	  "# $8 = total number of minima\n"
	  "# $9 = normalized frequency of minima\n"
	  "# $10 = counts of means per bin\n"
	  "# $11 = total number of means\n"
	  "# $12 = normalized frequency of means\n"
	  "# $13 = counts of sigmas per bin\n"
	  "# $14 = total number of sigmas\n"
	  "# $15 = normalized frequency of sigmas\n\n");

  /* Write data... */
  double nmax_max = 0, nmax_min = 0, nmax_mean = 0, nmax_sig = 0;
  for (int idx = 0; idx < IDXMAX; idx++) {
    nmax_max = MAX(hist_max[idx], nmax_max);
    nmax_min = MAX(hist_min[idx], nmax_min);
    nmax_mean = MAX(hist_mean[idx], nmax_mean);
    nmax_sig = MAX(hist_sig[idx], nmax_sig);
  }
  for (int idx = 0; idx < IDXMAX; idx++)
    fprintf(out,
	    "%g %g %g %d %d %g %d %d %g %d %d %g %d %d %g\n",
	    z_mean / np, lat_mean / np, (idx + .5) * DLAPSE + LAPSEMIN,
	    hist_max[idx], nhist_max,
	    (double) hist_max[idx] / (double) nmax_max, hist_min[idx],
	    nhist_min, (double) hist_min[idx] / (double) nmax_min,
	    hist_mean[idx], nhist_mean,
	    (double) hist_mean[idx] / (double) nmax_mean, hist_sig[idx],
	    nhist_sig, (double) hist_sig[idx] / (double) nmax_sig);

  /* Close file... */
  fclose(out);

  /* Free... */
  free(clim);
  free(met);
  free(dd);

  return EXIT_SUCCESS;
}
