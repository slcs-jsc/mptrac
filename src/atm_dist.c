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
  
  Copyright (C) 2013-2022 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Calculate transport deviations of trajectories.
*/

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  atm_t *atm1, *atm2;

  FILE *out;

  double *ahtd, *aqtd, *avtd, ahtdm, aqtdm[NQ], avtdm, lat0, lat1,
    *lat1_old, *lat2_old, *lh1, *lh2, lon0, lon1, *lon1_old, *lon2_old,
    *lv1, *lv2, p0, p1, *rhtd, *rqtd, *rvtd, rhtdm, rqtdm[NQ], rvtdm,
    t, t0 = 0, x0[3], x1[3], x2[3], z1, *z1_old, z2, *z2_old, *work, zscore;

  int ens, f, init = 0, ip, iq, np;

  /* Allocate... */
  ALLOC(atm1, atm_t, 1);
  ALLOC(atm2, atm_t, 1);
  ALLOC(lon1_old, double,
	NP);
  ALLOC(lat1_old, double,
	NP);
  ALLOC(z1_old, double,
	NP);
  ALLOC(lh1, double,
	NP);
  ALLOC(lv1, double,
	NP);
  ALLOC(lon2_old, double,
	NP);
  ALLOC(lat2_old, double,
	NP);
  ALLOC(z2_old, double,
	NP);
  ALLOC(lh2, double,
	NP);
  ALLOC(lv2, double,
	NP);
  ALLOC(ahtd, double,
	NP);
  ALLOC(avtd, double,
	NP);
  ALLOC(aqtd, double,
	NP * NQ);
  ALLOC(rhtd, double,
	NP);
  ALLOC(rvtd, double,
	NP);
  ALLOC(rqtd, double,
	NP * NQ);
  ALLOC(work, double,
	NP);

  /* Check arguments... */
  if (argc < 6)
    ERRMSG("Give parameters: <ctl> <dist.tab> <param> <atm1a> <atm1b>"
	   " [<atm2a> <atm2b> ...]");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  ens = (int) scan_ctl(argv[1], argc, argv, "DIST_ENS", -1, "-999", NULL);
  p0 = P(scan_ctl(argv[1], argc, argv, "DIST_Z0", -1, "-1000", NULL));
  p1 = P(scan_ctl(argv[1], argc, argv, "DIST_Z1", -1, "1000", NULL));
  lat0 = scan_ctl(argv[1], argc, argv, "DIST_LAT0", -1, "-1000", NULL);
  lat1 = scan_ctl(argv[1], argc, argv, "DIST_LAT1", -1, "1000", NULL);
  lon0 = scan_ctl(argv[1], argc, argv, "DIST_LON0", -1, "-1000", NULL);
  lon1 = scan_ctl(argv[1], argc, argv, "DIST_LON1", -1, "1000", NULL);
  zscore = scan_ctl(argv[1], argc, argv, "DIST_ZSCORE", -1, "-999", NULL);

  /* Write info... */
  LOG(1, "Write transport deviations: %s", argv[2]);

  /* Create output file... */
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = time difference [s]\n"
	  "# $3 = absolute horizontal distance (%s) [km]\n"
	  "# $4 = relative horizontal distance (%s) [%%]\n"
	  "# $5 = absolute vertical distance (%s) [km]\n"
	  "# $6 = relative vertical distance (%s) [%%]\n",
	  argv[3], argv[3], argv[3], argv[3]);
  for (iq = 0; iq < ctl.nq; iq++)
    fprintf(out,
	    "# $%d = %s absolute difference (%s) [%s]\n"
	    "# $%d = %s relative difference (%s) [%%]\n",
	    7 + 2 * iq, ctl.qnt_name[iq], argv[3], ctl.qnt_unit[iq],
	    8 + 2 * iq, ctl.qnt_name[iq], argv[3]);
  fprintf(out, "# $%d = number of particles\n\n", 7 + 2 * ctl.nq);

  /* Loop over file pairs... */
  for (f = 4; f < argc; f += 2) {

    /* Read atmopheric data... */
    if (!read_atm(argv[f], &ctl, atm1) || !read_atm(argv[f + 1], &ctl, atm2))
      continue;

    /* Check if structs match... */
    if (atm1->np != atm2->np)
      ERRMSG("Different numbers of particles!");

    /* Get time from filename... */
    t = time_from_filename(argv[f], ctl.atm_type < 2 ? 20 : 19);

    /* Save initial time... */
    if (!init) {
      init = 1;
      t0 = t;
    }

    /* Init... */
    np = 0;
    for (ip = 0; ip < atm1->np; ip++) {
      ahtd[ip] = avtd[ip] = rhtd[ip] = rvtd[ip] = 0;
      for (iq = 0; iq < ctl.nq; iq++)
	aqtd[iq * NP + ip] = rqtd[iq * NP + ip] = 0;
    }

    /* Loop over air parcels... */
    for (ip = 0; ip < atm1->np; ip++) {

      /* Check air parcel index... */
      if (ctl.qnt_idx > 0
	  && (atm1->q[ctl.qnt_idx][ip] != atm2->q[ctl.qnt_idx][ip]))
	ERRMSG("Air parcel index does not match!");

      /* Check ensemble index... */
      if (ctl.qnt_ens > 0
	  && (atm1->q[ctl.qnt_ens][ip] != ens
	      || atm2->q[ctl.qnt_ens][ip] != ens))
	continue;

      /* Check time... */
      if (!gsl_finite(atm1->time[ip]) || !gsl_finite(atm2->time[ip]))
	continue;

      /* Check spatial range... */
      if (atm1->p[ip] > p0 || atm1->p[ip] < p1
	  || atm1->lon[ip] < lon0 || atm1->lon[ip] > lon1
	  || atm1->lat[ip] < lat0 || atm1->lat[ip] > lat1)
	continue;
      if (atm2->p[ip] > p0 || atm2->p[ip] < p1
	  || atm2->lon[ip] < lon0 || atm2->lon[ip] > lon1
	  || atm2->lat[ip] < lat0 || atm2->lat[ip] > lat1)
	continue;

      /* Convert coordinates... */
      geo2cart(0, atm1->lon[ip], atm1->lat[ip], x1);
      geo2cart(0, atm2->lon[ip], atm2->lat[ip], x2);
      z1 = Z(atm1->p[ip]);
      z2 = Z(atm2->p[ip]);

      /* Calculate absolute transport deviations... */
      ahtd[np] = DIST(x1, x2);
      avtd[np] = z1 - z2;
      for (iq = 0; iq < ctl.nq; iq++)
	aqtd[iq * NP + np] = atm1->q[iq][ip] - atm2->q[iq][ip];

      /* Calculate relative transport deviations... */
      if (f > 4) {

	/* Get trajectory lengths... */
	geo2cart(0, lon1_old[ip], lat1_old[ip], x0);
	lh1[ip] += DIST(x0, x1);
	lv1[ip] += fabs(z1_old[ip] - z1);

	geo2cart(0, lon2_old[ip], lat2_old[ip], x0);
	lh2[ip] += DIST(x0, x2);
	lv2[ip] += fabs(z2_old[ip] - z2);

	/* Get relative transport deviations... */
	if (lh1[ip] + lh2[ip] > 0)
	  rhtd[np] = 200. * DIST(x1, x2) / (lh1[ip] + lh2[ip]);
	if (lv1[ip] + lv2[ip] > 0)
	  rvtd[np] = 200. * (z1 - z2) / (lv1[ip] + lv2[ip]);
      }

      /* Get relative transport deviations... */
      for (iq = 0; iq < ctl.nq; iq++)
	rqtd[iq * NP + np] = 200. * (atm1->q[iq][ip] - atm2->q[iq][ip])
	  / (fabs(atm1->q[iq][ip]) + fabs(atm2->q[iq][ip]));

      /* Save positions of air parcels... */
      lon1_old[ip] = atm1->lon[ip];
      lat1_old[ip] = atm1->lat[ip];
      z1_old[ip] = z1;

      lon2_old[ip] = atm2->lon[ip];
      lat2_old[ip] = atm2->lat[ip];
      z2_old[ip] = z2;

      /* Increment air parcel counter... */
      np++;
    }

    /* Filter data... */
    if (zscore > 0 && np > 1) {

      /* Get means and standard deviations of transport deviations... */
      size_t n = (size_t) np;
      double muh = gsl_stats_mean(ahtd, 1, n);
      double muv = gsl_stats_mean(avtd, 1, n);
      double sigh = gsl_stats_sd(ahtd, 1, n);
      double sigv = gsl_stats_sd(avtd, 1, n);

      /* Filter data... */
      np = 0;
      for (size_t i = 0; i < n; i++)
	if (fabs((ahtd[i] - muh) / sigh) < zscore
	    && fabs((avtd[i] - muv) / sigv) < zscore) {
	  ahtd[np] = ahtd[i];
	  rhtd[np] = rhtd[i];
	  avtd[np] = avtd[i];
	  rvtd[np] = rvtd[i];
	  for (iq = 0; iq < ctl.nq; iq++) {
	    aqtd[iq * NP + np] = aqtd[iq * NP + (int) i];
	    rqtd[iq * NP + np] = rqtd[iq * NP + (int) i];
	  }
	  np++;
	}
    }

    /* Get statistics... */
    if (strcasecmp(argv[3], "mean") == 0) {
      ahtdm = gsl_stats_mean(ahtd, 1, (size_t) np);
      rhtdm = gsl_stats_mean(rhtd, 1, (size_t) np);
      avtdm = gsl_stats_mean(avtd, 1, (size_t) np);
      rvtdm = gsl_stats_mean(rvtd, 1, (size_t) np);
      for (iq = 0; iq < ctl.nq; iq++) {
	aqtdm[iq] = gsl_stats_mean(&aqtd[iq * NP], 1, (size_t) np);
	rqtdm[iq] = gsl_stats_mean(&rqtd[iq * NP], 1, (size_t) np);
      }
    } else if (strcasecmp(argv[3], "stddev") == 0) {
      ahtdm = gsl_stats_sd(ahtd, 1, (size_t) np);
      rhtdm = gsl_stats_sd(rhtd, 1, (size_t) np);
      avtdm = gsl_stats_sd(avtd, 1, (size_t) np);
      rvtdm = gsl_stats_sd(rvtd, 1, (size_t) np);
      for (iq = 0; iq < ctl.nq; iq++) {
	aqtdm[iq] = gsl_stats_sd(&aqtd[iq * NP], 1, (size_t) np);
	rqtdm[iq] = gsl_stats_sd(&rqtd[iq * NP], 1, (size_t) np);
      }
    } else if (strcasecmp(argv[3], "min") == 0) {
      ahtdm = gsl_stats_min(ahtd, 1, (size_t) np);
      rhtdm = gsl_stats_min(rhtd, 1, (size_t) np);
      avtdm = gsl_stats_min(avtd, 1, (size_t) np);
      rvtdm = gsl_stats_min(rvtd, 1, (size_t) np);
      for (iq = 0; iq < ctl.nq; iq++) {
	aqtdm[iq] = gsl_stats_min(&aqtd[iq * NP], 1, (size_t) np);
	rqtdm[iq] = gsl_stats_min(&rqtd[iq * NP], 1, (size_t) np);
      }
    } else if (strcasecmp(argv[3], "max") == 0) {
      ahtdm = gsl_stats_max(ahtd, 1, (size_t) np);
      rhtdm = gsl_stats_max(rhtd, 1, (size_t) np);
      avtdm = gsl_stats_max(avtd, 1, (size_t) np);
      rvtdm = gsl_stats_max(rvtd, 1, (size_t) np);
      for (iq = 0; iq < ctl.nq; iq++) {
	aqtdm[iq] = gsl_stats_max(&aqtd[iq * NP], 1, (size_t) np);
	rqtdm[iq] = gsl_stats_max(&rqtd[iq * NP], 1, (size_t) np);
      }
    } else if (strcasecmp(argv[3], "skew") == 0) {
      ahtdm = gsl_stats_skew(ahtd, 1, (size_t) np);
      rhtdm = gsl_stats_skew(rhtd, 1, (size_t) np);
      avtdm = gsl_stats_skew(avtd, 1, (size_t) np);
      rvtdm = gsl_stats_skew(rvtd, 1, (size_t) np);
      for (iq = 0; iq < ctl.nq; iq++) {
	aqtdm[iq] = gsl_stats_skew(&aqtd[iq * NP], 1, (size_t) np);
	rqtdm[iq] = gsl_stats_skew(&rqtd[iq * NP], 1, (size_t) np);
      }
    } else if (strcasecmp(argv[3], "kurt") == 0) {
      ahtdm = gsl_stats_kurtosis(ahtd, 1, (size_t) np);
      rhtdm = gsl_stats_kurtosis(rhtd, 1, (size_t) np);
      avtdm = gsl_stats_kurtosis(avtd, 1, (size_t) np);
      rvtdm = gsl_stats_kurtosis(rvtd, 1, (size_t) np);
      for (iq = 0; iq < ctl.nq; iq++) {
	aqtdm[iq] = gsl_stats_kurtosis(&aqtd[iq * NP], 1, (size_t) np);
	rqtdm[iq] = gsl_stats_kurtosis(&rqtd[iq * NP], 1, (size_t) np);
      }
    } else if (strcasecmp(argv[3], "absdev") == 0) {
      ahtdm = gsl_stats_absdev_m(ahtd, 1, (size_t) np, 0.0);
      rhtdm = gsl_stats_absdev_m(rhtd, 1, (size_t) np, 0.0);
      avtdm = gsl_stats_absdev_m(avtd, 1, (size_t) np, 0.0);
      rvtdm = gsl_stats_absdev_m(rvtd, 1, (size_t) np, 0.0);
      for (iq = 0; iq < ctl.nq; iq++) {
	aqtdm[iq] = gsl_stats_absdev_m(&aqtd[iq * NP], 1, (size_t) np, 0.0);
	rqtdm[iq] = gsl_stats_absdev_m(&rqtd[iq * NP], 1, (size_t) np, 0.0);
      }
    } else if (strcasecmp(argv[3], "median") == 0) {
      ahtdm = gsl_stats_median(ahtd, 1, (size_t) np);
      rhtdm = gsl_stats_median(rhtd, 1, (size_t) np);
      avtdm = gsl_stats_median(avtd, 1, (size_t) np);
      rvtdm = gsl_stats_median(rvtd, 1, (size_t) np);
      for (iq = 0; iq < ctl.nq; iq++) {
	aqtdm[iq] = gsl_stats_median(&aqtd[iq * NP], 1, (size_t) np);
	rqtdm[iq] = gsl_stats_median(&rqtd[iq * NP], 1, (size_t) np);
      }
    } else if (strcasecmp(argv[3], "mad") == 0) {
      ahtdm = gsl_stats_mad0(ahtd, 1, (size_t) np, work);
      rhtdm = gsl_stats_mad0(rhtd, 1, (size_t) np, work);
      avtdm = gsl_stats_mad0(avtd, 1, (size_t) np, work);
      rvtdm = gsl_stats_mad0(rvtd, 1, (size_t) np, work);
      for (iq = 0; iq < ctl.nq; iq++) {
	aqtdm[iq] = gsl_stats_mad0(&aqtd[iq * NP], 1, (size_t) np, work);
	rqtdm[iq] = gsl_stats_mad0(&rqtd[iq * NP], 1, (size_t) np, work);
      }
    } else
      ERRMSG("Unknown parameter!");

    /* Write output... */
    fprintf(out, "%.2f %.2f %g %g %g %g", t, t - t0,
	    ahtdm, rhtdm, avtdm, rvtdm);
    for (iq = 0; iq < ctl.nq; iq++) {
      fprintf(out, " ");
      fprintf(out, ctl.qnt_format[iq], aqtdm[iq]);
      fprintf(out, " ");
      fprintf(out, ctl.qnt_format[iq], rqtdm[iq]);
    }
    fprintf(out, " %d\n", np);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(atm1);
  free(atm2);
  free(lon1_old);
  free(lat1_old);
  free(z1_old);
  free(lh1);
  free(lv1);
  free(lon2_old);
  free(lat2_old);
  free(z2_old);
  free(lh2);
  free(lv2);
  free(ahtd);
  free(avtd);
  free(aqtd);
  free(rhtd);
  free(rvtd);
  free(rqtd);
  free(work);

  return EXIT_SUCCESS;
}
