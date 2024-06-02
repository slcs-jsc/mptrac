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
  Calculate air parcel statistics.
*/

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  atm_t *atm, *atm_filt;

  FILE *out;

  double latm, lonm, t, t0 = NAN, qm[NQ], *work, zm, *zs;

  int f, init = 0, ip, iq;

  /* Allocate... */
  ALLOC(atm, atm_t, 1);
  ALLOC(atm_filt, atm_t, 1);
  ALLOC(work, double,
	NP);
  ALLOC(zs, double,
	NP);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <stat.tab> <param> <atm1> [<atm2> ...]");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  int ens = (int) scan_ctl(argv[1], argc, argv, "STAT_ENS", -1, "-999", NULL);
  double p0 = P(scan_ctl(argv[1], argc, argv, "STAT_Z0", -1, "-1000", NULL));
  double p1 = P(scan_ctl(argv[1], argc, argv, "STAT_Z1", -1, "1000", NULL));
  double lat0 = scan_ctl(argv[1], argc, argv, "STAT_LAT0", -1, "-1000", NULL);
  double lat1 = scan_ctl(argv[1], argc, argv, "STAT_LAT1", -1, "1000", NULL);
  double lon0 = scan_ctl(argv[1], argc, argv, "STAT_LON0", -1, "-1000", NULL);
  double lon1 = scan_ctl(argv[1], argc, argv, "STAT_LON1", -1, "1000", NULL);

  /* Write info... */
  LOG(1, "Write air parcel statistics: %s", argv[2]);

  /* Create output file... */
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = time difference [s]\n"
	  "# $3 = altitude (%s) [km]\n"
	  "# $4 = longitude (%s) [deg]\n"
	  "# $5 = latitude (%s) [deg]\n", argv[3], argv[3], argv[3]);
  for (iq = 0; iq < ctl.nq; iq++)
    fprintf(out, "# $%d = %s (%s) [%s]\n", iq + 6,
	    ctl.qnt_name[iq], argv[3], ctl.qnt_unit[iq]);
  fprintf(out, "# $%d = number of particles\n\n", ctl.nq + 6);

  /* Loop over files... */
  for (f = 4; f < argc; f++) {

    /* Read atmopheric data... */
    if (!read_atm(argv[f], &ctl, atm))
      continue;

    /* Get time from filename... */
    t = time_from_filename(argv[f], ctl.atm_type < 2 ? 20 : 19);

    /* Save initial time... */
    if (!init) {
      init = 1;
      t0 = t;
    }

    /* Filter data... */
    atm_filt->np = 0;
    for (ip = 0; ip < atm->np; ip++) {

      /* Check time... */
      if (!isfinite(atm->time[ip]))
	continue;

      /* Check ensemble index... */
      if (ctl.qnt_ens > 0 && atm->q[ctl.qnt_ens][ip] != ens)
	continue;

      /* Check spatial range... */
      if (atm->p[ip] > p0 || atm->p[ip] < p1
	  || atm->lon[ip] < lon0 || atm->lon[ip] > lon1
	  || atm->lat[ip] < lat0 || atm->lat[ip] > lat1)
	continue;

      /* Save data... */
      atm_filt->time[atm_filt->np] = atm->time[ip];
      atm_filt->p[atm_filt->np] = atm->p[ip];
      atm_filt->lon[atm_filt->np] = atm->lon[ip];
      atm_filt->lat[atm_filt->np] = atm->lat[ip];
      for (iq = 0; iq < ctl.nq; iq++)
	atm_filt->q[iq][atm_filt->np] = atm->q[iq][ip];
      atm_filt->np++;
    }

    /* Get heights... */
    for (ip = 0; ip < atm_filt->np; ip++)
      zs[ip] = Z(atm_filt->p[ip]);

    /* Get statistics... */
    if (strcasecmp(argv[3], "mean") == 0) {
      zm = gsl_stats_mean(zs, 1, (size_t) atm_filt->np);
      lonm = gsl_stats_mean(atm_filt->lon, 1, (size_t) atm_filt->np);
      latm = gsl_stats_mean(atm_filt->lat, 1, (size_t) atm_filt->np);
      for (iq = 0; iq < ctl.nq; iq++)
	qm[iq] = gsl_stats_mean(atm_filt->q[iq], 1, (size_t) atm_filt->np);
    } else if (strcasecmp(argv[3], "stddev") == 0) {
      zm = gsl_stats_sd(zs, 1, (size_t) atm_filt->np);
      lonm = gsl_stats_sd(atm_filt->lon, 1, (size_t) atm_filt->np);
      latm = gsl_stats_sd(atm_filt->lat, 1, (size_t) atm_filt->np);
      for (iq = 0; iq < ctl.nq; iq++)
	qm[iq] = gsl_stats_sd(atm_filt->q[iq], 1, (size_t) atm_filt->np);
    } else if (strcasecmp(argv[3], "min") == 0) {
      zm = gsl_stats_min(zs, 1, (size_t) atm_filt->np);
      lonm = gsl_stats_min(atm_filt->lon, 1, (size_t) atm_filt->np);
      latm = gsl_stats_min(atm_filt->lat, 1, (size_t) atm_filt->np);
      for (iq = 0; iq < ctl.nq; iq++)
	qm[iq] = gsl_stats_min(atm_filt->q[iq], 1, (size_t) atm_filt->np);
    } else if (strcasecmp(argv[3], "max") == 0) {
      zm = gsl_stats_max(zs, 1, (size_t) atm_filt->np);
      lonm = gsl_stats_max(atm_filt->lon, 1, (size_t) atm_filt->np);
      latm = gsl_stats_max(atm_filt->lat, 1, (size_t) atm_filt->np);
      for (iq = 0; iq < ctl.nq; iq++)
	qm[iq] = gsl_stats_max(atm_filt->q[iq], 1, (size_t) atm_filt->np);
    } else if (strcasecmp(argv[3], "skew") == 0) {
      zm = gsl_stats_skew(zs, 1, (size_t) atm_filt->np);
      lonm = gsl_stats_skew(atm_filt->lon, 1, (size_t) atm_filt->np);
      latm = gsl_stats_skew(atm_filt->lat, 1, (size_t) atm_filt->np);
      for (iq = 0; iq < ctl.nq; iq++)
	qm[iq] = gsl_stats_skew(atm_filt->q[iq], 1, (size_t) atm_filt->np);
    } else if (strcasecmp(argv[3], "kurt") == 0) {
      zm = gsl_stats_kurtosis(zs, 1, (size_t) atm_filt->np);
      lonm = gsl_stats_kurtosis(atm_filt->lon, 1, (size_t) atm_filt->np);
      latm = gsl_stats_kurtosis(atm_filt->lat, 1, (size_t) atm_filt->np);
      for (iq = 0; iq < ctl.nq; iq++)
	qm[iq] =
	  gsl_stats_kurtosis(atm_filt->q[iq], 1, (size_t) atm_filt->np);
    } else if (strcasecmp(argv[3], "median") == 0) {
      zm = gsl_stats_median(zs, 1, (size_t) atm_filt->np);
      lonm = gsl_stats_median(atm_filt->lon, 1, (size_t) atm_filt->np);
      latm = gsl_stats_median(atm_filt->lat, 1, (size_t) atm_filt->np);
      for (iq = 0; iq < ctl.nq; iq++)
	qm[iq] = gsl_stats_median(atm_filt->q[iq], 1, (size_t) atm_filt->np);
    } else if (strcasecmp(argv[3], "absdev") == 0) {
      zm = gsl_stats_absdev(zs, 1, (size_t) atm_filt->np);
      lonm = gsl_stats_absdev(atm_filt->lon, 1, (size_t) atm_filt->np);
      latm = gsl_stats_absdev(atm_filt->lat, 1, (size_t) atm_filt->np);
      for (iq = 0; iq < ctl.nq; iq++)
	qm[iq] = gsl_stats_absdev(atm_filt->q[iq], 1, (size_t) atm_filt->np);
    } else if (strcasecmp(argv[3], "mad") == 0) {
      zm = gsl_stats_mad0(zs, 1, (size_t) atm_filt->np, work);
      lonm = gsl_stats_mad0(atm_filt->lon, 1, (size_t) atm_filt->np, work);
      latm = gsl_stats_mad0(atm_filt->lat, 1, (size_t) atm_filt->np, work);
      for (iq = 0; iq < ctl.nq; iq++)
	qm[iq] =
	  gsl_stats_mad0(atm_filt->q[iq], 1, (size_t) atm_filt->np, work);
    } else
      ERRMSG("Unknown parameter!");

    /* Write data... */
    fprintf(out, "%.2f %.2f %g %g %g", t, t - t0, zm, lonm, latm);
    for (iq = 0; iq < ctl.nq; iq++) {
      fprintf(out, " ");
      fprintf(out, ctl.qnt_format[iq], qm[iq]);
    }
    fprintf(out, " %d\n", atm_filt->np);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(atm);
  free(atm_filt);
  free(work);
  free(zs);

  return EXIT_SUCCESS;
}
