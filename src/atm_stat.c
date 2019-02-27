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
  Calculate air parcel statistics.
*/

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  atm_t *atm;

  FILE *out;

  char tstr[LEN];

  double latm, lonm, t, qm[NQ], *work, zm, *zs;

  int f, ip, iq, year, mon, day, hour, min;

  /* Allocate... */
  ALLOC(atm, atm_t, 1);
  ALLOC(work, double,
	NP);
  ALLOC(zs, double,
	NP);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <outfile> <param> <atm1> [<atm2> ...]");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);

  /* Write info... */
  printf("Write air parcel statistics: %s\n", argv[2]);

  /* Create output file... */
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = altitude (%s) [km]\n"
	  "# $3 = longitude (%s) [deg]\n"
	  "# $4 = latitude (%s) [deg]\n", argv[3], argv[3], argv[3]);
  for (iq = 0; iq < ctl.nq; iq++)
    fprintf(out, "# $%d = %s (%s) [%s]\n", iq + 5,
	    ctl.qnt_name[iq], argv[3], ctl.qnt_unit[iq]);
  fprintf(out, "# $%d = number of particles\n\n", ctl.nq + 5);

  /* Loop over files... */
  for (f = 4; f < argc; f++) {

    /* Read atmopheric data... */
    read_atm(argv[f], &ctl, atm);

    /* Get time from filename... */
    sprintf(tstr, "%.4s", &argv[f][strlen(argv[f]) - 20]);
    year = atoi(tstr);
    sprintf(tstr, "%.2s", &argv[f][strlen(argv[f]) - 15]);
    mon = atoi(tstr);
    sprintf(tstr, "%.2s", &argv[f][strlen(argv[f]) - 12]);
    day = atoi(tstr);
    sprintf(tstr, "%.2s", &argv[f][strlen(argv[f]) - 9]);
    hour = atoi(tstr);
    sprintf(tstr, "%.2s", &argv[f][strlen(argv[f]) - 6]);
    min = atoi(tstr);
    time2jsec(year, mon, day, hour, min, 0, 0, &t);

    /* Get heights... */
    for (ip = 0; ip < atm->np; ip++)
      zs[ip] = Z(atm->p[ip]);

    /* Get statistics... */
    if (strcasecmp(argv[3], "mean") == 0) {
      zm = gsl_stats_mean(zs, 1, (size_t) atm->np);
      lonm = gsl_stats_mean(atm->lon, 1, (size_t) atm->np);
      latm = gsl_stats_mean(atm->lat, 1, (size_t) atm->np);
      for (iq = 0; iq < ctl.nq; iq++)
	qm[iq] = gsl_stats_mean(atm->q[iq], 1, (size_t) atm->np);
    } else if (strcasecmp(argv[3], "stddev") == 0) {
      zm = gsl_stats_sd(zs, 1, (size_t) atm->np);
      lonm = gsl_stats_sd(atm->lon, 1, (size_t) atm->np);
      latm = gsl_stats_sd(atm->lat, 1, (size_t) atm->np);
      for (iq = 0; iq < ctl.nq; iq++)
	qm[iq] = gsl_stats_sd(atm->q[iq], 1, (size_t) atm->np);
    } else if (strcasecmp(argv[3], "min") == 0) {
      zm = gsl_stats_min(zs, 1, (size_t) atm->np);
      lonm = gsl_stats_min(atm->lon, 1, (size_t) atm->np);
      latm = gsl_stats_min(atm->lat, 1, (size_t) atm->np);
      for (iq = 0; iq < ctl.nq; iq++)
	qm[iq] = gsl_stats_min(atm->q[iq], 1, (size_t) atm->np);
    } else if (strcasecmp(argv[3], "max") == 0) {
      zm = gsl_stats_max(zs, 1, (size_t) atm->np);
      lonm = gsl_stats_max(atm->lon, 1, (size_t) atm->np);
      latm = gsl_stats_max(atm->lat, 1, (size_t) atm->np);
      for (iq = 0; iq < ctl.nq; iq++)
	qm[iq] = gsl_stats_max(atm->q[iq], 1, (size_t) atm->np);
    } else if (strcasecmp(argv[3], "skew") == 0) {
      zm = gsl_stats_skew(zs, 1, (size_t) atm->np);
      lonm = gsl_stats_skew(atm->lon, 1, (size_t) atm->np);
      latm = gsl_stats_skew(atm->lat, 1, (size_t) atm->np);
      for (iq = 0; iq < ctl.nq; iq++)
	qm[iq] = gsl_stats_skew(atm->q[iq], 1, (size_t) atm->np);
    } else if (strcasecmp(argv[3], "kurt") == 0) {
      zm = gsl_stats_kurtosis(zs, 1, (size_t) atm->np);
      lonm = gsl_stats_kurtosis(atm->lon, 1, (size_t) atm->np);
      latm = gsl_stats_kurtosis(atm->lat, 1, (size_t) atm->np);
      for (iq = 0; iq < ctl.nq; iq++)
	qm[iq] = gsl_stats_kurtosis(atm->q[iq], 1, (size_t) atm->np);
    } else if (strcasecmp(argv[3], "median") == 0) {
      zm = gsl_stats_median(zs, 1, (size_t) atm->np);
      lonm = gsl_stats_median(atm->lon, 1, (size_t) atm->np);
      latm = gsl_stats_median(atm->lat, 1, (size_t) atm->np);
      for (iq = 0; iq < ctl.nq; iq++)
	qm[iq] = gsl_stats_median(atm->q[iq], 1, (size_t) atm->np);
    } else if (strcasecmp(argv[3], "absdev") == 0) {
      zm = gsl_stats_absdev(zs, 1, (size_t) atm->np);
      lonm = gsl_stats_absdev(atm->lon, 1, (size_t) atm->np);
      latm = gsl_stats_absdev(atm->lat, 1, (size_t) atm->np);
      for (iq = 0; iq < ctl.nq; iq++)
	qm[iq] = gsl_stats_absdev(atm->q[iq], 1, (size_t) atm->np);
    } else if (strcasecmp(argv[3], "mad") == 0) {
      zm = gsl_stats_mad0(zs, 1, (size_t) atm->np, work);
      lonm = gsl_stats_mad0(atm->lon, 1, (size_t) atm->np, work);
      latm = gsl_stats_mad0(atm->lat, 1, (size_t) atm->np, work);
      for (iq = 0; iq < ctl.nq; iq++)
	qm[iq] = gsl_stats_mad0(atm->q[iq], 1, (size_t) atm->np, work);
    } else
      ERRMSG("Unknown parameter!");

    /* Write data... */
    fprintf(out, "%.2f %g %g %g", t, zm, lonm, latm);
    for (iq = 0; iq < ctl.nq; iq++) {
      fprintf(out, " ");
      fprintf(out, ctl.qnt_format[iq], qm[iq]);
    }
    fprintf(out, " %d\n", atm->np);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(atm);
  free(work);
  free(zs);

  return EXIT_SUCCESS;
}
