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
  
  Copright (C) 2013-2018 Forschungszentrum Juelich GmbH
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

  char tstr[LEN];

  double ahtd, aqtd[NQ], atce1[NQ], atce2[NQ], avtd, lat0, lat1,
    *lat1_old, *lat2_old, *lh1, *lh2, lon0, lon1, *lon1_old, *lon2_old,
    *lv1, *lv2, p0, p1, *q1, *q2, rhtd, rqtd[NQ], rtce1[NQ], rtce2[NQ], rvtd,
    t, t0, x0[3], x1[3], x2[3], z1, *z1_old, z2, *z2_old;

  int ens, f, ip, iq, np, year, mon, day, hour, min;

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
  ALLOC(q1, double,
	NQ * NP);
  ALLOC(q2, double,
	NQ * NP);

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <outfile> <atm1a> <atm1b>"
	   " [<atm2a> <atm2b> ...]");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  ens = (int) scan_ctl(argv[1], argc, argv, "DIST_ENS", -1, "-1", NULL);
  p0 = P(scan_ctl(argv[1], argc, argv, "DIST_Z0", -1, "-1000", NULL));
  p1 = P(scan_ctl(argv[1], argc, argv, "DIST_Z1", -1, "1000", NULL));
  lat0 = scan_ctl(argv[1], argc, argv, "DIST_LAT0", -1, "-1000", NULL);
  lat1 = scan_ctl(argv[1], argc, argv, "DIST_LAT1", -1, "1000", NULL);
  lon0 = scan_ctl(argv[1], argc, argv, "DIST_LON0", -1, "-1000", NULL);
  lon1 = scan_ctl(argv[1], argc, argv, "DIST_LON1", -1, "1000", NULL);

  /* Write info... */
  printf("Write transport deviations: %s\n", argv[2]);

  /* Create output file... */
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = trajectory time [s]\n"
	  "# $3 = AHTD [km]\n"
	  "# $4 = RHTD [%%]\n" "# $5 = AVTD [km]\n" "# $6 = RVTD [%%]\n");
  for (iq = 0; iq < ctl.nq; iq++)
    fprintf(out,
	    "# $%d = AQTD (%s) [%s]\n"
	    "# $%d = RQTD (%s) [%%]\n",
	    7 + 2 * iq, ctl.qnt_name[iq], ctl.qnt_unit[iq],
	    8 + 2 * iq, ctl.qnt_name[iq]);
  for (iq = 0; iq < ctl.nq; iq++)
    fprintf(out,
	    "# $%d = ATCE_1 (%s) [%s]\n"
	    "# $%d = RTCE_1 (%s) [%%]\n",
	    7 + 2 * ctl.nq + 2 * iq, ctl.qnt_name[iq], ctl.qnt_unit[iq],
	    8 + 2 * ctl.nq + 2 * iq, ctl.qnt_name[iq]);
  for (iq = 0; iq < ctl.nq; iq++)
    fprintf(out,
	    "# $%d = ATCE_2 (%s) [%s]\n"
	    "# $%d = RTCE_2 (%s) [%%]\n",
	    7 + 4 * ctl.nq + 2 * iq, ctl.qnt_name[iq], ctl.qnt_unit[iq],
	    8 + 4 * ctl.nq + 2 * iq, ctl.qnt_name[iq]);
  fprintf(out, "# $%d = number of particles\n\n", 7 + 6 * ctl.nq);

  /* Loop over file pairs... */
  for (f = 3; f < argc; f += 2) {

    /* Read atmopheric data... */
    read_atm(argv[f], &ctl, atm1);
    read_atm(argv[f + 1], &ctl, atm2);

    /* Check if structs match... */
    if (atm1->np != atm2->np)
      ERRMSG("Different numbers of parcels!");
    for (ip = 0; ip < atm1->np; ip++)
      if (gsl_finite(atm1->time[ip]) && gsl_finite(atm2->time[ip])
	  && atm1->time[ip] != atm2->time[ip])
	ERRMSG("Times do not match!");
    
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

    /* Save initial data... */
    if (f == 3) {
      t0 = t;
      for (iq = 0; iq < ctl.nq; iq++)
	for (ip = 0; ip < atm1->np; ip++) {
	  q1[iq * NP + ip] = atm1->q[iq][ip];
	  q2[iq * NP + ip] = atm2->q[iq][ip];
	}
    }

    /* Init... */
    np = 0;
    ahtd = avtd = rhtd = rvtd = 0;
    for (iq = 0; iq < ctl.nq; iq++)
      aqtd[iq] = atce1[iq] = atce2[iq] = rqtd[iq] = rtce1[iq] = rtce2[iq] = 0;

    /* Loop over air parcels... */
    for (ip = 0; ip < atm1->np; ip++) {

      /* Check data... */
      if (!gsl_finite(atm1->time[ip]) || !gsl_finite(atm2->time[ip]))
	continue;

      /* Check ensemble ID... */
      if (ens >= 0 && ctl.qnt_ens >= 0 && atm1->q[ctl.qnt_ens][ip] != ens)
	continue;
      if (ens >= 0 && ctl.qnt_ens >= 0 && atm2->q[ctl.qnt_ens][ip] != ens)
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
      ahtd += DIST(x1, x2);
      avtd += fabs(z1 - z2);
      for (iq = 0; iq < ctl.nq; iq++)
	aqtd[iq] += fabs(atm1->q[iq][ip] - atm2->q[iq][ip]);

      /* Calculate relative transport deviations... */
      if (f > 3) {

	/* Get trajectory lengths... */
	geo2cart(0, lon1_old[ip], lat1_old[ip], x0);
	lh1[ip] += DIST(x0, x1);
	lv1[ip] += fabs(z1_old[ip] - z1);

	geo2cart(0, lon2_old[ip], lat2_old[ip], x0);
	lh2[ip] += DIST(x0, x2);
	lv2[ip] += fabs(z2_old[ip] - z2);

	/* Get relative transport deviations... */
	if (lh1[ip] + lh2[ip] > 0)
	  rhtd += 200. * DIST(x1, x2) / (lh1[ip] + lh2[ip]);
	if (lv1[ip] + lv2[ip] > 0)
	  rvtd += 200. * fabs(z1 - z2) / (lv1[ip] + lv2[ip]);
	for (iq = 0; iq < ctl.nq; iq++)
	  rqtd[iq] += 200. * fabs(atm1->q[iq][ip] - atm2->q[iq][ip])
	    / (fabs(atm1->q[iq][ip]) + fabs(atm2->q[iq][ip]));

	/* Get tracer conservation errors... */
	for (iq = 0; iq < ctl.nq; iq++) {
	  atce1[iq] += fabs(atm1->q[iq][ip] - q1[iq * NP + ip]);
	  rtce1[iq] += 200. * fabs(atm1->q[iq][ip] - q1[iq * NP + ip])
	    / (fabs(atm1->q[iq][ip]) + fabs(q1[iq * NP + ip]));
	  atce2[iq] += fabs(atm2->q[iq][ip] - q2[iq * NP + ip]);
	  rtce2[iq] += 200. * fabs(atm2->q[iq][ip] - q2[iq * NP + ip])
	    / (fabs(atm2->q[iq][ip]) + fabs(q2[iq * NP + ip]));
	}
      }

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

    /* Write output... */
    fprintf(out, "%.2f %.2f %g %g %g %g", t, t - t0,
	    ahtd / np, rhtd / np, avtd / np, rvtd / np);
    for (iq = 0; iq < ctl.nq; iq++) {
      fprintf(out, " ");
      fprintf(out, ctl.qnt_format[iq], aqtd[iq] / np);
      fprintf(out, " ");
      fprintf(out, ctl.qnt_format[iq], rqtd[iq] / np);
    }
    for (iq = 0; iq < ctl.nq; iq++) {
      fprintf(out, " ");
      fprintf(out, ctl.qnt_format[iq], atce1[iq] / np);
      fprintf(out, " ");
      fprintf(out, ctl.qnt_format[iq], rtce1[iq] / np);
    }
    for (iq = 0; iq < ctl.nq; iq++) {
      fprintf(out, " ");
      fprintf(out, ctl.qnt_format[iq], atce2[iq] / np);
      fprintf(out, " ");
      fprintf(out, ctl.qnt_format[iq], rtce2[iq] / np);
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

  return EXIT_SUCCESS;
}
