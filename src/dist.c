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

  double x0[3], x1[3], x2[3], *lon1, *lat1, *p1, *lh1, *lv1,
    *lon2, *lat2, *p2, *lh2, *lv2, ahtd, avtd, aqtd[NQ], rhtd, rvtd, t;

  int ens, f, ip, iq, np, year, mon, day, hour, min;

  /* Allocate... */
  ALLOC(atm1, atm_t, 1);
  ALLOC(atm2, atm_t, 1);
  ALLOC(lon1, double,
	NP);
  ALLOC(lat1, double,
	NP);
  ALLOC(p1, double,
	NP);
  ALLOC(lh1, double,
	NP);
  ALLOC(lv1, double,
	NP);
  ALLOC(lon2, double,
	NP);
  ALLOC(lat2, double,
	NP);
  ALLOC(p2, double,
	NP);
  ALLOC(lh2, double,
	NP);
  ALLOC(lv2, double,
	NP);

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <outfile> <atm1a> <atm1b>"
	   " [<atm2a> <atm2b> ...]");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  ens = (int) scan_ctl(argv[1], argc, argv, "DIST_ENS", -1, "-1", NULL);

  /* Write info... */
  printf("Write transport deviations: %s\n", argv[2]);

  /* Create output file... */
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = AHTD [km]\n"
	  "# $3 = RHTD [km]\n" "# $4 = AVTD [km]\n" "# $5 = RVTD [km]\n");
  for (iq = 0; iq < ctl.nq; iq++)
    fprintf(out,
	    "# $%d = AQTD (%s) [%s]\n",
	    6 + iq, ctl.qnt_name[iq], ctl.qnt_unit[iq]);
  fprintf(out, "\n");

  /* Loop over file pairs... */
  for (f = 3; f < argc; f += 2) {

    /* Read atmopheric data... */
    read_atm(argv[f], &ctl, atm1);
    read_atm(argv[f + 1], &ctl, atm2);

    /* Check if structs match... */
    if (atm1->np != atm2->np)
      ERRMSG("Different numbers of parcels!");
    for (ip = 0; ip < atm1->np; ip++)
      if (atm1->time[ip] != atm2->time[ip])
	ERRMSG("Times do not match!");

    /* Init... */
    np = 0;
    ahtd = avtd = rhtd = rvtd = 0;
    for (iq = 0; iq < ctl.nq; iq++)
      aqtd[iq] = 0;

    /* Loop over air parcels... */
    for (ip = 0; ip < atm1->np; ip++)
      if (ens < 0 || (ctl.qnt_ens >= 0 && atm1->q[ctl.qnt_ens][ip] == ens)) {

	/* Get Cartesian coordinates... */
	geo2cart(0, atm1->lon[ip], atm1->lat[ip], x1);
	geo2cart(0, atm2->lon[ip], atm2->lat[ip], x2);

	/* Calculate absolute transport deviations... */
	ahtd += DIST(x1, x2);
	avtd += fabs(Z(atm1->p[ip]) - Z(atm2->p[ip]));
	for (iq = 0; iq < ctl.nq; iq++)
	  aqtd[iq] += fabs(atm1->q[iq][ip] - atm2->q[iq][ip]);

	/* Calculate relative transport deviations... */
	if (f > 3) {

	  /* Get trajectory lengths... */
	  geo2cart(0, lon1[ip], lat1[ip], x0);
	  lh1[ip] += DIST(x0, x1);
	  lv1[ip] += fabs(Z(p1[ip]) - Z(atm1->p[ip]));

	  geo2cart(0, lon2[ip], lat2[ip], x0);
	  lh2[ip] += DIST(x0, x2);
	  lv2[ip] += fabs(Z(p2[ip]) - Z(atm2->p[ip]));

	  /* Get relative transport deviations... */
	  if (lh1[ip] + lh2[ip] > 0)
	    rhtd += 200. * DIST(x1, x2) / (lh1[ip] + lh2[ip]);
	  if (lv1[ip] + lv2[ip] > 0)
	    rvtd += 200. * fabs(Z(atm1->p[ip]) - Z(atm2->p[ip]))
	      / (lv1[ip] + lv2[ip]);
	}

	/* Save positions of air parcels... */
	lon1[ip] = atm1->lon[ip];
	lat1[ip] = atm1->lat[ip];
	p1[ip] = atm1->p[ip];

	lon2[ip] = atm2->lon[ip];
	lat2[ip] = atm2->lat[ip];
	p2[ip] = atm2->p[ip];

	/* Increment air parcel counter... */
	np++;
      }

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

    /* Write output... */
    fprintf(out, "%.2f %g %g %g %g", t,
	    ahtd / np, rhtd / np, avtd / np, rvtd / np);
    for (iq = 0; iq < ctl.nq; iq++) {
      fprintf(out, " ");
      fprintf(out, ctl.qnt_format[iq], aqtd[iq] / np);
    }
    fprintf(out, "\n");
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(atm1);
  free(atm2);
  free(lon1);
  free(lat1);
  free(p1);
  free(lh1);
  free(lv1);
  free(lon2);
  free(lat2);
  free(p2);
  free(lh2);
  free(lv2);

  return EXIT_SUCCESS;
}
