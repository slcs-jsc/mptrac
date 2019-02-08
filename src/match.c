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
  
  Copyright (C) 2013-2018 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Calculate deviations between two trajectories.
*/

#include "libtrac.h"

int main(
  int const argc,
  char const *argv[]) {

  ctl_t ctl;

  atm_t *atm1, *atm2, *atm3;

  FILE *out;

  char filename[LEN];

  double filter_dt, x1[3], x2[3], dh, dq[NQ], dv, lh = 0, lt = 0, lv = 0;

  int filter, ip1, ip2, iq, n;

  /* Allocate... */
  ALLOC(atm1, atm_t, 1);
  ALLOC(atm2, atm_t, 1);
  ALLOC(atm3, atm_t, 1);

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <atm_test> <atm_ref> <outfile>");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  filter = (int) scan_ctl(argv[1], argc, argv, "FILTER", -1, "0", NULL);
  filter_dt = scan_ctl(argv[1], argc, argv, "FILTER_DT", -1, "0", NULL);

  /* Read atmospheric data... */
  read_atm(argv[2], &ctl, atm1);
  read_atm(argv[3], &ctl, atm2);

  /* Write info... */
  printf("Write transport deviations: %s\n", argv[4]);

  /* Create output file... */
  if (!(out = fopen(argv[4], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = altitude [km]\n"
	  "# $3 = longitude [deg]\n" "# $4 = latitude [deg]\n");
  for (iq = 0; iq < ctl.nq; iq++)
    fprintf(out, "# $%i = %s [%s]\n", iq + 5, ctl.qnt_name[iq],
	    ctl.qnt_unit[iq]);
  fprintf(out,
	  "# $%d = trajectory time [s]\n"
	  "# $%d = vertical length of trajectory [km]\n"
	  "# $%d = horizontal length of trajectory [km]\n"
	  "# $%d = vertical deviation [km]\n"
	  "# $%d = horizontal deviation [km]\n",
	  5 + ctl.nq, 6 + ctl.nq, 7 + ctl.nq, 8 + ctl.nq, 9 + ctl.nq);
  for (iq = 0; iq < ctl.nq; iq++)
    fprintf(out, "# $%d = %s deviation [%s]\n", ctl.nq + iq + 10,
	    ctl.qnt_name[iq], ctl.qnt_unit[iq]);
  fprintf(out, "\n");

  /* Filtering of reference time series... */
  if (filter) {

    /* Copy data... */
    memcpy(atm3, atm2, sizeof(atm_t));

    /* Loop over data points... */
    for (ip1 = 0; ip1 < atm2->np; ip1++) {
      n = 0;
      atm2->p[ip1] = 0;
      for (iq = 0; iq < ctl.nq; iq++)
	atm2->q[iq][ip1] = 0;
      for (ip2 = 0; ip2 < atm2->np; ip2++)
	if (fabs(atm2->time[ip1] - atm2->time[ip2]) < filter_dt) {
	  atm2->p[ip1] += atm3->p[ip2];
	  for (iq = 0; iq < ctl.nq; iq++)
	    atm2->q[iq][ip1] += atm3->q[iq][ip2];
	  n++;
	}
      atm2->p[ip1] /= n;
      for (iq = 0; iq < ctl.nq; iq++)
	atm2->q[iq][ip1] /= n;
    }

    /* Write filtered data... */
    sprintf(filename, "%s.filt", argv[3]);
    write_atm(filename, &ctl, atm2, 0);
  }

  /* Loop over air parcels (reference data)... */
  for (ip2 = 0; ip2 < atm2->np; ip2++) {

    /* Get trajectory length... */
    if (ip2 > 0) {
      geo2cart(0, atm2->lon[ip2 - 1], atm2->lat[ip2 - 1], x1);
      geo2cart(0, atm2->lon[ip2], atm2->lat[ip2], x2);
      lh += DIST(x1, x2);
      lv += fabs(Z(atm2->p[ip2 - 1]) - Z(atm2->p[ip2]));
      lt = fabs(atm2->time[ip2] - atm2->time[0]);
    }

    /* Init... */
    n = 0;
    dh = 0;
    dv = 0;
    for (iq = 0; iq < ctl.nq; iq++)
      dq[iq] = 0;
    geo2cart(0, atm2->lon[ip2], atm2->lat[ip2], x2);

    /* Find corresponding time step (test data)... */
    for (ip1 = 0; ip1 < atm1->np; ip1++)
      if (fabs(atm1->time[ip1] - atm2->time[ip2])
	  < (filter ? filter_dt : 0.1)) {

	/* Calculate deviations... */
	geo2cart(0, atm1->lon[ip1], atm1->lat[ip1], x1);
	dh += DIST(x1, x2);
	dv += Z(atm1->p[ip1]) - Z(atm2->p[ip2]);
	for (iq = 0; iq < ctl.nq; iq++)
	  dq[iq] += atm1->q[iq][ip1] - atm2->q[iq][ip2];
	n++;
      }

    /* Write output... */
    if (n > 0) {
      fprintf(out, "%.2f %.4f %.4f %.4f",
	      atm2->time[ip2], Z(atm2->p[ip2]),
	      atm2->lon[ip2], atm2->lat[ip2]);
      for (iq = 0; iq < ctl.nq; iq++) {
	fprintf(out, " ");
	fprintf(out, ctl.qnt_format[iq], atm2->q[iq][ip2]);
      }
      fprintf(out, " %.2f %g %g %g %g", lt, lv, lh, dv / n, dh / n);
      for (iq = 0; iq < ctl.nq; iq++) {
	fprintf(out, " ");
	fprintf(out, ctl.qnt_format[iq], dq[iq] / n);
      }
      fprintf(out, "\n");
    }
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(atm1);
  free(atm2);
  free(atm3);

  return EXIT_SUCCESS;
}
