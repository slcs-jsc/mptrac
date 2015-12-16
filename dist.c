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
  
  Copright (C) 2013-2015 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Calculate transport deviations of trajectories.
*/

#include "libtrac.h"
#include <gsl/gsl_sort.h>

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  atm_t *atm1, *atm2;

  FILE *out;

  char *name, *year, *mon, *day, *hour, *min;

  double aux, x0[3], x1[3], x2[3], *lon1, *lat1, *p1, *lh1, *lv1,
    *lon2, *lat2, *p2, *lh2, *lv2, ahtd, avtd, ahtd2, avtd2,
    rhtd, rvtd, rhtd2, rvtd2, t, *dh, *dv;

  int f, i, ip, iph, ipv;

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
  ALLOC(dh, double,
	NP);
  ALLOC(dv, double,
	NP);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG
      ("Give parameters: <outfile> <atm1a> <atm1b> [<atm2a> <atm2b> ...]");

  /* Write info... */
  printf("Write transport deviations: %s\n", argv[1]);

  /* Create output file... */
  if (!(out = fopen(argv[1], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1  = time [s]\n"
	  "# $2  = AHTD (mean) [km]\n"
	  "# $3  = AHTD (sigma) [km]\n"
	  "# $4  = AHTD (minimum) [km]\n"
	  "# $5  = AHTD (10%% percentile) [km]\n"
	  "# $6  = AHTD (1st quartile) [km]\n"
	  "# $7  = AHTD (median) [km]\n"
	  "# $8  = AHTD (3rd quartile) [km]\n"
	  "# $9  = AHTD (90%% percentile) [km]\n"
	  "# $10 = AHTD (maximum) [km]\n"
	  "# $11 = AHTD (maximum trajectory index)\n"
	  "# $12 = RHTD (mean) [%%]\n" "# $13 = RHTD (sigma) [%%]\n");
  fprintf(out,
	  "# $14 = AVTD (mean) [km]\n"
	  "# $15 = AVTD (sigma) [km]\n"
	  "# $16 = AVTD (minimum) [km]\n"
	  "# $17 = AVTD (10%% percentile) [km]\n"
	  "# $18 = AVTD (1st quartile) [km]\n"
	  "# $19 = AVTD (median) [km]\n"
	  "# $20 = AVTD (3rd quartile) [km]\n"
	  "# $21 = AVTD (90%% percentile) [km]\n"
	  "# $22 = AVTD (maximum) [km]\n"
	  "# $23 = AVTD (maximum trajectory index)\n"
	  "# $24 = RVTD (mean) [%%]\n" "# $25 = RVTD (sigma) [%%]\n\n");

  /* Loop over file pairs... */
  for (f = 2; f < argc; f += 2) {

    /* Read atmopheric data... */
    read_atm(argv[f], atm1, &ctl);
    read_atm(argv[f + 1], atm2, &ctl);

    /* Check if structs match... */
    if (atm1->np != atm2->np)
      ERRMSG("Different numbers of parcels!");
    for (ip = 0; ip < atm1->np; ip++)
      if (atm1->time[ip] != atm2->time[ip])
	ERRMSG("Times do not match!");

    /* Init... */
    ahtd = ahtd2 = 0;
    avtd = avtd2 = 0;
    rhtd = rhtd2 = 0;
    rvtd = rvtd2 = 0;

    /* Loop over air parcels... */
    for (ip = 0; ip < atm1->np; ip++) {

      /* Get Cartesian coordinates... */
      geo2cart(0, atm1->lon[ip], atm1->lat[ip], x1);
      geo2cart(0, atm2->lon[ip], atm2->lat[ip], x2);

      /* Calculate absolute transport deviations... */
      dh[ip] = DIST(x1, x2);
      ahtd += dh[ip];
      ahtd2 += gsl_pow_2(dh[ip]);

      dv[ip] = fabs(Z(atm1->p[ip]) - Z(atm2->p[ip]));
      avtd += dv[ip];
      avtd2 += gsl_pow_2(dv[ip]);

      /* Calculate relative transport deviations... */
      if (f > 2) {

	/* Get trajectory lengths... */
	geo2cart(0, lon1[ip], lat1[ip], x0);
	lh1[ip] += DIST(x0, x1);
	lv1[ip] += fabs(Z(p1[ip]) - Z(atm1->p[ip]));

	geo2cart(0, lon2[ip], lat2[ip], x0);
	lh2[ip] += DIST(x0, x2);
	lv2[ip] += fabs(Z(p2[ip]) - Z(atm2->p[ip]));

	/* Get relative transport devations... */
	if (lh1[ip] + lh2[ip] > 0) {
	  aux = 200. * DIST(x1, x2) / (lh1[ip] + lh2[ip]);
	  rhtd += aux;
	  rhtd2 += gsl_pow_2(aux);
	}
	if (lv1[ip] + lv2[ip] > 0) {
	  aux =
	    200. * fabs(Z(atm1->p[ip]) - Z(atm2->p[ip])) / (lv1[ip] +
							    lv2[ip]);
	  rvtd += aux;
	  rvtd2 += gsl_pow_2(aux);
	}
      }

      /* Save positions of air parcels... */
      lon1[ip] = atm1->lon[ip];
      lat1[ip] = atm1->lat[ip];
      p1[ip] = atm1->p[ip];

      lon2[ip] = atm2->lon[ip];
      lat2[ip] = atm2->lat[ip];
      p2[ip] = atm2->p[ip];
    }

    /* Get indices of trajectories with maximum errors... */
    iph = (int) gsl_stats_max_index(dh, 1, (size_t) atm1->np);
    ipv = (int) gsl_stats_max_index(dv, 1, (size_t) atm1->np);

    /* Sort distances to calculate percentiles... */
    gsl_sort(dh, 1, (size_t) atm1->np);
    gsl_sort(dv, 1, (size_t) atm1->np);

    /* Get date from filename... */
    for (i = (int) strlen(argv[f]) - 1; argv[f][i] != '/' || i == 0; i--);
    name = strtok(&(argv[f][i]), "_");
    year = strtok(NULL, "_");
    mon = strtok(NULL, "_");
    day = strtok(NULL, "_");
    hour = strtok(NULL, "_");
    name = strtok(NULL, "_");	/* TODO: Why another "name" here? */
    min = strtok(name, ".");
    time2jsec(atoi(year), atoi(mon), atoi(day), atoi(hour), atoi(min), 0, 0,
	      &t);

    /* Write output... */
    fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %d %g %g"
	    " %g %g %g %g %g %g %g %g %g %d %g %g\n", t,
	    ahtd / atm1->np,
	    sqrt(ahtd2 / atm1->np - gsl_pow_2(ahtd / atm1->np)),
	    dh[0], dh[atm1->np / 10], dh[atm1->np / 4], dh[atm1->np / 2],
	    dh[atm1->np - atm1->np / 4], dh[atm1->np - atm1->np / 10],
	    dh[atm1->np - 1], iph, rhtd / atm1->np,
	    sqrt(rhtd2 / atm1->np - gsl_pow_2(rhtd / atm1->np)),
	    avtd / atm1->np,
	    sqrt(avtd2 / atm1->np - gsl_pow_2(avtd / atm1->np)),
	    dv[0], dv[atm1->np / 10], dv[atm1->np / 4], dv[atm1->np / 2],
	    dv[atm1->np - atm1->np / 4], dv[atm1->np - atm1->np / 10],
	    dv[atm1->np - 1], ipv, rvtd / atm1->np,
	    sqrt(rvtd2 / atm1->np - gsl_pow_2(rvtd / atm1->np)));
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
  free(dh);
  free(dv);

  return EXIT_SUCCESS;
}
