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
  Calculate accuracy of trajectories based on reference calculations.
*/

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {

  atm_t *atm1, *atm2;

  char *name;
  char *year, *mon, *day, *hour, *min;

  ctl_t ctl;

  FILE *out;

  double x0[3], x1[3], *lon1, *lat1, *p1, *dh1, *dv1, *lon2, *lat2, *p2,
    *dh2, *dv2, ahtd, avtd, rhtd, rvtd, t;

  int ip, f;

  /* Allocate... */
  ALLOC(atm1, atm_t, 1);
  ALLOC(atm2, atm_t, 1);
  ALLOC(lon1, double,
	NP);
  ALLOC(lat1, double,
	NP);
  ALLOC(p1, double,
	NP);
  ALLOC(dh1, double,
	NP);
  ALLOC(dv1, double,
	NP);
  ALLOC(lon2, double,
	NP);
  ALLOC(lat2, double,
	NP);
  ALLOC(p2, double,
	NP);
  ALLOC(dh2, double,
	NP);
  ALLOC(dv2, double,
	NP);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG
      ("Give parameters: <outfile> <atm1a> <atm1b> [<atm2a> <atm2b> ...]");

  if (argv[2][strlen(argv[2]) - 1] == 'c') {
    ctl.atm_iformat = 1;
  }

  /* Write info... */
  printf("Write trajectory analysis data: %s\n", argv[1]);

  /* Create output file... */
  if (!(out = fopen(argv[1], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = AHTD(t) [km]\n"
	  "# $3 = RHTD(t) [%%]\n"
	  "# $4 = AVTD(t) [km]\n" "# $5 = RVTD(t) [%%]\n\n");

  /* Loop over file pairs... */
  for (f = 2; f < argc; f += 2) {

    /* Read atmopheric data... */
    read_atm(NULL, argv[f], atm1, &ctl);
    read_atm(NULL, argv[f + 1], atm2, &ctl);

    /* Check if structs match... */
    if (atm1->np != atm2->np)
      ERRMSG("Different numbers of parcels!");
    for (ip = 0; ip < atm1->np; ip++)
      if (fabs(atm1->time[ip] - atm2->time[ip]) > 1)
	printf("WARNING! Times do not match! dt=%f\n",
	       fabs(atm1->time[ip] - atm2->time[ip]));

    /* Init... */
    ahtd = avtd = rhtd = rvtd = 0;

    /* Loop over air parcels... */
    for (ip = 0; ip < atm1->np; ip++) {

      /* Calculate total length of trajectories... */
      if (f > 2) {
	dv1[ip] += fabs(Z(p1[ip]) - Z(atm1->p[ip]));
	geo2cart(0, atm1->lon[ip], atm1->lat[ip], x0);
	geo2cart(0, lon1[ip], lat1[ip], x1);
	dh1[ip] += DIST(x0, x1);

	dv2[ip] += fabs(Z(p2[ip]) - Z(atm2->p[ip]));
	geo2cart(0, atm2->lon[ip], atm2->lat[ip], x0);
	geo2cart(0, lon2[ip], lat2[ip], x1);
	dh2[ip] += DIST(x0, x1);
      }

      /* Save last locations of air parcels... */
      lon1[ip] = atm1->lon[ip];
      lat1[ip] = atm1->lat[ip];
      p1[ip] = atm1->p[ip];
      lon2[ip] = atm2->lon[ip];
      lat2[ip] = atm2->lat[ip];
      p2[ip] = atm2->p[ip];

      /* Sum up the vertical error... */
      avtd += fabs(Z(atm1->p[ip]) - Z(atm2->p[ip])) / atm1->np;

      /* Sum up the horizontal error... */
      geo2cart(0, atm1->lon[ip], atm1->lat[ip], x0);
      geo2cart(0, atm2->lon[ip], atm2->lat[ip], x1);
      ahtd += DIST(x0, x1) / atm1->np;

      /* Sum up the relative transport devation... */
      if (f > 2) {
	t = ((200. * DIST(x0, x1)) / (dh1[ip] + dh2[ip]) / atm1->np);
	if (gsl_finite(t))
	  rhtd += t;
	t = ((200. * fabs(Z(atm1->p[ip]) - Z(atm2->p[ip]))) /
	     (dv1[ip] + dv2[ip]) / atm1->np);
	if (gsl_finite(t))
	  rvtd += t;
      }
    }

    /* Get date from filename... */
    for (ip = (int) strlen(argv[f]) - 1; argv[f][ip] != '/' || ip == 0; ip--);
    name = strtok(&(argv[f][ip]), "_");
    year = strtok(NULL, "_");
    mon = strtok(NULL, "_");
    day = strtok(NULL, "_");
    hour = strtok(NULL, "_");
    name = strtok(NULL, "_");
    min = strtok(name, ".");
    time2jsec(atoi(year), atoi(mon), atoi(day), atoi(hour), atoi(min), 0, 0,
	      &t);

    /* Write output... */
    fprintf(out, "%.2f %g %g %g %g\n", t, ahtd, rhtd, avtd, rvtd);
    printf("%.2f %g %g %g %g\n", t, ahtd, rhtd, avtd, rvtd);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(atm1);
  free(atm2);
  free(lon1);
  free(lat1);
  free(p1);
  free(dh1);
  free(dv1);
  free(lon2);
  free(lat2);
  free(p2);
  free(dh2);
  free(dv2);

  return EXIT_SUCCESS;
}
