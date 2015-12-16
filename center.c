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
  Calculate center of mass of air parcels.
*/

#include "libtrac.h"
#include <gsl/gsl_sort.h>

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  atm_t *atm;

  FILE *out;

  char *name, *year, *mon, *day, *hour, *min;

  double latm, lats, lonm, lons, t, zm, zs;

  int i, f, ip;

  /* Allocate... */
  ALLOC(atm, atm_t, 1);

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <outfile> <atm1> [<atm2> ...]");

  /* Write info... */
  printf("Write center of mass data: %s\n", argv[1]);

  /* Create output file... */
  if (!(out = fopen(argv[1], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1  = time [s]\n"
	  "# $2  = altitude (mean) [km]\n"
	  "# $3  = altitude (sigma) [km]\n"
	  "# $4  = altitude (minimum) [km]\n"
	  "# $5  = altitude (10%% percentile) [km]\n"
	  "# $6  = altitude (1st quarter) [km]\n"
	  "# $7  = altitude (median) [km]\n"
	  "# $8  = altitude (3rd quarter) [km]\n"
	  "# $9  = altitude (90%% percentile) [km]\n"
	  "# $10 = altitude (maximum) [km]\n");
  fprintf(out,
	  "# $11 = longitude (mean) [deg]\n"
	  "# $12 = longitude (sigma) [deg]\n"
	  "# $13 = longitude (minimum) [deg]\n"
	  "# $14 = longitude (10%% percentile) [deg]\n"
	  "# $15 = longitude (1st quarter) [deg]\n"
	  "# $16 = longitude (median) [deg]\n"
	  "# $17 = longitude (3rd quarter) [deg]\n"
	  "# $18 = longitude (90%% percentile) [deg]\n"
	  "# $19 = longitude (maximum) [deg]\n");
  fprintf(out,
	  "# $20 = latitude (mean) [deg]\n"
	  "# $21 = latitude (sigma) [deg]\n"
	  "# $22 = latitude (minimum) [deg]\n"
	  "# $23 = latitude (10%% percentile) [deg]\n"
	  "# $24 = latitude (1st quarter) [deg]\n"
	  "# $25 = latitude (median) [deg]\n"
	  "# $26 = latitude (3rd quarter) [deg]\n"
	  "# $27 = latitude (90%% percentile) [deg]\n"
	  "# $28 = latitude (maximum) [deg]\n\n");

  /* Loop over files... */
  for (f = 2; f < argc; f++) {

    /* Read atmopheric data... */
    read_atm(argv[f], atm, &ctl);

    /* Initialize... */
    zm = zs = 0;
    lonm = lons = 0;
    latm = lats = 0;

    /* Calculate mean and standard deviation... */
    for (ip = 0; ip < atm->np; ip++) {
      zm += Z(atm->p[ip]) / atm->np;
      lonm += atm->lon[ip] / atm->np;
      latm += atm->lat[ip] / atm->np;
      zs += gsl_pow_2(Z(atm->p[ip])) / atm->np;
      lons += gsl_pow_2(atm->lon[ip]) / atm->np;
      lats += gsl_pow_2(atm->lat[ip]) / atm->np;
    }

    /* Normalize... */
    zs = sqrt(zs - gsl_pow_2(zm));
    lons = sqrt(lons - gsl_pow_2(lonm));
    lats = sqrt(lats - gsl_pow_2(latm));

    /* Sort arrays... */
    gsl_sort(atm->p, 1, (size_t) atm->np);
    gsl_sort(atm->lon, 1, (size_t) atm->np);
    gsl_sort(atm->lat, 1, (size_t) atm->np);

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

    /* Write data... */
    fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g "
	    "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	    t, zm, zs, Z(atm->p[atm->np - 1]),
	    Z(atm->p[atm->np - atm->np / 10]),
	    Z(atm->p[atm->np - atm->np / 4]),
	    Z(atm->p[atm->np / 2]), Z(atm->p[atm->np / 4]),
	    Z(atm->p[atm->np / 10]), Z(atm->p[0]),
	    lonm, lons, atm->lon[0], atm->lon[atm->np / 10],
	    atm->lon[atm->np / 4], atm->lon[atm->np / 2],
	    atm->lon[atm->np - atm->np / 4],
	    atm->lon[atm->np - atm->np / 10],
	    atm->lon[atm->np - 1],
	    latm, lats, atm->lat[0], atm->lat[atm->np / 10],
	    atm->lat[atm->np / 4], atm->lat[atm->np / 2],
	    atm->lat[atm->np - atm->np / 4],
	    atm->lat[atm->np - atm->np / 10], atm->lat[atm->np - 1]);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(atm);

  return EXIT_SUCCESS;
}
