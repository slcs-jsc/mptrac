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
  
  Copyright (C) 2013-2024 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Create atmospheric data file with initial air parcel positions.
*/

#include "mptrac.h"

int main(
  int argc,
  char *argv[]) {

  atm_t *atm;

  ctl_t ctl;

  gsl_rng *rng;

  /* Allocate... */
  ALLOC(atm, atm_t, 1);

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <atm_out>");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  double t0 = scan_ctl(argv[1], argc, argv, "INIT_T0", -1, "0", NULL);
  double t1 = scan_ctl(argv[1], argc, argv, "INIT_T1", -1, "0", NULL);
  double dt = scan_ctl(argv[1], argc, argv, "INIT_DT", -1, "1", NULL);
  double z0 = scan_ctl(argv[1], argc, argv, "INIT_Z0", -1, "0", NULL);
  double z1 = scan_ctl(argv[1], argc, argv, "INIT_Z1", -1, "0", NULL);
  double dz = scan_ctl(argv[1], argc, argv, "INIT_DZ", -1, "1", NULL);
  double lon0 = scan_ctl(argv[1], argc, argv, "INIT_LON0", -1, "0", NULL);
  double lon1 = scan_ctl(argv[1], argc, argv, "INIT_LON1", -1, "0", NULL);
  double dlon = scan_ctl(argv[1], argc, argv, "INIT_DLON", -1, "1", NULL);
  double lat0 = scan_ctl(argv[1], argc, argv, "INIT_LAT0", -1, "0", NULL);
  double lat1 = scan_ctl(argv[1], argc, argv, "INIT_LAT1", -1, "0", NULL);
  double dlat = scan_ctl(argv[1], argc, argv, "INIT_DLAT", -1, "1", NULL);
  double st = scan_ctl(argv[1], argc, argv, "INIT_ST", -1, "0", NULL);
  double sz = scan_ctl(argv[1], argc, argv, "INIT_SZ", -1, "0", NULL);
  double slon = scan_ctl(argv[1], argc, argv, "INIT_SLON", -1, "0", NULL);
  double slat = scan_ctl(argv[1], argc, argv, "INIT_SLAT", -1, "0", NULL);
  double sx = scan_ctl(argv[1], argc, argv, "INIT_SX", -1, "0", NULL);
  double ut = scan_ctl(argv[1], argc, argv, "INIT_UT", -1, "0", NULL);
  double uz = scan_ctl(argv[1], argc, argv, "INIT_UZ", -1, "0", NULL);
  double ulon = scan_ctl(argv[1], argc, argv, "INIT_ULON", -1, "0", NULL);
  double ulat = scan_ctl(argv[1], argc, argv, "INIT_ULAT", -1, "0", NULL);
  int even =
    (int) scan_ctl(argv[1], argc, argv, "INIT_EVENLY", -1, "0", NULL);
  int rep = (int) scan_ctl(argv[1], argc, argv, "INIT_REP", -1, "1", NULL);
  double m = scan_ctl(argv[1], argc, argv, "INIT_MASS", -1, "0", NULL);
  double vmr = scan_ctl(argv[1], argc, argv, "INIT_VMR", -1, "0", NULL);
  double bellrad =
    scan_ctl(argv[1], argc, argv, "INIT_BELLRAD", -1, "0", NULL);

  /* Initialize random number generator... */
  gsl_rng_env_setup();
  rng = gsl_rng_alloc(gsl_rng_default);

  /* Create grid... */
  for (double t = t0; t <= t1; t += dt)
    for (double z = z0; z <= z1; z += dz)
      for (double lon = lon0; lon <= lon1; lon += dlon)
	for (double lat = lat0; lat <= lat1; lat += dlat)
	  for (int irep = 0; irep < rep; irep++) {

	    /* Set position... */
	    double rg = gsl_ran_gaussian_ziggurat(rng, st / 2.3548);
	    double ru = ut * (gsl_rng_uniform(rng) - 0.5);
	    atm->time[atm->np] = (t + rg + ru);

	    rg = gsl_ran_gaussian_ziggurat(rng, sz / 2.3548);
	    ru = uz * (gsl_rng_uniform(rng) - 0.5);
	    atm->p[atm->np] = P(z + rg + ru);

	    rg = gsl_ran_gaussian_ziggurat(rng, slon / 2.3548);
	    double rx =
	      gsl_ran_gaussian_ziggurat(rng, DX2DEG(sx, lat) / 2.3548);
	    ru = ulon * (gsl_rng_uniform(rng) - 0.5);
	    atm->lon[atm->np] = (lon + rg + rx + ru);

	    do {
	      rg = gsl_ran_gaussian_ziggurat(rng, slat / 2.3548);
	      rx = gsl_ran_gaussian_ziggurat(rng, DY2DEG(sx) / 2.3548);
	      ru = ulat * (gsl_rng_uniform(rng) - 0.5);
	      atm->lat[atm->np] = (lat + rg + rx + ru);
	    } while (even && gsl_rng_uniform(rng) >
		     fabs(cos(atm->lat[atm->np] * M_PI / 180.)));

	    /* Apply cosine bell (Williamson et al., 1992)... */
	    if (bellrad > 0) {
	      double x0[3], x1[3];
	      geo2cart(0.0, 0.5 * (lon0 + lon1), 0.5 * (lat0 + lat1), x0);
	      geo2cart(0.0, atm->lon[atm->np], atm->lat[atm->np], x1);
	      double rad = RE * acos(DOTP(x0, x1) / NORM(x0) / NORM(x1));
	      if (rad > bellrad)
		continue;
	      if (ctl.qnt_m >= 0)
		atm->q[ctl.qnt_m][atm->np] =
		  0.5 * (1. + cos(M_PI * rad / bellrad));
	      if (ctl.qnt_vmr >= 0)
		atm->q[ctl.qnt_vmr][atm->np] =
		  0.5 * (1. + cos(M_PI * rad / bellrad));
	    }

	    /* Set particle counter... */
	    if ((++atm->np) > NP)
	      ERRMSG("Too many particles!");
	  }

  /* Check number of air parcels... */
  if (atm->np <= 0)
    ERRMSG("Did not create any air parcels!");

  /* Initialize mass... */
  if (ctl.qnt_m >= 0 && bellrad <= 0)
    for (int ip = 0; ip < atm->np; ip++)
      atm->q[ctl.qnt_m][ip] = m / atm->np;

  /* Initialize volume mixing ratio... */
  if (ctl.qnt_vmr >= 0 && bellrad <= 0)
    for (int ip = 0; ip < atm->np; ip++)
      atm->q[ctl.qnt_vmr][ip] = vmr;

  /* Initialize air parcel index... */
  if (ctl.qnt_idx >= 0)
    for (int ip = 0; ip < atm->np; ip++)
      atm->q[ctl.qnt_idx][ip] = ip;

  /* Initialize age of air... */
  if (ctl.qnt_aoa >= 0)
    for (int ip = 0; ip < atm->np; ip++)
      atm->q[ctl.qnt_aoa][ip] = atm->time[ip];

  /* Save data... */
  write_atm(argv[2], &ctl, atm, 0);

  /* Free... */
  gsl_rng_free(rng);
  free(atm);

  return EXIT_SUCCESS;
}
