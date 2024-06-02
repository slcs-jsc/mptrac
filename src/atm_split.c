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
  Split air parcels into a larger set of parcels.
*/

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {

  atm_t *atm, *atm2;

  ctl_t ctl;

  gsl_rng *rng;

  char kernel[LEN];

  double k, kw[EP], kz[EP], mmax = 0, mtot = 0, z, zmin = 0, zmax = 0;

  int ip, nk = 0;

  /* Allocate... */
  ALLOC(atm, atm_t, 1);
  ALLOC(atm2, atm_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <atm_in> <atm_out>");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  int n = (int) scan_ctl(argv[1], argc, argv, "SPLIT_N", -1, "", NULL);
  double m = scan_ctl(argv[1], argc, argv, "SPLIT_M", -1, "-999", NULL);
  double um = scan_ctl(argv[1], argc, argv, "SPLIT_UM", -1, "0", NULL);
  double dt = scan_ctl(argv[1], argc, argv, "SPLIT_DT", -1, "0", NULL);
  double t0 = scan_ctl(argv[1], argc, argv, "SPLIT_T0", -1, "0", NULL);
  double t1 = scan_ctl(argv[1], argc, argv, "SPLIT_T1", -1, "0", NULL);
  double dz = scan_ctl(argv[1], argc, argv, "SPLIT_DZ", -1, "0", NULL);
  double z0 = scan_ctl(argv[1], argc, argv, "SPLIT_Z0", -1, "0", NULL);
  double z1 = scan_ctl(argv[1], argc, argv, "SPLIT_Z1", -1, "0", NULL);
  double dx = scan_ctl(argv[1], argc, argv, "SPLIT_DX", -1, "0", NULL);
  double lon0 = scan_ctl(argv[1], argc, argv, "SPLIT_LON0", -1, "0", NULL);
  double lon1 = scan_ctl(argv[1], argc, argv, "SPLIT_LON1", -1, "0", NULL);
  double lat0 = scan_ctl(argv[1], argc, argv, "SPLIT_LAT0", -1, "0", NULL);
  double lat1 = scan_ctl(argv[1], argc, argv, "SPLIT_LAT1", -1, "0", NULL);
  scan_ctl(argv[1], argc, argv, "SPLIT_KERNEL", -1, "-", kernel);

  /* Init random number generator... */
  gsl_rng_env_setup();
  rng = gsl_rng_alloc(gsl_rng_default);

  /* Read atmospheric data... */
  if (!read_atm(argv[2], &ctl, atm))
    ERRMSG("Cannot open file!");

  /* Read kernel function... */
  if (kernel[0] != '-') {
    read_kernel(kernel, kz, kw, &nk);
    zmax = gsl_stats_max(kz, 1, (size_t) nk);
    zmin = gsl_stats_min(kz, 1, (size_t) nk);
  }

  /* Get total and maximum mass... */
  if (ctl.qnt_m >= 0)
    for (ip = 0; ip < atm->np; ip++) {
      mtot += atm->q[ctl.qnt_m][ip];
      mmax = MAX(mmax, atm->q[ctl.qnt_m][ip]);
    }
  if (m >= 0)
    mtot = m;

  /* Loop over air parcels... */
  for (int i = 0; i < n; i++) {

    /* Select air parcel... */
    if (ctl.qnt_m >= 0)
      do {
	ip = (int) gsl_rng_uniform_int(rng, (long unsigned int) atm->np);
      } while (gsl_rng_uniform(rng) > atm->q[ctl.qnt_m][ip] / mmax);
    else
      ip = (int) gsl_rng_uniform_int(rng, (long unsigned int) atm->np);

    /* Set time... */
    if (t1 > t0)
      atm2->time[atm2->np] = t0 + (t1 - t0) * gsl_rng_uniform_pos(rng);
    else
      atm2->time[atm2->np] = atm->time[ip]
	+ gsl_ran_gaussian_ziggurat(rng, dt / 2.3548);

    /* Set vertical position... */
    do {
      if (nk > 0) {
	do {
	  z = zmin + (zmax - zmin) * gsl_rng_uniform_pos(rng);
	  k = kernel_weight(kz, kw, nk, P(z));
	} while (gsl_rng_uniform(rng) > k);
	atm2->p[atm2->np] = P(z);
      } else if (z1 > z0)
	atm2->p[atm2->np] = P(z0 + (z1 - z0) * gsl_rng_uniform_pos(rng));
      else
	atm2->p[atm2->np] = atm->p[ip]
	  + DZ2DP(gsl_ran_gaussian_ziggurat(rng, dz / 2.3548), atm->p[ip]);
    } while (atm2->p[atm2->np] < P(100.) || atm2->p[atm2->np] > P(-1.));

    /* Set horizontal position... */
    if (lon1 > lon0 && lat1 > lat0) {
      atm2->lon[atm2->np] = lon0 + (lon1 - lon0) * gsl_rng_uniform_pos(rng);
      atm2->lat[atm2->np] = lat0 + (lat1 - lat0) * gsl_rng_uniform_pos(rng);
    } else {
      atm2->lon[atm2->np] = atm->lon[ip]
	+ gsl_ran_gaussian_ziggurat(rng, DX2DEG(dx, atm->lat[ip]) / 2.3548);
      atm2->lat[atm2->np] = atm->lat[ip]
	+ gsl_ran_gaussian_ziggurat(rng, DY2DEG(dx) / 2.3548);
    }

    /* Copy quantities... */
    for (int iq = 0; iq < ctl.nq; iq++)
      atm2->q[iq][atm2->np] = atm->q[iq][ip];

    /* Adjust mass... */
    if (ctl.qnt_m >= 0)
      atm2->q[ctl.qnt_m][atm2->np]
	= (mtot + (um > 0 ? um * (gsl_rng_uniform_pos(rng) - 0.5) : 0.0)) / n;

    /* Adjust air parcel index... */
    if (ctl.qnt_idx >= 0)
      atm2->q[ctl.qnt_idx][atm2->np] = atm2->np;

    /* Increment particle counter... */
    if ((++atm2->np) > NP)
      ERRMSG("Too many air parcels!");
  }

  /* Save data and close file... */
  write_atm(argv[3], &ctl, atm2, 0);

  /* Free... */
  free(atm);
  free(atm2);

  return EXIT_SUCCESS;
}
