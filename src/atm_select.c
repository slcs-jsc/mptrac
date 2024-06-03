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
  Extract subsets of air parcels from atmospheric data files.
*/

#include "mptrac.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  atm_t *atm, *atm2;

  /* Allocate... */
  ALLOC(atm, atm_t, 1);
  ALLOC(atm2, atm_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <atm_select> <atm1> [<atm2> ...]");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  int stride =
    (int) scan_ctl(argv[1], argc, argv, "SELECT_STRIDE", -1, "1", NULL);
  int idx0 =
    (int) scan_ctl(argv[1], argc, argv, "SELECT_IDX0", -1, "-999", NULL);
  int idx1 =
    (int) scan_ctl(argv[1], argc, argv, "SELECT_IDX1", -1, "-999", NULL);
  int ip0 =
    (int) scan_ctl(argv[1], argc, argv, "SELECT_IP0", -1, "-999", NULL);
  int ip1 =
    (int) scan_ctl(argv[1], argc, argv, "SELECT_IP1", -1, "-999", NULL);
  double t0 = scan_ctl(argv[1], argc, argv, "SELECT_T0", -1, "0", NULL);
  double t1 = scan_ctl(argv[1], argc, argv, "SELECT_T1", -1, "0", NULL);
  double p0 = P(scan_ctl(argv[1], argc, argv, "SELECT_Z0", -1, "0", NULL));
  double p1 = P(scan_ctl(argv[1], argc, argv, "SELECT_Z1", -1, "0", NULL));
  double theta0 =
    scan_ctl(argv[1], argc, argv, "SELECT_THETA0", -1, "0", NULL);
  double theta1 =
    scan_ctl(argv[1], argc, argv, "SELECT_THETA1", -1, "0", NULL);
  double lon0 = scan_ctl(argv[1], argc, argv, "SELECT_LON0", -1, "0", NULL);
  double lon1 = scan_ctl(argv[1], argc, argv, "SELECT_LON1", -1, "0", NULL);
  double lat0 = scan_ctl(argv[1], argc, argv, "SELECT_LAT0", -1, "0", NULL);
  double lat1 = scan_ctl(argv[1], argc, argv, "SELECT_LAT1", -1, "0", NULL);
  double r0 = scan_ctl(argv[1], argc, argv, "SELECT_R0", -1, "0", NULL);
  double r1 = scan_ctl(argv[1], argc, argv, "SELECT_R1", -1, "0", NULL);
  double rlon = scan_ctl(argv[1], argc, argv, "SELECT_RLON", -1, "0", NULL);
  double rlat = scan_ctl(argv[1], argc, argv, "SELECT_RLAT", -1, "0", NULL);

  /* Get Cartesian coordinates... */
  double x0[3], x1[3];
  geo2cart(0, rlon, rlat, x0);

  /* Loop over files... */
  for (int f = 3; f < argc; f++) {

    /* Read atmopheric data... */
    if (!read_atm(argv[f], &ctl, atm))
      continue;

    /* Adjust range of air parcels... */
    if (ip0 < 0)
      ip0 = 0;
    ip0 = MIN(ip0, atm->np - 1);
    if (ip1 < 0)
      ip1 = atm->np - 1;
    ip1 = MIN(ip1, atm->np - 1);
    if (ip1 < ip0)
      ip1 = ip0;

    /* Loop over air parcels... */
    for (int ip = ip0; ip <= ip1; ip += stride) {

      /* Check air parcel index... */
      if (ctl.qnt_idx >= 0 && idx0 >= 0 && idx1 >= 0)
	if (atm->q[ctl.qnt_idx][ip] < idx0 || atm->q[ctl.qnt_idx][ip] > idx1)
	  continue;

      /* Check time... */
      if (t0 != t1)
	if ((t1 > t0 && (atm->time[ip] < t0 || atm->time[ip] > t1))
	    || (t1 < t0 && (atm->time[ip] < t0 && atm->time[ip] > t1)))
	  continue;

      /* Check vertical distance... */
      if (p0 != p1)
	if ((p0 > p1 && (atm->p[ip] > p0 || atm->p[ip] < p1))
	    || (p0 < p1 && (atm->p[ip] > p0 && atm->p[ip] < p1)))
	  continue;

      /* Check potential temperature... */
      if (theta0 != theta1) {
	double theta;
	if (ctl.qnt_theta >= 0)
	  theta = atm->q[ctl.qnt_theta][ip];
	else if (ctl.qnt_t >= 0)
	  theta = THETA(atm->p[ip], atm->q[ctl.qnt_t][ip]);
	else
	  ERRMSG
	    ("Filtering requires temperature or potential temperature data!");
	if ((theta1 > theta0 && (theta < theta0 || theta > theta1))
	    || (theta1 < theta0 && (theta < theta0 && theta > theta1)))
	  continue;
      }

      /* Check longitude... */
      if (lon0 != lon1)
	if ((lon1 > lon0 && (atm->lon[ip] < lon0 || atm->lon[ip] > lon1))
	    || (lon1 < lon0 && (atm->lon[ip] < lon0 && atm->lon[ip] > lon1)))
	  continue;

      /* Check latitude... */
      if (lat0 != lat1)
	if ((lat1 > lat0 && (atm->lat[ip] < lat0 || atm->lat[ip] > lat1))
	    || (lat1 < lat0 && (atm->lat[ip] < lat0 && atm->lat[ip] > lat1)))
	  continue;

      /* Check horizontal distace... */
      if (r0 != r1) {
	geo2cart(0, atm->lon[ip], atm->lat[ip], x1);
	double r = DIST(x0, x1);
	if ((r1 > r0 && (r < r0 || r > r1))
	    || (r1 < r0 && (r < r0 && r > r1)))
	  continue;
      }

      /* Copy data... */
      atm2->time[atm2->np] = atm->time[ip];
      atm2->p[atm2->np] = atm->p[ip];
      atm2->lon[atm2->np] = atm->lon[ip];
      atm2->lat[atm2->np] = atm->lat[ip];
      for (int iq = 0; iq < ctl.nq; iq++)
	atm2->q[iq][atm2->np] = atm->q[iq][ip];
      if ((++atm2->np) > NP)
	ERRMSG("Too many air parcels!");
    }
  }

  /* Close file... */
  write_atm(argv[2], &ctl, atm2, 0);

  /* Free... */
  free(atm);
  free(atm2);

  return EXIT_SUCCESS;
}
