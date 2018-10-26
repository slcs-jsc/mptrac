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
  Clustering of trajectories.
*/

#include "libtrac.h"

/* ------------------------------------------------------------
   Defines...
   ------------------------------------------------------------ */

/*! Maximum number of seeds. */
#define NS 100

/*! Maximum number of timesteps. */
#define NT 1000

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  atm_t *atm;

  gsl_rng *rng;

  FILE *out;

  static double d2, *dist, lat, lon, rmsd[NS],
    x[3], xs[NT][NS][3], z, zs[NT][NS];

  static int *cluster, f, idx[NS], ip, is, it, itmax, np[NS], ns;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <cluster.log> <atm1> [<atm2> <atm3> ...]");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  ns = (int) scan_ctl(argv[1], argc, argv, "CLUSTER_NS", -1, "7", NULL);
  if (ns > NS)
    ERRMSG("Too many seeds!");
  itmax =
    (int) scan_ctl(argv[1], argc, argv, "CLUSTER_ITMAX", -1, "10", NULL);

  /* Initialize random number generator... */
  gsl_rng_env_setup();
  rng = gsl_rng_alloc(gsl_rng_default);

  /* Allocate... */
  ALLOC(atm, atm_t, 1);
  ALLOC(cluster, int,
	NP);
  ALLOC(dist, double,
	NP * NS);

  /* Create output file... */
  printf("Write cluster data: %s\n", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = iteration index\n"
	  "# $2 = seed index\n"
	  "# $3 = time step index\n"
	  "# $4 = mean altitude [km]\n"
	  "# $5 = mean longitude [deg]\n"
	  "# $6 = mean latitude [deg]\n"
	  "# $7 = number of points\n" "# $8 = RMSD [km^2]\n");

  /* Get seeds (random selection of trajectories)... */
  for (f = 3; f < argc; f++) {

    /* Check number of timesteps... */
    if (f - 3 > NT)
      ERRMSG("Too many timesteps!");

    /* Read atmopheric data... */
    read_atm(argv[f], &ctl, atm);

    /* Pick seeds (random selection)... */
    if (f == 3)
      for (is = 0; is < ns; is++)
	idx[is] = (int) gsl_rng_uniform_int(rng, (long unsigned int) atm->np);

    /* Save seeds... */
    for (is = 0; is < ns; is++) {
      geo2cart(0, atm->lon[idx[is]], atm->lat[idx[is]], xs[f - 3][is]);
      zs[f - 3][is] = Z(atm->p[idx[is]]);
    }
  }

  /* Iterations... */
  for (it = 0; it < itmax; it++) {

    /* Write output... */
    for (is = 0; is < ns; is++) {
      fprintf(out, "\n");
      for (f = 3; f < argc; f++) {
	cart2geo(xs[f - 3][is], &z, &lon, &lat);
	fprintf(out, "%d %d %d %g %g %g %d %g\n",
		it, is, f - 3, zs[f - 3][is], lon, lat, np[is], rmsd[is]);
      }
    }

    /* Init... */
    for (ip = 0; ip < atm->np; ip++)
      for (is = 0; is < ns; is++) {
	dist[ip * NS + is] = 0;
	rmsd[is] = 0;
      }

    /* Get distances between seeds and trajectories... */
    for (f = 3; f < argc; f++) {

      /* Read atmopheric data... */
      read_atm(argv[f], &ctl, atm);

      /* Get distances... */
      for (ip = 0; ip < atm->np; ip++) {
	geo2cart(0, atm->lon[ip], atm->lat[ip], x);
	z = Z(atm->p[ip]);
	for (is = 0; is < ns; is++) {
	  d2 =
	    DIST2(x, xs[f - 3][is]) + gsl_pow_2((z - zs[f - 3][is]) * 200.);
	  dist[ip * NS + is] += d2;
	  rmsd[is] += d2;
	}
      }
    }

    /* Assign clusters... */
    for (ip = 0; ip < atm->np; ip++)
      cluster[ip] = (int) gsl_stats_min_index(&dist[ip * NS], 1, (size_t) ns);

    /* Recalculate seeds (mean trajectories)... */
    for (f = 3; f < argc; f++) {

      /* Read atmopheric data... */
      read_atm(argv[f], &ctl, atm);

      /* Calculate new seeds... */
      for (is = 0; is < ns; is++) {
	xs[f - 3][is][0] = 0;
	xs[f - 3][is][1] = 0;
	xs[f - 3][is][2] = 0;
	zs[f - 3][is] = 0;
	np[is] = 0;
      }
      for (ip = 0; ip < atm->np; ip++) {
	geo2cart(0, atm->lon[ip], atm->lat[ip], x);
	xs[f - 3][cluster[ip]][0] += x[0];
	xs[f - 3][cluster[ip]][1] += x[1];
	xs[f - 3][cluster[ip]][2] += x[2];
	zs[f - 3][cluster[ip]] += Z(atm->p[ip]);
	np[cluster[ip]]++;
      }
      for (is = 0; is < ns; is++) {
	xs[f - 3][is][0] /= np[is];
	xs[f - 3][is][1] /= np[is];
	xs[f - 3][is][2] /= np[is];
	zs[f - 3][is] /= np[is];
      }
    }
  }

  /* Write output... */
  for (is = 0; is < ns; is++) {
    fprintf(out, "\n");
    for (f = 3; f < argc; f++) {
      cart2geo(xs[f - 3][is], &z, &lon, &lat);
      fprintf(out, "%d %d %d %g %g %g %d %g\n",
	      it, is, f - 3, zs[f - 3][is], lon, lat, np[is], rmsd[is]);
    }
  }

  /* Close output file... */
  fclose(out);

  /* Write clustering results... */
  if (ctl.qnt_ens >= 0)

    /* Recalculate seeds (mean trajectories)... */
    for (f = 3; f < argc; f++) {

      /* Read atmopheric data... */
      read_atm(argv[f], &ctl, atm);

      /* Set ensemble ID... */
      for (ip = 0; ip < atm->np; ip++)
	atm->q[ctl.qnt_ens][ip] = cluster[ip];

      /* Write atmospheric data... */
      write_atm(argv[f], &ctl, atm, 0);
    }

  /* Free... */
  gsl_rng_free(rng);
  free(atm);
  free(cluster);
  free(dist);

  return EXIT_SUCCESS;
}
