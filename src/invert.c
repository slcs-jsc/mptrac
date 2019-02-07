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
  
  Copyright (C) 2013-2019 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Calculate emission estimates.
*/

#include "libtrac.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Maximum number of observations. */
#define MMAX 100000

/*! Maximum number of parameters. */
#define NMAX 1000

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Add observations. */
void add_obs(
  double rt_old,
  double rt_all[MMAX],
  double rlon_old,
  double rlon_all[MMAX],
  double rlat_old,
  double rlat_all[MMAX],
  double robs_old,
  double robs_all[MMAX],
  double rsig_old,
  double rsig_all[MMAX],
  double rsim_old,
  double rsim_all[MMAX][NMAX],
  int *m,
  int *n);

/* Convolute vmr profile with kernel. */
double conv_kernel(
  const char *filename,
  double *rz,
  double *rvmr,
  int nz);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  static FILE *in, *in2, *out;

  static char line[LEN], filename[LEN], wrkdir[NMAX][LEN];

  static double rt, rt_old, rt_all[MMAX], rz[GZ], rlon, rlon_old,
    rlon_all[MMAX], rlat, rlat_old, rlat_all[MMAX], rvmr[GZ], rdum, robs,
    robs_old, robs_all[MMAX], rsig, rsig_old, rsig_all[MMAX],
    rsim_old, rsim_all[MMAX][NMAX], sigma, chisq, tol,
    rcond, rnorm, snorm, lambda, G_gcv;

  static int i, j, m, n, nz, method, offset;

  static size_t rank, npoints = 200, reg_idx;

  gsl_multifit_linear_workspace *W;

  gsl_matrix *C, *X;

  gsl_vector *c, *r, *w, *y, *reg_param, *rho, *eta, *G;

  /* Check arguments... */
  if (argc < 8)
    ERRMSG
      ("Give parameters: <ctl> <dirlist> <prof.tab> <kernel.tab>"
       " <fit.tab> <corr.tab> <res.tab>");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  method =
    (int) scan_ctl(argv[1], argc, argv, "INVERT_METHOD", -1, "2", NULL);
  offset =
    (int) scan_ctl(argv[1], argc, argv, "INVERT_OFFSET", -1, "1", NULL);
  sigma = scan_ctl(argv[1], argc, argv, "INVERT_SIGMA", -1, "0.5", NULL);
  tol = scan_ctl(argv[1], argc, argv, "INVERT_TOL", -1, "0.001", NULL);

  /* Open directory list... */
  if (!(in = fopen(argv[2], "r")))
    ERRMSG("Cannot open directory list!");

  /* Loop over directories... */
  while (fscanf(in, "%s", wrkdir[n]) != EOF) {

    /* Open profile data file... */
    sprintf(filename, "%s/%s", wrkdir[n], argv[3]);
    printf("Read profile data: %s\n", filename);
    if (!(in2 = fopen(filename, "r")))
      ERRMSG("Cannot open file!");

    /* Read profile data... */
    nz = 0;
    while (fgets(line, LEN, in2))
      if (sscanf(line, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
		 &rt, &rz[nz], &rlon, &rlat, &rdum, &rdum,
		 &rvmr[nz], &rdum, &rdum, &robs, &rsig) == 11) {
	if (nz > 0 && (rt != rt_old || rlon != rlon_old || rlat != rlat_old)) {

	  /* Convolute with kernel... */
	  rsim_old = conv_kernel(argv[4], rz, rvmr, nz);

	  /* Add observations... */
	  add_obs(rt_old, rt_all, rlon_old, rlon_all, rlat_old,
		  rlat_all, robs_old, robs_all, rsig_old, rsig_all,
		  rsim_old, rsim_all, &m, &n);

	  /* Reset profile counter... */
	  rz[0] = rz[nz];
	  rvmr[0] = rvmr[nz];
	  nz = 1;

	} else if ((++nz) >= GZ)
	  ERRMSG("Too many levels!");

	/* Save old position... */
	rt_old = rt;
	rlon_old = rlon;
	rlat_old = rlat;
	robs_old = robs;
	rsig_old = rsig;
      }

    /* Add last profile... */
    if (nz > 0) {

      /* Convolute with kernel... */
      rsim_old = conv_kernel(argv[4], rz, rvmr, nz);

      /* Add observations... */
      add_obs(rt_old, rt_all, rlon_old, rlon_all, rlat_old,
	      rlat_all, robs_old, robs_all, rsig_old, rsig_all,
	      rsim_old, rsim_all, &m, &n);
    }

    /* Increment parameter counter... */
    if ((++n) >= NMAX)
      ERRMSG("Too many parameters!");

    /* Close file... */
    fclose(in2);
  }

  /* Close dirlist... */
  fclose(in);

  /* Add offset... */
  if (offset) {
    if ((++n) >= NMAX)
      ERRMSG("Too many parameters!");
    sprintf(wrkdir[n - 1], "OFFSET");
  }

  /* Write info... */
  printf("Calculate least square fit...\n");
  printf("  m= %d\n", m);
  printf("  n= %d\n", n);

  /* Allocate... */
  W = gsl_multifit_linear_alloc((size_t) m, (size_t) n);
  C = gsl_matrix_alloc((size_t) n, (size_t) n);
  X = gsl_matrix_alloc((size_t) m, (size_t) n);
  c = gsl_vector_alloc((size_t) n);
  r = gsl_vector_alloc((size_t) m);
  w = gsl_vector_alloc((size_t) m);
  y = gsl_vector_alloc((size_t) m);
  reg_param = gsl_vector_alloc(npoints);
  rho = gsl_vector_alloc(npoints);
  eta = gsl_vector_alloc(npoints);
  G = gsl_vector_alloc(npoints);

  /* Set weights, observation vector, and prediction matrix... */
  for (i = 0; i < m; i++) {
    gsl_vector_set(w, (size_t) i,
		   (!gsl_finite(rsig_all[i]) || rsig_all[i] < sigma)
		   ? 1. / gsl_pow_2(sigma) : 1. / gsl_pow_2(rsig_all[i]));
    gsl_vector_set(y, (size_t) i, robs_all[i]);
    for (j = 0; j < n; j++)
      if (offset && j == n - 1)
	gsl_matrix_set(X, (size_t) i, (size_t) j, 1.0);
      else
	gsl_matrix_set(X, (size_t) i, (size_t) j, rsim_all[i][j]);
  }

  /* Truncated SVD... */
  if (method == 1) {

    /* Compute fit... */
    gsl_multifit_wlinear_tsvd(X, w, y, tol, c, C, &chisq, &rank, W);

    /* Write info... */
    printf("Truncated SVD:");
    printf("  chisq / dof = %g\n", chisq / (m - n));
    printf("  effective rank = %lu\n", rank);
  }

  /* Regularized LSQ... */
  else if (method == 2) {

    /* Compute SVD of X... */
    printf("Compute SVD...\n");
    gsl_multifit_linear_svd(X, W);

    /* Get condition number... */
    rcond = gsl_multifit_linear_rcond(W);
    printf("  matrix condition number = %g\n", 1.0 / rcond);

    /* Unregularized LSQ (lambda = 0)... */
    gsl_multifit_linear_solve(0.0, X, y, c, &rnorm, &snorm, W);
    chisq = pow(rnorm, 2.0);

    /* Write info... */
    printf("Unregularized LSQ:");
    printf("  optimal lambda = %g\n", lambda);
    printf("  residual norm = %g\n", rnorm);
    printf("  solution norm = %g\n", snorm);
    printf("  chisq / dof = %g\n", chisq / (m - n));

    /* Calculate L-curve and find its corner... */
    gsl_multifit_linear_lcurve(y, reg_param, rho, eta, W);
    gsl_multifit_linear_lcorner(rho, eta, &reg_idx);

    /* Store optimal regularization parameter... */
    lambda = gsl_vector_get(reg_param, reg_idx);

    /* Regularize with lambda_l... */
    gsl_multifit_linear_solve(lambda, X, y, c, &rnorm, &snorm, W);
    chisq = pow(rnorm, 2.0) + pow(lambda * snorm, 2.0);

    /* Write info... */
    printf("Regularization with L-curve:");
    printf("  optimal lambda = %g\n", lambda);
    printf("  residual norm = %g\n", rnorm);
    printf("  solution norm = %g\n", snorm);
    printf("  chisq / dof = %g\n", chisq / (m - n));
    for (i = 0; i < (int) npoints; ++i)
      printf("    l_curve: reg_param= %g rho= %g eta= %g\n",
	     gsl_vector_get(reg_param, (size_t) i),
	     gsl_vector_get(rho, (size_t) i),
	     gsl_vector_get(eta, (size_t) i));
    printf("    l_corner_point: rho= %g eta= %g\n",
	   gsl_vector_get(rho, reg_idx), gsl_vector_get(eta, reg_idx));

    /* Calculate GCV curve and find its minimum... */
    gsl_multifit_linear_gcv(y, reg_param, G, &lambda, &G_gcv, W);

    /* Regularize with lambda_gcv... */
    gsl_multifit_linear_solve(lambda, X, y, c, &rnorm, &snorm, W);
    chisq = pow(rnorm, 2.0) + pow(lambda * snorm, 2.0);

    /* Write info... */
    printf("Regularization with GCV curve:");
    printf("  optimal lambda = %g\n", lambda);
    printf("  residual norm = %g\n", rnorm);
    printf("  solution norm = %g\n", snorm);
    printf("  chisq / dof = %g\n", chisq / (m - n));
    for (i = 0; i < (int) npoints; ++i)
      printf("    gcv_curve: reg_param= %g G= %g\n",
	     gsl_vector_get(reg_param, (size_t) i),
	     gsl_vector_get(G, (size_t) i));
    printf("    gcv_minimum: lambda= %g G= %g\n", lambda, G_gcv);
  }

  /* Error... */
  else
    ERRMSG("Unknown inversion method!");

  /* Write fit parameters... */
  if (argv[5][0] != '-') {

    /* Create file... */
    printf("Write fit parameters: %s\n", argv[5]);
    if (!(out = fopen(argv[5], "w")))
      ERRMSG("Cannot create file!");

    /* Write header... */
    fprintf(out,
	    "# $1 = parameter index\n"
	    "# $2 = parameter name\n"
	    "# $3 = parameter best fit\n"
	    "# $4 = parameter error\n"
	    "# $5 = chisq / (m - n)\n" "# $6 = rank\n\n");
    for (i = 0; i < n; i++)
      fprintf(out, "%d %s %g %g %g %lu\n", i, wrkdir[i],
	      gsl_vector_get(c, (size_t) i),
	      sqrt(gsl_matrix_get(C, (size_t) i, (size_t) i)),
	      chisq / (double) (m - n), rank);

    /* Close file... */
    fclose(out);
  }

  /* Write correlation matrix... */
  if (argv[6][0] != '-') {

    /* Create file... */
    printf("Write correlation matrix: %s\n", argv[6]);
    if (!(out = fopen(argv[6], "w")))
      ERRMSG("Cannot create file!");

    /* Write header... */
    fprintf(out,
	    "# $1 = row index\n"
	    "# $2 = row parameter name\n"
	    "# $3 = column index\n"
	    "# $4 = column parameter name\n"
	    "# $5 = correlation coefficient\n");
    for (i = 0; i < n; i++) {
      fprintf(out, "\n");
      for (j = 0; j < n; j++)
	fprintf(out, "%d %s %d %s %g\n",
		i, wrkdir[i], j, wrkdir[j],
		gsl_matrix_get(C, (size_t) i, (size_t) j)
		/ sqrt(gsl_matrix_get(C, (size_t) i, (size_t) i))
		/ sqrt(gsl_matrix_get(C, (size_t) j, (size_t) j)));
    }
  }

  /* Write residuals... */
  if (argv[7][0] != '-') {

    /* Calculate residuals... */
    gsl_multifit_linear_residuals(X, y, c, r);

    /* Create file... */
    printf("Write residuals: %s\n", argv[7]);
    if (!(out = fopen(argv[7], "w")))
      ERRMSG("Cannot create file!");

    /* Write header... */
    fprintf(out,
	    "# $1 = observation index\n"
	    "# $2 = observation time [s]\n"
	    "# $3 = observation longitude [deg]\n"
	    "# $4 = observation latitude [deg]\n"
	    "# $5 = observation [K]\n"
	    "# $6 = observation error [K]\n"
	    "# $7 = residual (obs - sim) [K]\n\n");

    /* Write data... */
    for (i = 0; i < m; i++)
      fprintf(out, "%d %.2f %g %g %g %g %g\n",
	      i, rt_all[i], rlon_all[i], rlat_all[i], robs_all[i],
	      rsig_all[i], gsl_vector_get(r, (size_t) i));

    /* Close file... */
    fclose(out);
  }

  /* Free... */
  gsl_multifit_linear_free(W);
  gsl_matrix_free(C);
  gsl_matrix_free(X);
  gsl_vector_free(c);
  gsl_vector_free(r);
  gsl_vector_free(w);
  gsl_vector_free(y);
  gsl_vector_free(reg_param);
  gsl_vector_free(rho);
  gsl_vector_free(eta);
  gsl_vector_free(G);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void add_obs(
  double rt_old,
  double rt_all[MMAX],
  double rlon_old,
  double rlon_all[MMAX],
  double rlat_old,
  double rlat_all[MMAX],
  double robs_old,
  double robs_all[MMAX],
  double rsig_old,
  double rsig_all[MMAX],
  double rsim_old,
  double rsim_all[MMAX][NMAX],
  int *m,
  int *n) {

  int i;

  /* Find observation index... */
  for (i = 0; i < *m; i++)
    if (rt_old == rt_all[i]
	&& rlon_old == rlon_all[i]
	&& rlat_old == rlat_all[i])
      break;
  if (i >= *m)
    if ((++(*m)) >= MMAX)
      ERRMSG("Too many observations!");

  /* Save data... */
  rt_all[i] = rt_old;
  rlon_all[i] = rlon_old;
  rlat_all[i] = rlat_old;
  robs_all[i] = robs_old;
  rsig_all[i] = rsig_old;
  rsim_all[i][*n] = rsim_old;
}

/*****************************************************************************/

double conv_kernel(
  const char *filename,
  double *rz,
  double *rvmr,
  int nz) {

  static FILE *in;

  static char line[LEN];

  static double sum, z[GZ], k[GZ];

  static int init, ik, iz, nk;

  /* Init... */
  if (!init) {
    init = 1;

    /* Read kernel file... */
    printf("Read kernel file: %s\n", filename);
    if (!(in = fopen(filename, "r")))
      ERRMSG("Cannot open file!");
    while (fgets(line, LEN, in))
      if (sscanf(line, "%lg %lg", &z[nk], &k[nk]) == 2)
	if ((++nk) >= GZ)
	  ERRMSG("Too many levels!");
    fclose(in);
  }

  /* Convolute with kernel... */
  sum = 0;
  for (ik = 0; ik < nk; ik++) {
    iz = locate(rz, nz, z[ik]);
    sum += k[ik]
      * LIN(rz[iz], rvmr[iz], rz[iz + 1], rvmr[iz + 1], z[ik]);
  }
  return sum;
}
