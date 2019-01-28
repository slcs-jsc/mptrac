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

/*! Maximum number of state variables. */
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
    robs_old, robs_all[MMAX], rsim_old, rsim_all[MMAX][NMAX], sigma, chisq,
    tol;

  static int i, j, m, n, nz, offset;

  static size_t rank;

  gsl_multifit_linear_workspace *W;

  gsl_matrix *C, *X;

  gsl_vector *c, *r, *w, *y;

  /* Check arguments... */
  if (argc < 8)
    ERRMSG
      ("Give parameters: <ctl> <dirlist> <prof.tab> <kernel.tab>"
       " <fit.tab> <corr.tab> <res.tab>");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);
  offset = (int)scan_ctl(argv[1], argc, argv, "INVERT_OFFSET", -1, "1", NULL);
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
      if (sscanf(line, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
		 &rt, &rz[nz], &rlon, &rlat, &rdum,
		 &rdum, &rvmr[nz], &rdum, &rdum, &robs) == 10) {
	if (nz > 0 && (rt != rt_old || rlon != rlon_old || rlat != rlat_old)) {

	  /* Convolute with kernel... */
	  rsim_old = conv_kernel(argv[4], rz, rvmr, nz);

	  /* Add observations... */
	  add_obs(rt_old, rt_all, rlon_old, rlon_all, rlat_old,
		  rlat_all, robs_old, robs_all, rsim_old, rsim_all,
		  &m, &n);
	  
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
      }

    /* Add last profile... */
    if (nz > 0) {

      /* Convolute with kernel... */
      rsim_old = conv_kernel(argv[4], rz, rvmr, nz);

      /* Add observations... */
      add_obs(rt_old, rt_all, rlon_old, rlon_all, rlat_old, rlat_all,
	      robs_old, robs_all, rsim_old, rsim_all, &m, &n);
    }

    /* Increment stae vector counter... */
    if ((++n) >= NMAX)
      ERRMSG("Too many state vector elements!");

    /* Close file... */
    fclose(in2);
  }

  /* Close dirlist... */
  fclose(in);

  /* Add offset... */
  if(offset) {
    if ((++n) >= NMAX)
      ERRMSG("Too many state vector elements!");
    sprintf(wrkdir[n-1], "OFFSET");
  }
  
  /* Write info... */
  printf("Calculate least square fit (m= %d, n= %d)...\n", m, n);
  
  /* Allocate... */
  W = gsl_multifit_linear_alloc((size_t) m, (size_t) n);
  C = gsl_matrix_alloc((size_t) n, (size_t) n);
  X = gsl_matrix_alloc((size_t) m, (size_t) n);
  c = gsl_vector_alloc((size_t) n);
  r = gsl_vector_alloc((size_t) m);
  w = gsl_vector_alloc((size_t) m);
  y = gsl_vector_alloc((size_t) m);

  /* Set prediction matrix and observation vector... */
  for (i = 0; i < m; i++) {
    gsl_vector_set(w, (size_t) i, sigma);
    gsl_vector_set(y, (size_t) i, robs_all[i]);
    for (j = 0; j < n; j++)
      if(offset && j==n-1)
	gsl_matrix_set(X, (size_t) i, (size_t) j, 1.0);
      else
	gsl_matrix_set(X, (size_t) i, (size_t) j, rsim_all[i][j]);
  }
  
  /* Least square fit... */
  if(sigma>0) {
    if(tol>0)
      gsl_multifit_wlinear_tsvd(X, w, y, tol, c, C, &chisq, &rank, W);
    else
      gsl_multifit_wlinear(X, w, y, c, C, &chisq, W);
  } else {
    if(tol>0)
      gsl_multifit_linear_tsvd(X, y, tol, c, C, &chisq, &rank, W);
    else
      gsl_multifit_linear(X, y, c, C, &chisq, W);
  }
  
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
	    "# $5 = chisq / m\n" "# $6 = rank\n\n");
    for (i = 0; i < n; i++)
      fprintf(out, "%d %s %g %g %g %lu\n", i, wrkdir[i],
	      gsl_vector_get(c, (size_t) i),
	      sqrt(gsl_matrix_get(C, (size_t) i, (size_t) i)),
	      chisq / (double) m, rank);

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
	    "# $5 = observation [K]\n" "# $6 = residual [K]\n\n");

    /* Write data... */
    for (i = 0; i < m; i++)
      fprintf(out, "%d %.2f %g %g %g %g\n",
	      i, rt_all[i], rlon_all[i], rlat_all[i], robs_all[i],
	      gsl_vector_get(r, (size_t) i));

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
  double rsim_old,
  double rsim_all[MMAX][NMAX],
  int *m,
  int *n) {

  int i;

  /* Find observation index... */
  for (i = 0; i < *m; i++)
    if (rt_old == rt_all[i]
	&& rlon_old == rlon_all[i]
	&& rlat_old == rlat_all[i]
	&& robs_old == robs_all[i])
      break;
  if (i >= *m)
    if ((++(*m)) >= MMAX)
      ERRMSG("Too many observations!");

  /* Save data... */
  rt_all[i] = rt_old;
  rlon_all[i] = rlon_old;
  rlat_all[i] = rlat_old;
  robs_all[i] = robs_old;
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
