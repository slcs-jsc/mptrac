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
  
  Copyright (C) 2013-2026 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Check model time step based on given meteorological data.
*/

#include "mptrac.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Compute PBL weight for a pressure level. */
static double pbl_weight_dt(
  const double p,
  const double pbl,
  const double ps,
  const double pbl_trans);

/*! Compute tropospheric weight for a pressure level. */
static double tropo_weight_dt(
  const double p,
  const double pt);

/*! Print command-line help. */
void usage(
  void);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  clim_t *clim;

  met_t *met;

  dd_t *dd;

  /* Print usage information... */
  USAGE;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: met_check_dt <ctl> <dt_file> <met> [KEY VALUE ...]\n\n"
	   "Use -h for full help.");

  /* Allocate... */
  ALLOC(clim, clim_t, 1);
  ALLOC(met, met_t, 1);
  ALLOC(dd, dd_t, 1);

  /* Read control parameters... */
  mptrac_read_ctl(argv[1], argc, argv, &ctl);
  const double pbl_trans = ctl.turb_pbl_trans;
  const double dx = 1e3 * scan_ctl(argv[1], argc, argv, "DX", -1, "", NULL);
  const double c_max = scan_ctl(argv[1], argc, argv, "CMAX", -1, "0.5", NULL);
  const double n_max = scan_ctl(argv[1], argc, argv, "NMAX", -1, "0.3", NULL);

  /* Read climatological data... */
  mptrac_read_clim(&ctl, clim);

  /* Read meteo data... */
  if (!mptrac_read_met(argv[3], &ctl, clim, met, dd))
    ERRMSG("Cannot open file!");

  /* Create output file... */
  FILE *out;
  LOG(1, "Write time step data file: %s", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = height [km]\n"
	  "# $2 = time step for horizontal advection [s]\n"
	  "# $3 = time step for vertical advection [s]\n"
	  "# $4 = time step for horizontal diffusion [s]\n"
	  "# $5 = time step for vertical diffusion [s]\n"
	  "# $6 = time step for PBL transition diffusion [s]\n"
	  "# $7 = time step for PBL depth diffusion [s]\n\n");

  /* Loop over pressure levels... */
  for (int ip = 1; ip < met->np - 1; ip++) {

    /* Init... */
    double dt_x_min = 1e100;
    double dt_p_min = 1e100;
    double dt_dx_min = 1e100;
    double dt_dp_min = 1e100;
    double dt_pbl_min = 1e100;
    double dt_pbl_depth_min = 1e100;

    /* Loop over columns... */
#pragma omp parallel for default(shared) collapse(2) reduction(min:dt_x_min,dt_p_min,dt_dx_min,dt_dp_min,dt_pbl_min,dt_pbl_depth_min)
    for (int ix = 0; ix < met->nx; ix++)
      for (int iy = 1; iy < met->ny - 1; iy++) {

	/* Check advection... */
	const double vh =
	  sqrt(SQR(met->u[ix][iy][ip]) + SQR(met->v[ix][iy][ip]));
	const double dt_x = fabs(c_max * dx / vh);
	if (vh != 0)
	  dt_x_min = MIN(dt_x, dt_x_min);

	const double dp = 0.5 * fabs(met->p[ip + 1] - met->p[ip - 1]);
	const double dt_p = fabs(c_max * dp / met->w[ix][iy][ip]);
	if (met->w[ix][iy][ip] != 0)
	  dt_p_min = MIN(dt_p, dt_p_min);

	/* Get layer weights... */
	const double pt = clim_tropo(clim, met->time,
				     ctl.met_coord_type ==
				     0 ? met->lat[iy] : ctl.met_utm_ref_lat);
	const double wpbl =
	  pbl_weight_dt(met->p[ip], met->pbl[ix][iy], met->ps[ix][iy],
			pbl_trans);
	const double wtrop = tropo_weight_dt(met->p[ip], pt) * (1.0 - wpbl);
	const double wstrat = 1.0 - wpbl - wtrop;

	/* Check layer-aware diffusion... */
	const double dx_loc = wpbl * ctl.turb_dx_pbl
	  + wtrop * ctl.turb_dx_trop + wstrat * ctl.turb_dx_strat;
	if (dx_loc > 0) {
	  const double dt_dx = 0.5 * SQR(n_max * dx) / dx_loc;
	  dt_dx_min = MIN(dt_dx, dt_dx_min);
	}

	const double dz_loc = wpbl * ctl.turb_dz_pbl
	  + wtrop * ctl.turb_dz_trop + wstrat * ctl.turb_dz_strat;
	if (dz_loc > 0) {
	  const double dt_dp = 0.5 * SQR(n_max * dp)
	    / (SQR(met->p[ip] / (100. * H0)) * dz_loc);
	  dt_dp_min = MIN(dt_dp, dt_dp_min);
	}

	/* Check PBL transition diffusion... */
	if (pbl_trans > 0 && met->ps[ix][iy] > met->pbl[ix][iy]
	    && met->p[ip] <= met->pbl[ix][iy]
	    && met->p[ip] >= met->pbl[ix][iy]
	    - pbl_trans * (met->ps[ix][iy] - met->pbl[ix][iy])) {
	  const double p0 = met->pbl[ix][iy];
	  const double p1 = p0 - pbl_trans * (met->ps[ix][iy] - p0);
	  const double dz_trans = 1e3 * fabs(Z(p1) - Z(p0));
	  const double dz_max = MAX(ctl.turb_dz_pbl, ctl.turb_dz_trop);
	  if (dz_trans > 0 && dz_max > 0) {
	    const double dt_trans = 0.5 * SQR(n_max * dz_trans) / dz_max;
	    dt_pbl_min = MIN(dt_trans, dt_pbl_min);
	  }
	}

	/* Check PBL depth diffusion on the full PBL scale. */
	if (met->ps[ix][iy] > met->pbl[ix][iy]
	    && met->p[ip] >= met->pbl[ix][iy]) {
	  const double dz_pbl =
	    1e3 * fabs(Z(met->pbl[ix][iy]) - Z(met->ps[ix][iy]));
	  if (dz_pbl > 0 && ctl.turb_dz_pbl > 0) {
	    const double dt_pbl_depth =
	      0.5 * SQR(n_max * dz_pbl) / ctl.turb_dz_pbl;
	    dt_pbl_depth_min = MIN(dt_pbl_depth, dt_pbl_depth_min);
	  }
	}
      }

    /* Write output... */
    const double out_dt_dx = dt_dx_min < 1e99 ? dt_dx_min : NAN;
    const double out_dt_dp = dt_dp_min < 1e99 ? dt_dp_min : NAN;
    const double out_dt_pbl = dt_pbl_min < 1e99 ? dt_pbl_min : NAN;
    const double out_dt_pbl_depth =
      dt_pbl_depth_min < 1e99 ? dt_pbl_depth_min : NAN;
    fprintf(out, "%g %g %g %g %g %g %g\n", Z(met->p[ip]), dt_x_min, dt_p_min,
	    out_dt_dx, out_dt_dp, out_dt_pbl, out_dt_pbl_depth);
  }

  /* Close output file... */
  fclose(out);

  /* Free... */
  free(clim);
  free(met);
  free(dd);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

/*! Print command-line help. */
void usage(
  void) {

  printf("\nMPTRAC met_check_dt tool.\n\n");
  printf("Check model time-step constraints for meteorological data.\n");
  printf("\n");
  printf("Usage:\n");
  printf("  met_check_dt <ctl> <dt_file> <met> [KEY VALUE ...]\n");
  printf("\n");
  printf("Arguments:\n");
  printf("  <ctl>      Control file.\n");
  printf("  <dt_file>  Output table for time-step diagnostics.\n");
  printf("  <met>      Meteorological input file.\n");
  printf("  [KEY VALUE]  Optional control parameters.\n");
  printf("\nFurther information:\n");
  printf("  Manual: https://slcs-jsc.github.io/mptrac/\n");
}

/*****************************************************************************/

static double pbl_weight_dt(
  const double p,
  const double pbl,
  const double ps,
  const double pbl_trans) {

  const double p0 = pbl;
  const double p1 = pbl - pbl_trans * (ps - pbl);

  if (p > p0)
    return 1.0;
  else if (p < p1)
    return 0.0;
  else
    return LIN(p0, 1.0, p1, 0.0, p);
}

/*****************************************************************************/

static double tropo_weight_dt(
  const double p,
  const double pt) {

  const double p1 = pt * 0.866877899;
  const double p0 = pt / 0.866877899;

  if (p > p0)
    return 1.0;
  else if (p < p1)
    return 0.0;
  else
    return LIN(p0, 1.0, p1, 0.0, p);
}
