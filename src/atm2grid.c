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
  Convert atmosphric data file to grid data file.
*/

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  atm_t *atm;

  static double kz[EP], kw[EP];

  static int nk;

  double *cd, *mean[NQ], *stddev[NQ], *vmr_impl, *z, *lon, *lat, *area,
    *press;

  int *ixs, *iys, *izs, *np;
  
  /* Allocate... */
  ALLOC(atm, atm_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <atm_in> <grid_out>");

  read_ctl(argv[1], argc, argv, &ctl);

  /* Read atmospheric data... */
  if (!read_atm(argv[2], &ctl, atm))
    ERRMSG("Cannot open file!");

  /* Write info... */
  LOG(1, "Write grid data: %s", argv[3]);


  /* Read kernel data... */
  if (ctl.grid_kernel[0] != '-')
    read_kernel(ctl.grid_kernel, kz, kw, &nk);

  /* Allocate... */
  ALLOC(cd, double,
	ctl.grid_nx * ctl.grid_ny * ctl.grid_nz);
  for (int iq = 0; iq < ctl.nq; iq++) {
    ALLOC(mean[iq], double,
	  ctl.grid_nx * ctl.grid_ny * ctl.grid_nz);
    ALLOC(stddev[iq], double,
	  ctl.grid_nx * ctl.grid_ny * ctl.grid_nz);
  }
  ALLOC(vmr_impl, double,
	ctl.grid_nx * ctl.grid_ny * ctl.grid_nz);
  ALLOC(z, double,
	ctl.grid_nz);
  ALLOC(lon, double,
	ctl.grid_nx);
  ALLOC(lat, double,
	ctl.grid_ny);
  ALLOC(area, double,
	ctl.grid_ny);
  ALLOC(press, double,
	ctl.grid_nz);
  ALLOC(np, int,
	ctl.grid_nx * ctl.grid_ny * ctl.grid_nz);
  ALLOC(ixs, int,
	atm->np);
  ALLOC(iys, int,
	atm->np);
  ALLOC(izs, int,
	atm->np);

  /* Set grid box size... */
  double dz = (ctl.grid_z1 - ctl.grid_z0) / ctl.grid_nz;
  double dlon = (ctl.grid_lon1 - ctl.grid_lon0) / ctl.grid_nx;
  double dlat = (ctl.grid_lat1 - ctl.grid_lat0) / ctl.grid_ny;

  /* Set vertical coordinates... */
  for (int iz = 0; iz < ctl.grid_nz; iz++) {
    z[iz] = ctl.grid_z0 + dz * (iz + 0.5);
    press[iz] = P(z[iz]);
  }

  /* Set horizontal coordinates... */
  for (int ix = 0; ix < ctl.grid_nx; ix++)
    lon[ix] = ctl.grid_lon0 + dlon * (ix + 0.5);
  for (int iy = 0; iy < ctl.grid_ny; iy++) {
    lat[iy] = ctl.grid_lat0 + dlat * (iy + 0.5);
    area[iy] = dlat * dlon * SQR(RE * M_PI / 180.)
      * cos(lat[iy] * M_PI / 180.);
  }

  /* Set time interval for output... */
  int year, month, day, hour, minute;
  double t;
  sscanf(argv[2], "atm_%d_%d_%d_%d_%d.tab", &year, &month, &day, &hour, &minute);
  time2jsec(year, month, day, hour, minute, 0, 0, &t);

  double t0 = t - 0.5 * ctl.dt_mod;
  double t1 = t + 0.5 * ctl.dt_mod;

  /* Get grid box indices... */
  for (int ip = 0; ip < atm->np; ip++) {
    ixs[ip] = (int) ((atm->lon[ip] - ctl.grid_lon0) / dlon);
    iys[ip] = (int) ((atm->lat[ip] - ctl.grid_lat0) / dlat);
    izs[ip] = (int) ((Z(atm->p[ip]) - ctl.grid_z0) / dz);
    if (atm->time[ip] < t0 || atm->time[ip] > t1
	|| ixs[ip] < 0 || ixs[ip] >= ctl.grid_nx
	|| iys[ip] < 0 || iys[ip] >= ctl.grid_ny
	|| izs[ip] < 0 || izs[ip] >= ctl.grid_nz)
      izs[ip] = -1;
  }

  /* Average data... */
  for (int ip = 0; ip < atm->np; ip++)
    if (izs[ip] >= 0) {
      int idx =
	ARRAY_3D(ixs[ip], iys[ip], ctl.grid_ny, izs[ip], ctl.grid_nz);
      double kernel = kernel_weight(kz, kw, nk, atm->p[ip]);
      np[idx]++;
      for (int iq = 0; iq < ctl.nq; iq++) {
	mean[iq][idx] += kernel * atm->q[iq][ip];
	stddev[iq][idx] += SQR(kernel * atm->q[iq][ip]);
      }
    }

  /* Calculate column density and volume mixing ratio... */
  for (int ix = 0; ix < ctl.grid_nx; ix++)
    for (int iy = 0; iy < ctl.grid_ny; iy++)
      for (int iz = 0; iz < ctl.grid_nz; iz++) {

	/* Get grid index... */
	int idx = ARRAY_3D(ix, iy, ctl.grid_ny, iz, ctl.grid_nz);

	/* Calculate column density... */
	cd[idx] = GSL_NAN;
	if (ctl.qnt_m >= 0)
	  cd[idx] = mean[ctl.qnt_m][idx] / (1e6 * area[iy]);

	/* Calculate volume mixing ratio (implicit)... */
	vmr_impl[idx] = GSL_NAN;
	if (ctl.qnt_m >= 0 && ctl.molmass > 0 && ctl.qnt_t >= 0) {
	  vmr_impl[idx] = 0;
	  if (mean[ctl.qnt_m][idx] > 0) {

	    /* Calculate volume mixing ratio... */
	    vmr_impl[idx] = MA / ctl.molmass * mean[ctl.qnt_m][idx]
	      / (RHO(press[iz], mean[ctl.qnt_t][idx]) * 1e6 * area[iy] * 1e3 * dz);
	  }
	}

	/* Calculate mean... */
	if (np[idx] > 0)
	  for (int iq = 0; iq < ctl.nq; iq++) {
	    mean[iq][idx] /= np[idx];
	    double var = stddev[iq][idx] / np[idx] - SQR(mean[iq][idx]);
	    stddev[iq][idx] = (var > 0 ? sqrt(var) : 0);
	} else
	  for (int iq = 0; iq < ctl.nq; iq++) {
	    mean[iq][idx] = GSL_NAN;
	    stddev[iq][idx] = GSL_NAN;
	  }
      }

  /* Write ASCII data... */
  if (ctl.grid_type == 0)
    write_grid_asc(argv[3], &ctl, cd, mean, stddev, vmr_impl,
		   t, z, lon, lat, area, dz, np);

  /* Write netCDF data... */
  else if (ctl.grid_type == 1)
    write_grid_nc(argv[3], &ctl, cd, mean, stddev, vmr_impl,
		  t, z, lon, lat, area, dz, np);

  /* Error message... */
  else
    ERRMSG("Grid data format GRID_TYPE unknown!");

  /* Free... */
  free(cd);
  for (int iq = 0; iq < ctl.nq; iq++) {
    free(mean[iq]);
    free(stddev[iq]);
  }
  free(vmr_impl);
  free(z);
  free(lon);
  free(lat);
  free(area);
  free(press);
  free(np);
  free(ixs);
  free(iys);
  free(izs);
  free(atm);

  return EXIT_SUCCESS;

}