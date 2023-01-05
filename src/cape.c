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
  
  Copyright (C) 2013-2022 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Add CAPE data to netCDF file.
*/

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  clim_t *clim;

  met_t *met;

  float help[EX * EY];

  int dims[10], ncid, varid;

  size_t start[10], count[10];

  /* Allocate... */
  ALLOC(clim, clim_t, 1);
  ALLOC(met, met_t, 1);

  /* Check arguments... */
  if (argc < 2)
    ERRMSG("Give parameters: <ctl> <met.nc>");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);

  /* Read climatological data... */
  read_clim(&ctl, clim);

  /* Read meteorological data... */
  if (!read_met(argv[2], &ctl, clim, met))
    ERRMSG("Cannot open file!");

  /* Open netCDF file... */
  if (nc_open(argv[2], NC_WRITE, &ncid) != NC_NOERR)
    ERRMSG("Cannot open file!");

  /* Get dimensions... */
  NC_INQ_DIM("time", &dims[0], 1, 1);
  NC_INQ_DIM("lat", &dims[1], met->ny, met->ny);
  NC_INQ_DIM("lon", &dims[2], met->nx - 1, met->nx - 1);
  NC(nc_inq_dimid(ncid, "time", &dims[0]));
  NC(nc_inq_dimid(ncid, "lat", &dims[1]));
  NC(nc_inq_dimid(ncid, "lon", &dims[2]));

  /* Set define mode... */
  NC(nc_redef(ncid));

  /* Create variables... */
  NC_DEF_VAR("CAPE_MPT", NC_FLOAT, 3, dims,
	     "convective available potential energy", "J kg**-1");
  NC_DEF_VAR("CIN_MPT", NC_FLOAT, 3, dims,
	     "convective inhibition", "J kg**-1");
  NC_DEF_VAR("PEL_MPT", NC_FLOAT, 3, dims,
	     "pressure at equilibrium level", "hPa");

  /* End define mode... */
  NC(nc_enddef(ncid));

  /* Write data... */
  for (int ix = 0; ix < met->nx - 1; ix++)
    for (int iy = 0; iy < met->ny; iy++)
      help[ARRAY_2D(iy, ix, met->nx - 1)] = met->cape[ix][iy];
  NC_PUT_FLOAT("CAPE_MPT", help, 0);

  for (int ix = 0; ix < met->nx - 1; ix++)
    for (int iy = 0; iy < met->ny; iy++)
      help[ARRAY_2D(iy, ix, met->nx - 1)] = met->cin[ix][iy];
  NC_PUT_FLOAT("CIN_MPT", help, 0);

  for (int ix = 0; ix < met->nx - 1; ix++)
    for (int iy = 0; iy < met->ny; iy++)
      help[ARRAY_2D(iy, ix, met->nx - 1)] = met->pel[ix][iy];
  NC_PUT_FLOAT("PEL_MPT", help, 0);

  /* Close file... */
  nc_close(ncid);

  /* Free... */
  free(clim);
  free(met);

  return EXIT_SUCCESS;
}
