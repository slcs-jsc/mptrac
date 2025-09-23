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
  
  Copyright (C) 2013-2025 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Create meteorological data files with synthetic wind fields.
*/

#include "mptrac.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  static char filename[LEN];

  static double r, dataLon[EX], dataLat[EY], dataZ[EP];

  static float *dataT, *dataU, *dataV, *dataW;

  static int ncid, varid, dims[4], year, mon, day, hour, min, sec;

  static size_t start[4], count[4];

  /* Allocate... */
  ALLOC(dataT, float,
	EP * EY * EX);
  ALLOC(dataU, float,
	EP * EY * EX);
  ALLOC(dataV, float,
	EP * EY * EX);
  ALLOC(dataW, float,
	EP * EY * EX);

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <metbase>");

  /* Read control parameters... */
  mptrac_read_ctl(argv[1], argc, argv, &ctl);
  double t0 = scan_ctl(argv[1], argc, argv, "WIND_T0", -1, "0", NULL);
  const int nx =
    (int) scan_ctl(argv[1], argc, argv, "WIND_NX", -1, "360", NULL);
  const int ny =
    (int) scan_ctl(argv[1], argc, argv, "WIND_NY", -1, "181", NULL);
  const int nz =
    (int) scan_ctl(argv[1], argc, argv, "WIND_NZ", -1, "61", NULL);
  const double z0 = scan_ctl(argv[1], argc, argv, "WIND_Z0", -1, "0", NULL);
  const double z1 = scan_ctl(argv[1], argc, argv, "WIND_Z1", -1, "60", NULL);
  const double u0 =
    scan_ctl(argv[1], argc, argv, "WIND_U0", -1, "38.587660177302", NULL);
  const double u1 =
    scan_ctl(argv[1], argc, argv, "WIND_U1", -1, "38.587660177302", NULL);
  const double w0 = scan_ctl(argv[1], argc, argv, "WIND_W0", -1, "0", NULL);
  const double alpha =
    scan_ctl(argv[1], argc, argv, "WIND_ALPHA", -1, "0.0", NULL);
  const int lat_reverse =
    (int) scan_ctl(argv[1], argc, argv, "WIND_LAT_REVERSE", -1, "0", NULL);

  /* Check dimensions... */
  if (nx < 1 || nx > EX)
    ERRMSG("Set 1 <= NX <= MAX!");
  if (ny < 1 || ny > EY)
    ERRMSG("Set 1 <= NY <= MAX!");
  if (nz < 1 || nz > EP)
    ERRMSG("Set 1 <= NZ <= MAX!");

  /* Get time... */
  jsec2time(t0, &year, &mon, &day, &hour, &min, &sec, &r);
  t0 = year * 10000. + mon * 100. + day + hour / 24.;

  /* Set filename... */
  sprintf(filename, "%s_%d_%02d_%02d_%02d.nc", argv[2], year, mon, day, hour);

  /* Create netCDF file... */
  NC(nc_create(filename, NC_NETCDF4, &ncid));

  /* Create dimensions... */
  NC(nc_def_dim(ncid, "time", 1, &dims[0]));
  NC(nc_def_dim(ncid, "lev", (size_t) nz, &dims[1]));
  NC(nc_def_dim(ncid, "lat", (size_t) ny, &dims[2]));
  NC(nc_def_dim(ncid, "lon", (size_t) nx, &dims[3]));

  /* Create variables... */
  NC_DEF_VAR("time", NC_DOUBLE, 1, &dims[0], "time", "day as %Y%m%d.%f", 0,
	     0);
  NC_DEF_VAR("lev", NC_DOUBLE, 1, &dims[1], "air_pressure", "Pa", 0, 0);
  NC_DEF_VAR("lat", NC_DOUBLE, 1, &dims[2], "latitude", "degrees_north", 0,
	     0);
  NC_DEF_VAR("lon", NC_DOUBLE, 1, &dims[3], "longitude", "degrees_east", 0,
	     0);
  NC_DEF_VAR("T", NC_FLOAT, 4, &dims[0], "Temperature", "K", 0, 0);
  NC_DEF_VAR("U", NC_FLOAT, 4, &dims[0], "zonal wind", "m s**-1", 0, 0);
  NC_DEF_VAR("V", NC_FLOAT, 4, &dims[0], "meridional wind", "m s**-1", 0, 0);
  NC_DEF_VAR("W", NC_FLOAT, 4, &dims[0], "vertical velocity", "Pa s**-1", 0,
	     0);

  /* End definition... */
  NC(nc_enddef(ncid));

  /* Set coordinates... */
  for (int ix = 0; ix < nx; ix++)
    dataLon[ix] = 360.0 / nx * (double) ix;
  for (int iy = 0; iy < ny; iy++)
    dataLat[iy] = (lat_reverse ? - (180.0 / (ny - 1) * (double) iy - 90) : (180.0 / (ny - 1) * (double) iy - 90));
  for (int iz = 0; iz < nz; iz++)
    dataZ[iz] = 100. * P(LIN(0.0, z0, nz - 1.0, z1, iz));

  /* Write coordinates... */
  NC_PUT_DOUBLE("time", &t0, 0);
  NC_PUT_DOUBLE("lev", dataZ, 0);
  NC_PUT_DOUBLE("lat", dataLat, 0);
  NC_PUT_DOUBLE("lon", dataLon, 0);

  /* Create wind fields (Williamson et al., 1992)... */
  for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++)
      for (int iz = 0; iz < nz; iz++) {
	int idx = (iz * ny + iy) * nx + ix;
	dataU[idx] = (float) (LIN(0.0, u0, nz - 1.0, u1, iz)
			      * (cos(DEG2RAD(dataLat[iy]))
				 * cos(DEG2RAD(alpha))
				 + sin(DEG2RAD(dataLat[iy]))
				 * cos(DEG2RAD(dataLon[ix]))
				 * sin(DEG2RAD(alpha))));
	dataV[idx] = (float) (-LIN(0.0, u0, nz - 1.0, u1, iz)
			      * sin(DEG2RAD(dataLon[ix]))
			      * sin(DEG2RAD(alpha)));
	dataW[idx] = (float) DZ2DP(1e-3 * w0, dataZ[iz]);
      }

  /* Write data... */
  NC_PUT_FLOAT("T", dataT, 0);
  NC_PUT_FLOAT("U", dataU, 0);
  NC_PUT_FLOAT("V", dataV, 0);
  NC_PUT_FLOAT("W", dataW, 0);

  /* Close file... */
  NC(nc_close(ncid));

  /* Free... */
  free(dataT);
  free(dataU);
  free(dataV);
  free(dataW);

  return EXIT_SUCCESS;
}
