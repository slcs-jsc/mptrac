/* Encode global x/y/level coordinates in a synthetic wind file for the
 * direct mode-0/mode-1 DD reader comparison. Modify a private copy only. */
#include <math.h>
#include <netcdf.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(call) do { \
    const int status = (call); \
    if (status != NC_NOERR) { \
      fprintf(stderr, "%s: %s\n", #call, nc_strerror(status)); \
      return EXIT_FAILURE; \
    } \
  } while (0)

int main(int argc, char **argv) {

  if (argc != 2 && argc != 3) {
    fprintf(stderr, "Usage: met_reader_encode private_wind.nc "
	    "[--duplicate-endpoint]\n");
    return EXIT_FAILURE;
  }
  const int endpoint = argc == 3;
  if (endpoint && strcmp(argv[2], "--duplicate-endpoint") != 0)
    return EXIT_FAILURE;
  int ncid, xdim, ydim, zdim, varid;
  size_t nx, ny, nz;
  CHECK(nc_open(argv[1], NC_WRITE, &ncid));
  CHECK(nc_inq_dimid(ncid, "lon", &xdim));
  CHECK(nc_inq_dimid(ncid, "lat", &ydim));
  CHECK(nc_inq_dimid(ncid, "lev", &zdim));
  CHECK(nc_inq_dimlen(ncid, xdim, &nx));
  CHECK(nc_inq_dimlen(ncid, ydim, &ny));
  CHECK(nc_inq_dimlen(ncid, zdim, &nz));
  if (nx > SIZE_MAX / ny || nx * ny > SIZE_MAX / sizeof(float))
    return EXIT_FAILURE;
  float *slice = malloc(nx * ny * sizeof(float));
  if (slice == NULL)
    return EXIT_FAILURE;

  if (endpoint) {
    double *lon = malloc(nx * sizeof(double));
    if (lon == NULL)
      return EXIT_FAILURE;
    CHECK(nc_inq_varid(ncid, "lon", &varid));
    for (size_t x = 0; x < nx; x++)
      lon[x] = 360.0 * (double) x / (double) (nx - 1);
    CHECK(nc_put_var_double(ncid, varid, lon));
    free(lon);
  }

  CHECK(nc_inq_varid(ncid, "u", &varid));
  for (size_t k = 0; k < nz; k++) {
    for (size_t y = 0; y < ny; y++)
      for (size_t x = 0; x < nx; x++)
	slice[y * nx + x] = (float) (10000 * k + 100 * y + x + 1);
    if (endpoint)
      for (size_t y = 0; y < ny; y++)
	slice[y * nx + nx - 1] = slice[y * nx];
    if (k == 1) {
      slice[nx + 1] = NAN;
      slice[2 * nx + 2] = 1e20f;
    }
    const size_t start[4] = { 0, k, 0, 0 };
    const size_t count[4] = { 1, 1, ny, nx };
    CHECK(nc_put_vara_float(ncid, varid, start, count, slice));
  }

  CHECK(nc_inq_varid(ncid, "sp", &varid));
  for (size_t y = 0; y < ny; y++)
    for (size_t x = 0; x < nx; x++)
      slice[y * nx + x] = 100000.f + (float) (100 * y + x);
  if (endpoint)
    for (size_t y = 0; y < ny; y++)
      slice[y * nx + nx - 1] = slice[y * nx];
  const size_t start[3] = { 0, 0, 0 };
  const size_t count[3] = { 1, ny, nx };
  CHECK(nc_put_vara_float(ncid, varid, start, count, slice));

  free(slice);
  CHECK(nc_close(ncid));
  return EXIT_SUCCESS;
}
