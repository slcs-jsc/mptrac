/* Check the dimensions and representative finite values of synthetic wind files. */
#include <math.h>
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>

static void check(int status) {
  if (status != NC_NOERR) {
    fprintf(stderr, "%s\n", nc_strerror(status));
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char **argv) {
  if (argc < 5) {
    fprintf(stderr, "Usage: validate_wind_grid NX NY NZ FILE...\n");
    return EXIT_FAILURE;
  }
  const size_t expected[] = {
    (size_t) atoi(argv[1]), (size_t) atoi(argv[2]),
    (size_t) atoi(argv[3])
  };
  const char *dims[] = { "lon", "lat", "lev" };
  const char *fields[] = { "t", "u", "v", "w" };
  for (int file = 4; file < argc; ++file) {
    int ncid;
    check(nc_open(argv[file], NC_NOWRITE, &ncid));
    for (int dim = 0; dim < 3; ++dim) {
      int dimid;
      size_t length;
      check(nc_inq_dimid(ncid, dims[dim], &dimid));
      check(nc_inq_dimlen(ncid, dimid, &length));
      if (length != expected[dim]) {
        fprintf(stderr, "Wrong %s dimension in %s: %zu != %zu\n",
                dims[dim], argv[file], length, expected[dim]);
        return EXIT_FAILURE;
      }
    }
    for (int field = 0; field < 4; ++field) {
      int varid;
      check(nc_inq_varid(ncid, fields[field], &varid));
      for (int sample = 0; sample < 3; ++sample) {
        size_t index[] = { 0, sample * (expected[2] - 1) / 2,
                           sample * (expected[1] - 1) / 2,
                           sample * (expected[0] - 1) / 2 };
        float value;
        check(nc_get_var1_float(ncid, varid, index, &value));
        if (!isfinite(value)) {
          fprintf(stderr, "Non-finite %s sample in %s\n",
                  fields[field], argv[file]);
          return EXIT_FAILURE;
        }
      }
    }
    check(nc_close(ncid));
  }
  printf("Validated dimensions and finite wind samples in %d files\n", argc - 4);
  return EXIT_SUCCESS;
}
