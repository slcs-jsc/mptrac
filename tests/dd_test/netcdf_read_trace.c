/* Optional LD_PRELOAD profiler for the DD NetCDF I/O benchmark.
 * It observes parallel meteorological opens, float hyperslab reads and closes.
 * No MPTRAC production code or NetCDF access pattern is changed. */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <mpi.h>
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <time.h>

static int active_ncid = -1;
static char active_name[512];
static double open_seconds, core_seconds, halo_seconds, empty_seconds,
  close_seconds;
static unsigned long long core_bytes, halo_bytes;
static unsigned long long core_calls, halo_calls, empty_calls;
static size_t core_shape[4], halo_shape[4];

static double now(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double) ts.tv_sec + 1e-9 * (double) ts.tv_nsec;
}

int nc_open_par(const char *path, int mode, MPI_Comm comm, MPI_Info info,
                int *ncidp) {
  static int (*real_open)(const char *, int, MPI_Comm, MPI_Info, int *);
  if (!real_open) {
    void *symbol = dlsym(RTLD_NEXT, "nc_open_par");
    memcpy(&real_open, &symbol, sizeof(real_open));
  }
  double start = now();
  int status = real_open(path, mode, comm, info, ncidp);
  double elapsed = now() - start;
  if (status == NC_NOERR && getenv("DD_IO_TRACE_DIR")) {
    active_ncid = *ncidp;
    snprintf(active_name, sizeof(active_name), "%s", path);
    open_seconds = elapsed;
    core_seconds = halo_seconds = empty_seconds = close_seconds = 0;
    core_bytes = halo_bytes = core_calls = halo_calls = empty_calls = 0;
    memset(core_shape, 0, sizeof(core_shape));
    memset(halo_shape, 0, sizeof(halo_shape));
  }
  return status;
}

int nc_get_vara_float(int ncid, int varid, const size_t *startp,
                      const size_t *countp, float *ip) {
  static int (*real_read)(int, int, const size_t *, const size_t *, float *);
  if (!real_read) {
    void *symbol = dlsym(RTLD_NEXT, "nc_get_vara_float");
    memcpy(&real_read, &symbol, sizeof(real_read));
  }
  if (ncid != active_ncid)
    return real_read(ncid, varid, startp, countp, ip);
  int ndims = 0;
  if (nc_inq_varndims(ncid, varid, &ndims) != NC_NOERR || ndims < 3 || ndims > 4)
    return real_read(ncid, varid, startp, countp, ip);
  double start = now();
  int status = real_read(ncid, varid, startp, countp, ip);
  double elapsed = now() - start;
  if (status == NC_NOERR) {
    unsigned long long elements = 1;
    for (int i = 0; i < ndims; ++i)
      elements *= (unsigned long long) countp[i];
    if (elements == 0) {
      empty_seconds += elapsed;
      ++empty_calls;
      return status;
    }
    int halo = countp[ndims - 1] <= 2;
    double *seconds = halo ? &halo_seconds : &core_seconds;
    unsigned long long *bytes = halo ? &halo_bytes : &core_bytes;
    unsigned long long *calls = halo ? &halo_calls : &core_calls;
    size_t *shape = halo ? halo_shape : core_shape;
    *seconds += elapsed;
    *bytes += elements * sizeof(float);
    ++*calls;
    if (ndims == 4)
      memcpy(shape, countp, 4 * sizeof(size_t));
  }
  return status;
}

int nc_close(int ncid) {
  static int (*real_close)(int);
  if (!real_close) {
    void *symbol = dlsym(RTLD_NEXT, "nc_close");
    memcpy(&real_close, &symbol, sizeof(real_close));
  }
  if (ncid != active_ncid)
    return real_close(ncid);
  double start = now();
  int status = real_close(ncid);
  close_seconds = now() - start;
  const char *directory = getenv("DD_IO_TRACE_DIR");
  const char *rank = getenv("SLURM_PROCID");
  if (!rank)
    rank = getenv("OMPI_COMM_WORLD_RANK");
  if (directory && rank) {
    char filename[1024];
    snprintf(filename, sizeof(filename), "%s/rank_%s.csv", directory, rank);
    FILE *stream = fopen(filename, "a");
    if (stream) {
      struct rusage usage;
      getrusage(RUSAGE_SELF, &usage);
      if (ftell(stream) == 0)
        fprintf(stream, "file,core_seconds,halo_seconds,close_seconds,core_calls,halo_calls,core_bytes,halo_bytes,core_lev,core_lat,core_lon,halo_lev,halo_lat,halo_lon,peak_rss_kib,open_seconds,empty_seconds,empty_calls\n");
      fprintf(stream, "%s,%.9f,%.9f,%.9f,%llu,%llu,%llu,%llu,%zu,%zu,%zu,%zu,%zu,%zu,%ld,%.9f,%.9f,%llu\n",
              active_name, core_seconds, halo_seconds, close_seconds,
              core_calls, halo_calls, core_bytes, halo_bytes,
              core_shape[1], core_shape[2], core_shape[3],
              halo_shape[1], halo_shape[2], halo_shape[3],
              usage.ru_maxrss, open_seconds, empty_seconds, empty_calls);
      fclose(stream);
    }
  }
  active_ncid = -1;
  return status;
}
