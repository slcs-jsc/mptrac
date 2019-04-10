
#include "libtrac.h"
#include "rng.h"

static gsl_rng *rng[NTHREADS];

void init_rng()
{
    /* Initialize random number generators... */
    gsl_rng_env_setup();
    if (omp_get_max_threads() > NTHREADS)
      ERRMSG("Too many threads!");
    for (int i = 0; i < NTHREADS; i++) {
      rng[i] = gsl_rng_alloc(gsl_rng_default);
      gsl_rng_set(rng[i], gsl_rng_default_seed + (long unsigned) i);
    }
}

void fill_rng_buffer(double* buffer, size_t n)
{
    for(size_t i = 0; i < n; ++i)
        buffer[i] = gsl_ran_gaussian_ziggurat(rng[0], 1.0);
}
