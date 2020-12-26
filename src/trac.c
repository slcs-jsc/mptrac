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
  
  Copyright (C) 2013-2020 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Lagrangian particle dispersion model.
*/

/* ------------------------------------------------------------
   Includes...
   ------------------------------------------------------------ */

#include "libtrac.h"

#ifdef MPI
#include "mpi.h"
#endif

#ifdef _OPENACC
#include "openacc.h"
#include "curand.h"
#endif

/* ------------------------------------------------------------
   NVIDIA Tools Extension (NVTX)...
   ------------------------------------------------------------ */

#ifdef USE_NVTX
#include "nvToolsExt.h"

/*! Red color code. */
#define RED 0xFFFF0000

/*! Blue color code. */
#define BLUE 0xFF0000FF

/*! Yellow color code. */
#define YELLOW 0xFFFFFF00

/*! Gray color code. */
#define GRAY 0xFF808080

/*! Macro for calling nvtxRangePushEx. */
#define RANGE_PUSH(range_title, range_color) {		\
    nvtxEventAttributes_t eventAttrib = {0};		\
    eventAttrib.version = NVTX_VERSION;			\
    eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;	\
    eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;	\
    eventAttrib.colorType = NVTX_COLOR_ARGB;		\
    eventAttrib.color = range_color;			\
    eventAttrib.message.ascii = range_title;		\
    nvtxRangePushEx(&eventAttrib);			\
  }

/*! Macro for calling nvtxRangePop. */
#define RANGE_POP {				\
    nvtxRangePop();				\
  }
#else

/* Empty definitions of RANGE_PUSH and RANGE_POP... */
#define RANGE_PUSH(range_title, range_color) {}
#define RANGE_POP {}
#endif

/* ------------------------------------------------------------
   Global variables...
   ------------------------------------------------------------ */

#ifdef _OPENACC
curandGenerator_t rng;
#else
static gsl_rng *rng[NTHREADS];
#endif

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Calculate advection of air parcels. */
void module_advection(
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/*! Calculate exponential decay of particle mass. */
void module_decay(
  ctl_t * ctl,
  atm_t * atm,
  double *dt);

/*! Initialize random number generator... */
void module_diffusion_init(
  void);

/*! Calculate mesoscale diffusion. */
void module_diffusion_meso(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  cache_t * cache,
  double *dt,
  double *rs);

/*! Generate random numbers. */
void module_diffusion_rng(
  double *rs,
  size_t n);

/*! Calculate turbulent diffusion. */
void module_diffusion_turb(
  ctl_t * ctl,
  atm_t * atm,
  double *dt,
  double *rs);

/*! Initialize isosurface module. */
void module_isosurf_init(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  cache_t * cache);

/*! Force air parcels to stay on isosurface. */
void module_isosurf(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  cache_t * cache);

/*! Interpolate meteorological data for air parcel positions. */
void module_meteo(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm);

/*! Check position of air parcels. */
void module_position(
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/*! Calculate sedimentation of air parcels. */
void module_sedi(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/*! Calculate OH chemistry. */
void module_oh_chem(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/*! Calculate wet deposition. */
void module_wet_deposition(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt);

/*! Write simulation output. */
void write_output(
  const char *dirname,
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double t);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  atm_t *atm;

  cache_t *cache;

  met_t *met0, *met1;

  FILE *dirlist;

  char dirname[LEN], filename[2 * LEN];

  double *dt, *rs, t;

  int ntask = -1, rank = 0, size = 1;

  /* Initialize MPI... */
#ifdef MPI
  RANGE_PUSH("Initialize MPI", GRAY);
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  RANGE_POP;
#endif

  /* Initialize GPUs... */
#ifdef _OPENACC
  RANGE_PUSH("Initialize GPUs", GRAY);
  acc_device_t device_type = acc_get_device_type();
  int num_devices = acc_get_num_devices(acc_device_nvidia);
  int device_num = rank % num_devices;
  acc_set_device_num(device_num, acc_device_nvidia);
  acc_init(device_type);
  RANGE_POP;
#endif

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <dirlist> <ctl> <atm_in> <metbase>");

  /* Open directory list... */
  if (!(dirlist = fopen(argv[1], "r")))
    ERRMSG("Cannot open directory list!");

  /* Loop over directories... */
  while (fscanf(dirlist, "%s", dirname) != EOF) {

    /* MPI parallelization... */
    if ((++ntask) % size != rank)
      continue;

    /* ------------------------------------------------------------
       Initialize model run...
       ------------------------------------------------------------ */

    /* Set timers... */
    START_TIMER(TIMER_ZERO);
    START_TIMER(TIMER_TOTAL);
    START_TIMER(TIMER_INIT);

    /* Allocate... */
    RANGE_PUSH("Allocate", GRAY);
    ALLOC(atm, atm_t, 1);
    ALLOC(cache, cache_t, 1);
    ALLOC(met0, met_t, 1);
    ALLOC(met1, met_t, 1);
    ALLOC(dt, double,
	  NP);
    ALLOC(rs, double,
	  3 * NP);
    RANGE_POP;

    /* Read control parameters... */
    RANGE_PUSH("Read ctl", GRAY);
    sprintf(filename, "%s/%s", dirname, argv[2]);
    read_ctl(filename, argc, argv, &ctl);
    RANGE_POP;

    /* Read atmospheric data... */
    RANGE_PUSH("Read atm", GRAY);
    sprintf(filename, "%s/%s", dirname, argv[3]);
    if (!read_atm(filename, &ctl, atm))
      ERRMSG("Cannot open file!");
    RANGE_POP;

    /* Set start time... */
    if (ctl.direction == 1) {
      ctl.t_start = gsl_stats_min(atm->time, 1, (size_t) atm->np);
      if (ctl.t_stop > 1e99)
	ctl.t_stop = gsl_stats_max(atm->time, 1, (size_t) atm->np);
    } else {
      ctl.t_start = gsl_stats_max(atm->time, 1, (size_t) atm->np);
      if (ctl.t_stop > 1e99)
	ctl.t_stop = gsl_stats_min(atm->time, 1, (size_t) atm->np);
    }

    /* Check time interval... */
    if (ctl.direction * (ctl.t_stop - ctl.t_start) <= 0)
      ERRMSG("Nothing to do!");

    /* Round start time... */
    if (ctl.direction == 1)
      ctl.t_start = floor(ctl.t_start / ctl.dt_mod) * ctl.dt_mod;
    else
      ctl.t_start = ceil(ctl.t_start / ctl.dt_mod) * ctl.dt_mod;

    /* Copy to GPU... */
    RANGE_PUSH("Copy to GPU", GRAY);
#ifdef _OPENACC
#pragma acc enter data copyin(ctl)
#pragma acc enter data create(atm[:1],cache[:1],met0[:1],met1[:1],dt[:NP],rs[:3*NP])
#pragma acc update device(atm[:1],cache[:1],ctl)
#endif
    RANGE_POP;

    /* Initialize random number generator... */
    module_diffusion_init();

    /* Set timers... */
    STOP_TIMER(TIMER_INIT);

    /* Initialize meteorological data... */
    START_TIMER(TIMER_INPUT);
    get_met(&ctl, argv[4], ctl.t_start, &met0, &met1);
    if (ctl.dt_mod > fabs(met0->lon[1] - met0->lon[0]) * 111132. / 150.)
      WARN("Violation of CFL criterion! Check DT_MOD!");
    STOP_TIMER(TIMER_INPUT);

    /* Initialize isosurface... */
    START_TIMER(TIMER_ISOSURF);
    if (ctl.isosurf >= 1 && ctl.isosurf <= 4)
      module_isosurf_init(&ctl, met0, met1, atm, cache);
    STOP_TIMER(TIMER_ISOSURF);

    /* ------------------------------------------------------------
       Loop over timesteps...
       ------------------------------------------------------------ */

    /* Loop over timesteps... */
    for (t = ctl.t_start; ctl.direction * (t - ctl.t_stop) < ctl.dt_mod;
	 t += ctl.direction * ctl.dt_mod) {

      /* Adjust length of final time step... */
      if (ctl.direction * (t - ctl.t_stop) > 0)
	t = ctl.t_stop;

      /* Set time steps for air parcels... */
      RANGE_PUSH("Set time steps", GRAY);
#ifdef _OPENACC
#pragma acc parallel loop independent gang vector present(ctl,atm,atm->time,dt)
#endif
      for (int ip = 0; ip < atm->np; ip++) {
	double atmtime = atm->time[ip];
	double tstart = ctl.t_start;
	double tstop = ctl.t_stop;
	int dir = ctl.direction;
	if ((dir * (atmtime - tstart) >= 0 && dir * (atmtime - tstop) <= 0
	     && dir * (atmtime - t) < 0))
	  dt[ip] = t - atmtime;
	else
	  dt[ip] = 0;
      }
      RANGE_POP;

      /* Get meteorological data... */
      RANGE_PUSH("Get met data", GRAY);
      START_TIMER(TIMER_INPUT);
      if (t != ctl.t_start)
	get_met(&ctl, argv[4], t, &met0, &met1);
      STOP_TIMER(TIMER_INPUT);
      RANGE_POP;

      /* Check initial position... */
      RANGE_PUSH("Check init pos", GRAY);
      START_TIMER(TIMER_POSITION);
      module_position(met0, met1, atm, dt);
      STOP_TIMER(TIMER_POSITION);
      RANGE_POP;

      /* Advection... */
      RANGE_PUSH("Advection", GRAY);
      START_TIMER(TIMER_ADVECT);
      module_advection(met0, met1, atm, dt);
      STOP_TIMER(TIMER_ADVECT);
      RANGE_POP;

      /* Turbulent diffusion... */
      RANGE_PUSH("Turbulent diffusion", GRAY);
      START_TIMER(TIMER_DIFFTURB);
      if (ctl.turb_dx_trop > 0 || ctl.turb_dz_trop > 0
	  || ctl.turb_dx_strat > 0 || ctl.turb_dz_strat > 0) {
	module_diffusion_rng(rs, 3 * (size_t) atm->np);
	module_diffusion_turb(&ctl, atm, dt, rs);
      }
      STOP_TIMER(TIMER_DIFFTURB);
      RANGE_POP;

      /* Mesoscale diffusion... */
      RANGE_PUSH("Mesoscale diffusion", GRAY);
      START_TIMER(TIMER_DIFFMESO);
      if (ctl.turb_mesox > 0 || ctl.turb_mesoz > 0) {
	module_diffusion_rng(rs, 3 * (size_t) atm->np);
	module_diffusion_meso(&ctl, met0, met1, atm, cache, dt, rs);
      }
      STOP_TIMER(TIMER_DIFFMESO);
      RANGE_POP;

      /* Sedimentation... */
      RANGE_PUSH("Sedimentation", GRAY);
      START_TIMER(TIMER_SEDI);
      if (ctl.qnt_r >= 0 && ctl.qnt_rho >= 0)
	module_sedi(&ctl, met0, met1, atm, dt);
      STOP_TIMER(TIMER_SEDI);
      RANGE_POP;

      /* Isosurface... */
      RANGE_PUSH("Isosurface", GRAY);
      START_TIMER(TIMER_ISOSURF);
      if (ctl.isosurf >= 1 && ctl.isosurf <= 4)
	module_isosurf(&ctl, met0, met1, atm, cache);
      STOP_TIMER(TIMER_ISOSURF);
      RANGE_POP;

      /* Check final position... */
      RANGE_PUSH("Check final pos", GRAY);
      START_TIMER(TIMER_POSITION);
      module_position(met0, met1, atm, dt);
      STOP_TIMER(TIMER_POSITION);
      RANGE_POP;

      /* Interpolate meteorological data... */
      RANGE_PUSH("Interpolate met data", GRAY);
      START_TIMER(TIMER_METEO);
      if (ctl.met_dt_out > 0
	  && (ctl.met_dt_out < ctl.dt_mod || fmod(t, ctl.met_dt_out) == 0))
	module_meteo(&ctl, met0, met1, atm);
      STOP_TIMER(TIMER_METEO);
      RANGE_POP;

      /* Decay of particle mass... */
      RANGE_PUSH("Decay of particle mass", GRAY);
      START_TIMER(TIMER_DECAY);
      if (ctl.tdec_trop > 0 && ctl.tdec_strat > 0)
	module_decay(&ctl, atm, dt);
      STOP_TIMER(TIMER_DECAY);
      RANGE_POP;

      /* OH chemistry... */
      RANGE_PUSH("OH chem", GRAY);
      START_TIMER(TIMER_OHCHEM);
      if (ctl.oh_chem[0] > 0 && ctl.oh_chem[2] > 0)
	module_oh_chem(&ctl, met0, met1, atm, dt);
      STOP_TIMER(TIMER_OHCHEM);
      RANGE_POP;

      /* Wet deposition... */
      RANGE_PUSH("Wet deposition", GRAY);
      START_TIMER(TIMER_WETDEPO);
      if (ctl.wet_depo[0] > 0 && ctl.wet_depo[1] > 0
	  && ctl.wet_depo[2] > 0 && ctl.wet_depo[3] > 0)
	module_wet_deposition(&ctl, met0, met1, atm, dt);
      STOP_TIMER(TIMER_WETDEPO);
      RANGE_POP;

      /* Write output... */
      RANGE_PUSH("Write output", GRAY);
      START_TIMER(TIMER_OUTPUT);
      write_output(dirname, &ctl, met0, met1, atm, t);
      STOP_TIMER(TIMER_OUTPUT);
      RANGE_POP;
    }

    /* ------------------------------------------------------------
       Finalize model run...
       ------------------------------------------------------------ */

    /* Report problem size... */
    printf("SIZE_NP = %d\n", atm->np);
    printf("SIZE_TASKS = %d\n", size);
    printf("SIZE_THREADS = %d\n", omp_get_max_threads());

    /* Report memory usage... */
    printf("MEMORY_ATM = %g MByte\n", sizeof(atm_t) / 1024. / 1024.);
    printf("MEMORY_CACHE = %g MByte\n", sizeof(cache_t) / 1024. / 1024.);
    printf("MEMORY_METEO = %g MByte\n", 2 * sizeof(met_t) / 1024. / 1024.);
    printf("MEMORY_DYNAMIC = %g MByte\n", (sizeof(met_t)
					   + 4 * NP * sizeof(double)
					   + EX * EY * EP * sizeof(float)) /
	   1024. / 1024.);
    printf("MEMORY_STATIC = %g MByte\n", (EX * EY * sizeof(double)
					  + EX * EY * EP * sizeof(float)
					  + 4 * GX * GY * GZ * sizeof(double)
					  + 2 * GX * GY * GZ * sizeof(int)
					  + 2 * GX * GY * sizeof(double)
					  + GX * GY * sizeof(int)) / 1024. /
	   1024.);

    /* Report timers... */
    STOP_TIMER(TIMER_ZERO);
    PRINT_TIMER(TIMER_INIT);
    PRINT_TIMER(TIMER_INPUT);
    PRINT_TIMER(TIMER_OUTPUT);
    PRINT_TIMER(TIMER_ADVECT);
    PRINT_TIMER(TIMER_DECAY);
    PRINT_TIMER(TIMER_DIFFMESO);
    PRINT_TIMER(TIMER_DIFFTURB);
    PRINT_TIMER(TIMER_ISOSURF);
    PRINT_TIMER(TIMER_METEO);
    PRINT_TIMER(TIMER_POSITION);
    PRINT_TIMER(TIMER_SEDI);
    PRINT_TIMER(TIMER_OHCHEM);
    PRINT_TIMER(TIMER_WETDEPO);
    STOP_TIMER(TIMER_TOTAL);
    PRINT_TIMER(TIMER_TOTAL);

    /* Free... */
    RANGE_PUSH("Deallocations", GRAY);
    free(atm);
    free(cache);
    free(met0);
    free(met1);
    free(dt);
    free(rs);
#ifdef _OPENACC
#pragma acc exit data delete(ctl,atm,cache,met0,met1,dt,rs)
#endif
    RANGE_POP;
  }

#ifdef MPI
  /* Finalize MPI... */
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void module_advection(
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt) {

#ifdef _OPENACC
#pragma acc data present(met0,met1,atm,dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++)
    if (dt[ip] != 0) {

      int ci[3] = { 0 };

      double dtm = 0.0, v[3] = { 0.0 }, xm[3] = {
      0.0};
      double cw[3] = { 0.0 };

      /* Interpolate meteorological data... */
      intpol_met_time_3d(met0, met0->u, met1, met1->u, atm->time[ip],
			 atm->p[ip], atm->lon[ip], atm->lat[ip], &v[0], ci,
			 cw, 1);
      intpol_met_time_3d(met0, met0->v, met1, met1->v, atm->time[ip],
			 atm->p[ip], atm->lon[ip], atm->lat[ip], &v[1], ci,
			 cw, 0);
      intpol_met_time_3d(met0, met0->w, met1, met1->w, atm->time[ip],
			 atm->p[ip], atm->lon[ip], atm->lat[ip], &v[2], ci,
			 cw, 0);

      /* Get position of the mid point... */
      dtm = atm->time[ip] + 0.5 * dt[ip];
      xm[0] =
	atm->lon[ip] + DX2DEG(0.5 * dt[ip] * v[0] / 1000., atm->lat[ip]);
      xm[1] = atm->lat[ip] + DY2DEG(0.5 * dt[ip] * v[1] / 1000.);
      xm[2] = atm->p[ip] + 0.5 * dt[ip] * v[2];

      /* Interpolate meteorological data for mid point... */
      intpol_met_time_3d(met0, met0->u, met1, met1->u, dtm, xm[2], xm[0],
			 xm[1], &v[0], ci, cw, 1);
      intpol_met_time_3d(met0, met0->v, met1, met1->v, dtm, xm[2], xm[0],
			 xm[1], &v[1], ci, cw, 0);
      intpol_met_time_3d(met0, met0->w, met1, met1->w, dtm, xm[2], xm[0],
			 xm[1], &v[2], ci, cw, 0);

      /* Save new position... */
      atm->time[ip] += dt[ip];
      atm->lon[ip] += DX2DEG(dt[ip] * v[0] / 1000., xm[1]);
      atm->lat[ip] += DY2DEG(dt[ip] * v[1] / 1000.);
      atm->p[ip] += dt[ip] * v[2];
    }
}

/*****************************************************************************/

void module_decay(
  ctl_t * ctl,
  atm_t * atm,
  double *dt) {

  /* Check quantity flags... */
  if (ctl->qnt_m < 0)
    ERRMSG("Module needs quantity mass!");

#ifdef _OPENACC
#pragma acc data present(ctl,atm,dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++)
    if (dt[ip] != 0) {

      double p0, p1, pt, tdec, w;

      /* Get tropopause pressure... */
      pt = clim_tropo(atm->time[ip], atm->lat[ip]);

      /* Get weighting factor... */
      p1 = pt * 0.866877899;
      p0 = pt / 0.866877899;
      if (atm->p[ip] > p0)
	w = 1;
      else if (atm->p[ip] < p1)
	w = 0;
      else
	w = LIN(p0, 1.0, p1, 0.0, atm->p[ip]);

      /* Set lifetime... */
      tdec = w * ctl->tdec_trop + (1 - w) * ctl->tdec_strat;

      /* Calculate exponential decay... */
      atm->q[ctl->qnt_m][ip] *= exp(-dt[ip] / tdec);
    }
}

/*****************************************************************************/

void module_diffusion_init(
  void) {

  /* Initialize random number generator... */
#ifdef _OPENACC

  if (curandCreateGenerator(&rng, CURAND_RNG_PSEUDO_DEFAULT)
      != CURAND_STATUS_SUCCESS)
    ERRMSG("Cannot create random number generator!");
  if (curandSetStream(rng, (cudaStream_t) acc_get_cuda_stream(acc_async_sync))
      != CURAND_STATUS_SUCCESS)
    ERRMSG("Cannot set stream for random number generator!");

#else

  gsl_rng_env_setup();
  if (omp_get_max_threads() > NTHREADS)
    ERRMSG("Too many threads!");
  for (int i = 0; i < NTHREADS; i++) {
    rng[i] = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng[i], gsl_rng_default_seed + (long unsigned) i);
  }

#endif
}

/*****************************************************************************/

void module_diffusion_meso(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  cache_t * cache,
  double *dt,
  double *rs) {

#ifdef _OPENACC
#pragma acc data present(ctl,met0,met1,atm,cache,dt,rs)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++)
    if (dt[ip] != 0) {

      double u[16], v[16], w[16];

      /* Get indices... */
      int ix = locate_reg(met0->lon, met0->nx, atm->lon[ip]);
      int iy = locate_reg(met0->lat, met0->ny, atm->lat[ip]);
      int iz = locate_irr(met0->p, met0->np, atm->p[ip]);

      /* Caching of wind standard deviations... */
      if (cache->tsig[ix][iy][iz] != met0->time) {

	/* Collect local wind data... */
	u[0] = met0->u[ix][iy][iz];
	u[1] = met0->u[ix + 1][iy][iz];
	u[2] = met0->u[ix][iy + 1][iz];
	u[3] = met0->u[ix + 1][iy + 1][iz];
	u[4] = met0->u[ix][iy][iz + 1];
	u[5] = met0->u[ix + 1][iy][iz + 1];
	u[6] = met0->u[ix][iy + 1][iz + 1];
	u[7] = met0->u[ix + 1][iy + 1][iz + 1];

	v[0] = met0->v[ix][iy][iz];
	v[1] = met0->v[ix + 1][iy][iz];
	v[2] = met0->v[ix][iy + 1][iz];
	v[3] = met0->v[ix + 1][iy + 1][iz];
	v[4] = met0->v[ix][iy][iz + 1];
	v[5] = met0->v[ix + 1][iy][iz + 1];
	v[6] = met0->v[ix][iy + 1][iz + 1];
	v[7] = met0->v[ix + 1][iy + 1][iz + 1];

	w[0] = met0->w[ix][iy][iz];
	w[1] = met0->w[ix + 1][iy][iz];
	w[2] = met0->w[ix][iy + 1][iz];
	w[3] = met0->w[ix + 1][iy + 1][iz];
	w[4] = met0->w[ix][iy][iz + 1];
	w[5] = met0->w[ix + 1][iy][iz + 1];
	w[6] = met0->w[ix][iy + 1][iz + 1];
	w[7] = met0->w[ix + 1][iy + 1][iz + 1];

	/* Collect local wind data... */
	u[8] = met1->u[ix][iy][iz];
	u[9] = met1->u[ix + 1][iy][iz];
	u[10] = met1->u[ix][iy + 1][iz];
	u[11] = met1->u[ix + 1][iy + 1][iz];
	u[12] = met1->u[ix][iy][iz + 1];
	u[13] = met1->u[ix + 1][iy][iz + 1];
	u[14] = met1->u[ix][iy + 1][iz + 1];
	u[15] = met1->u[ix + 1][iy + 1][iz + 1];

	v[8] = met1->v[ix][iy][iz];
	v[9] = met1->v[ix + 1][iy][iz];
	v[10] = met1->v[ix][iy + 1][iz];
	v[11] = met1->v[ix + 1][iy + 1][iz];
	v[12] = met1->v[ix][iy][iz + 1];
	v[13] = met1->v[ix + 1][iy][iz + 1];
	v[14] = met1->v[ix][iy + 1][iz + 1];
	v[15] = met1->v[ix + 1][iy + 1][iz + 1];

	w[8] = met1->w[ix][iy][iz];
	w[9] = met1->w[ix + 1][iy][iz];
	w[10] = met1->w[ix][iy + 1][iz];
	w[11] = met1->w[ix + 1][iy + 1][iz];
	w[12] = met1->w[ix][iy][iz + 1];
	w[13] = met1->w[ix + 1][iy][iz + 1];
	w[14] = met1->w[ix][iy + 1][iz + 1];
	w[15] = met1->w[ix + 1][iy + 1][iz + 1];

	/* Get standard deviations of local wind data... */
	cache->usig[ix][iy][iz] = (float) stddev(u, 16);
	cache->vsig[ix][iy][iz] = (float) stddev(v, 16);
	cache->wsig[ix][iy][iz] = (float) stddev(w, 16);
	cache->tsig[ix][iy][iz] = met0->time;
      }

      /* Set temporal correlations for mesoscale fluctuations... */
      double r = 1 - 2 * fabs(dt[ip]) / ctl->dt_met;
      double r2 = sqrt(1 - r * r);

      /* Calculate horizontal mesoscale wind fluctuations... */
      if (ctl->turb_mesox > 0) {
	cache->up[ip] = (float)
	  (r * cache->up[ip]
	   + r2 * rs[3 * ip] * ctl->turb_mesox * cache->usig[ix][iy][iz]);
	atm->lon[ip] += DX2DEG(cache->up[ip] * dt[ip] / 1000., atm->lat[ip]);

	cache->vp[ip] = (float)
	  (r * cache->vp[ip]
	   + r2 * rs[3 * ip + 1] * ctl->turb_mesox * cache->vsig[ix][iy][iz]);
	atm->lat[ip] += DY2DEG(cache->vp[ip] * dt[ip] / 1000.);
      }

      /* Calculate vertical mesoscale wind fluctuations... */
      if (ctl->turb_mesoz > 0) {
	cache->wp[ip] = (float)
	  (r * cache->wp[ip]
	   + r2 * rs[3 * ip + 2] * ctl->turb_mesoz * cache->wsig[ix][iy][iz]);
	atm->p[ip] += cache->wp[ip] * dt[ip];
      }
    }
}

/*****************************************************************************/

void module_diffusion_rng(
  double *rs,
  size_t n) {

#ifdef _OPENACC

#pragma acc host_data use_device(rs)
  {
    if (curandGenerateNormalDouble(rng, rs, n, 0.0, 1.0)
	!= CURAND_STATUS_SUCCESS)
      ERRMSG("Cannot create random numbers!");
  }

#else

#pragma omp parallel for default(shared)
  for (size_t i = 0; i < n; ++i)
    rs[i] = gsl_ran_gaussian_ziggurat(rng[omp_get_thread_num()], 1.0);

#endif

}

/*****************************************************************************/

void module_diffusion_turb(
  ctl_t * ctl,
  atm_t * atm,
  double *dt,
  double *rs) {

#ifdef _OPENACC
#pragma acc data present(ctl,atm,dt,rs)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++)
    if (dt[ip] != 0) {

      double w;

      /* Get tropopause pressure... */
      double pt = clim_tropo(atm->time[ip], atm->lat[ip]);

      /* Get weighting factor... */
      double p1 = pt * 0.866877899;
      double p0 = pt / 0.866877899;
      if (atm->p[ip] > p0)
	w = 1;
      else if (atm->p[ip] < p1)
	w = 0;
      else
	w = LIN(p0, 1.0, p1, 0.0, atm->p[ip]);

      /* Set diffusivity... */
      double dx = w * ctl->turb_dx_trop + (1 - w) * ctl->turb_dx_strat;
      double dz = w * ctl->turb_dz_trop + (1 - w) * ctl->turb_dz_strat;

      /* Horizontal turbulent diffusion... */
      if (dx > 0) {
	double sigma = sqrt(2.0 * dx * fabs(dt[ip]));
	atm->lon[ip] += DX2DEG(rs[3 * ip] * sigma / 1000., atm->lat[ip]);
	atm->lat[ip] += DY2DEG(rs[3 * ip + 1] * sigma / 1000.);
      }

      /* Vertical turbulent diffusion... */
      if (dz > 0) {
	double sigma = sqrt(2.0 * dz * fabs(dt[ip]));
	atm->p[ip]
	  += DZ2DP(rs[3 * ip + 2] * sigma / 1000., atm->p[ip]);
      }
    }
}

/*****************************************************************************/

void module_isosurf_init(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  cache_t * cache) {

  FILE *in;

  char line[LEN];

  double t, cw[3];

  int ci[3];

  /* Save pressure... */
  if (ctl->isosurf == 1)
    for (int ip = 0; ip < atm->np; ip++)
      cache->iso_var[ip] = atm->p[ip];

  /* Save density... */
  else if (ctl->isosurf == 2)
    for (int ip = 0; ip < atm->np; ip++) {
      intpol_met_time_3d(met0, met0->t, met1, met1->t, atm->time[ip],
			 atm->p[ip], atm->lon[ip], atm->lat[ip], &t, ci, cw,
			 1);
      cache->iso_var[ip] = atm->p[ip] / t;
    }

  /* Save potential temperature... */
  else if (ctl->isosurf == 3)
    for (int ip = 0; ip < atm->np; ip++) {
      intpol_met_time_3d(met0, met0->t, met1, met1->t, atm->time[ip],
			 atm->p[ip], atm->lon[ip], atm->lat[ip], &t, ci, cw,
			 1);
      cache->iso_var[ip] = THETA(atm->p[ip], t);
    }

  /* Read balloon pressure data... */
  else if (ctl->isosurf == 4) {

    /* Write info... */
    printf("Read balloon pressure data: %s\n", ctl->balloon);

    /* Open file... */
    if (!(in = fopen(ctl->balloon, "r")))
      ERRMSG("Cannot open file!");

    /* Read pressure time series... */
    while (fgets(line, LEN, in))
      if (sscanf(line, "%lg %lg", &(cache->iso_ts[cache->iso_n]),
		 &(cache->iso_ps[cache->iso_n])) == 2)
	if ((++cache->iso_n) > NP)
	  ERRMSG("Too many data points!");

    /* Check number of points... */
    if (cache->iso_n < 1)
      ERRMSG("Could not read any data!");

    /* Close file... */
    fclose(in);
  }
}

/*****************************************************************************/

void module_isosurf(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  cache_t * cache) {

#ifdef _OPENACC
#pragma acc data present(ctl,met0,met1,atm,cache)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++) {

    double t, cw[3];

    int ci[3];

    /* Restore pressure... */
    if (ctl->isosurf == 1)
      atm->p[ip] = cache->iso_var[ip];

    /* Restore density... */
    else if (ctl->isosurf == 2) {
      intpol_met_time_3d(met0, met0->t, met1, met1->t, atm->time[ip],
			 atm->p[ip], atm->lon[ip], atm->lat[ip], &t, ci, cw,
			 1);
      atm->p[ip] = cache->iso_var[ip] * t;
    }

    /* Restore potential temperature... */
    else if (ctl->isosurf == 3) {
      intpol_met_time_3d(met0, met0->t, met1, met1->t, atm->time[ip],
			 atm->p[ip], atm->lon[ip], atm->lat[ip], &t, ci, cw,
			 1);
      atm->p[ip] = 1000. * pow(cache->iso_var[ip] / t, -1. / 0.286);
    }

    /* Interpolate pressure... */
    else if (ctl->isosurf == 4) {
      if (atm->time[ip] <= cache->iso_ts[0])
	atm->p[ip] = cache->iso_ps[0];
      else if (atm->time[ip] >= cache->iso_ts[cache->iso_n - 1])
	atm->p[ip] = cache->iso_ps[cache->iso_n - 1];
      else {
	int idx = locate_irr(cache->iso_ts, cache->iso_n, atm->time[ip]);
	atm->p[ip] = LIN(cache->iso_ts[idx], cache->iso_ps[idx],
			 cache->iso_ts[idx + 1], cache->iso_ps[idx + 1],
			 atm->time[ip]);
      }
    }
  }
}

/*****************************************************************************/

void module_meteo(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm) {

  /* Check quantity flags... */
  if (ctl->qnt_tsts >= 0)
    if (ctl->qnt_tice < 0 || ctl->qnt_tnat < 0)
      ERRMSG("Need T_ice and T_NAT to calculate T_STS!");

#ifdef _OPENACC
#pragma acc data present(ctl,met0,met1,atm)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++) {

    double ps, pt, pc, pv, t, u, v, w, h2o, o3, lwc, iwc, z, cw[3];

    int ci[3];

    /* Interpolate meteorological data... */
    intpol_met_time_3d(met0, met0->z, met1, met1->z, atm->time[ip],
		       atm->p[ip], atm->lon[ip], atm->lat[ip], &z, ci, cw, 1);
    intpol_met_time_3d(met0, met0->t, met1, met1->t, atm->time[ip],
		       atm->p[ip], atm->lon[ip], atm->lat[ip], &t, ci, cw, 0);
    intpol_met_time_3d(met0, met0->u, met1, met1->u, atm->time[ip],
		       atm->p[ip], atm->lon[ip], atm->lat[ip], &u, ci, cw, 0);
    intpol_met_time_3d(met0, met0->v, met1, met1->v, atm->time[ip],
		       atm->p[ip], atm->lon[ip], atm->lat[ip], &v, ci, cw, 0);
    intpol_met_time_3d(met0, met0->w, met1, met1->w, atm->time[ip],
		       atm->p[ip], atm->lon[ip], atm->lat[ip], &w, ci, cw, 0);
    intpol_met_time_3d(met0, met0->pv, met1, met1->pv, atm->time[ip],
		       atm->p[ip], atm->lon[ip], atm->lat[ip], &pv, ci, cw,
		       0);
    intpol_met_time_3d(met0, met0->h2o, met1, met1->h2o, atm->time[ip],
		       atm->p[ip], atm->lon[ip], atm->lat[ip], &h2o, ci, cw,
		       0);
    intpol_met_time_3d(met0, met0->o3, met1, met1->o3, atm->time[ip],
		       atm->p[ip], atm->lon[ip], atm->lat[ip], &o3, ci, cw,
		       0);
    intpol_met_time_3d(met0, met0->lwc, met1, met1->lwc, atm->time[ip],
		       atm->p[ip], atm->lon[ip], atm->lat[ip], &lwc, ci, cw,
		       0);
    intpol_met_time_3d(met0, met0->iwc, met1, met1->iwc, atm->time[ip],
		       atm->p[ip], atm->lon[ip], atm->lat[ip], &iwc, ci, cw,
		       0);
    intpol_met_time_2d(met0, met0->ps, met1, met1->ps, atm->time[ip],
		       atm->lon[ip], atm->lat[ip], &ps, ci, cw, 0);
    intpol_met_time_2d(met0, met0->pt, met1, met1->pt, atm->time[ip],
		       atm->lon[ip], atm->lat[ip], &pt, ci, cw, 0);
    intpol_met_time_2d(met0, met0->pc, met1, met1->pc, atm->time[ip],
		       atm->lon[ip], atm->lat[ip], &pc, ci, cw, 0);

    /* Set surface pressure... */
    if (ctl->qnt_ps >= 0)
      atm->q[ctl->qnt_ps][ip] = ps;

    /* Set tropopause pressure... */
    if (ctl->qnt_pt >= 0)
      atm->q[ctl->qnt_pt][ip] = pt;

    /* Set pressure... */
    if (ctl->qnt_p >= 0)
      atm->q[ctl->qnt_p][ip] = atm->p[ip];

    /* Set geopotential height... */
    if (ctl->qnt_z >= 0)
      atm->q[ctl->qnt_z][ip] = z;

    /* Set temperature... */
    if (ctl->qnt_t >= 0)
      atm->q[ctl->qnt_t][ip] = t;

    /* Set zonal wind... */
    if (ctl->qnt_u >= 0)
      atm->q[ctl->qnt_u][ip] = u;

    /* Set meridional wind... */
    if (ctl->qnt_v >= 0)
      atm->q[ctl->qnt_v][ip] = v;

    /* Set vertical velocity... */
    if (ctl->qnt_w >= 0)
      atm->q[ctl->qnt_w][ip] = w;

    /* Set water vapor vmr... */
    if (ctl->qnt_h2o >= 0)
      atm->q[ctl->qnt_h2o][ip] = h2o;

    /* Set ozone vmr... */
    if (ctl->qnt_o3 >= 0)
      atm->q[ctl->qnt_o3][ip] = o3;

    /* Set cloud liquid water content... */
    if (ctl->qnt_lwc >= 0)
      atm->q[ctl->qnt_lwc][ip] = lwc;

    /* Set cloud ice water content... */
    if (ctl->qnt_iwc >= 0)
      atm->q[ctl->qnt_iwc][ip] = iwc;

    /* Set cloud top pressure... */
    if (ctl->qnt_pc >= 0)
      atm->q[ctl->qnt_pc][ip] = pc;

    /* Set nitric acid vmr... */
    if (ctl->qnt_hno3 >= 0)
      atm->q[ctl->qnt_hno3][ip] =
	clim_hno3(atm->time[ip], atm->lat[ip], atm->p[ip]);

    /* Set hydroxyl number concentration... */
    if (ctl->qnt_oh >= 0)
      atm->q[ctl->qnt_oh][ip] =
	clim_oh(atm->time[ip], atm->lat[ip], atm->p[ip]);

    /* Calculate horizontal wind... */
    if (ctl->qnt_vh >= 0)
      atm->q[ctl->qnt_vh][ip] = sqrt(u * u + v * v);

    /* Calculate vertical velocity... */
    if (ctl->qnt_vz >= 0)
      atm->q[ctl->qnt_vz][ip] = -1e3 * H0 / atm->p[ip] * w;

    /* Calculate relative humidty... */
    if (ctl->qnt_rh >= 0)
      atm->q[ctl->qnt_rh][ip] = RH(atm->p[ip], t, h2o);

    /* Calculate potential temperature... */
    if (ctl->qnt_theta >= 0)
      atm->q[ctl->qnt_theta][ip] = THETA(atm->p[ip], t);

    /* Set potential vorticity... */
    if (ctl->qnt_pv >= 0)
      atm->q[ctl->qnt_pv][ip] = pv;

    /* Calculate T_ice (Marti and Mauersberger, 1993)... */
    if (ctl->qnt_tice >= 0)
      atm->q[ctl->qnt_tice][ip] =
	-2663.5 /
	(log10((ctl->psc_h2o > 0 ? ctl->psc_h2o : h2o) * atm->p[ip] * 100.) -
	 12.537);

    /* Calculate T_NAT (Hanson and Mauersberger, 1988)... */
    if (ctl->qnt_tnat >= 0) {
      double p_hno3;
      if (ctl->psc_hno3 > 0)
	p_hno3 = ctl->psc_hno3 * atm->p[ip] / 1.333224;
      else
	p_hno3 = clim_hno3(atm->time[ip], atm->lat[ip], atm->p[ip])
	  * 1e-9 * atm->p[ip] / 1.333224;
      double p_h2o =
	(ctl->psc_h2o > 0 ? ctl->psc_h2o : h2o) * atm->p[ip] / 1.333224;
      double a = 0.009179 - 0.00088 * log10(p_h2o);
      double b = (38.9855 - log10(p_hno3) - 2.7836 * log10(p_h2o)) / a;
      double c = -11397.0 / a;
      double x1 = (-b + sqrt(b * b - 4. * c)) / 2.;
      double x2 = (-b - sqrt(b * b - 4. * c)) / 2.;
      if (x1 > 0)
	atm->q[ctl->qnt_tnat][ip] = x1;
      if (x2 > 0)
	atm->q[ctl->qnt_tnat][ip] = x2;
    }

    /* Calculate T_STS (mean of T_ice and T_NAT)... */
    if (ctl->qnt_tsts >= 0)
      atm->q[ctl->qnt_tsts][ip] = 0.5 * (atm->q[ctl->qnt_tice][ip]
					 + atm->q[ctl->qnt_tnat][ip]);
  }
}

/*****************************************************************************/

void module_position(
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt) {

#ifdef _OPENACC
#pragma acc data present(met0,met1,atm,dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++)
    if (dt[ip] != 0) {

      double ps, cw[3];

      int ci[3];

      /* Calculate modulo... */
      atm->lon[ip] = FMOD(atm->lon[ip], 360.);
      atm->lat[ip] = FMOD(atm->lat[ip], 360.);

      /* Check latitude... */
      while (atm->lat[ip] < -90 || atm->lat[ip] > 90) {
	if (atm->lat[ip] > 90) {
	  atm->lat[ip] = 180 - atm->lat[ip];
	  atm->lon[ip] += 180;
	}
	if (atm->lat[ip] < -90) {
	  atm->lat[ip] = -180 - atm->lat[ip];
	  atm->lon[ip] += 180;
	}
      }

      /* Check longitude... */
      while (atm->lon[ip] < -180)
	atm->lon[ip] += 360;
      while (atm->lon[ip] >= 180)
	atm->lon[ip] -= 360;

      /* Check pressure... */
      if (atm->p[ip] < met0->p[met0->np - 1])
	atm->p[ip] = met0->p[met0->np - 1];
      else if (atm->p[ip] > 300.) {
	intpol_met_time_2d(met0, met0->ps, met1, met1->ps, atm->time[ip],
			   atm->lon[ip], atm->lat[ip], &ps, ci, cw, 1);
	if (atm->p[ip] > ps)
	  atm->p[ip] = ps;
      }
    }
}

/*****************************************************************************/

void module_sedi(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt) {

#ifdef _OPENACC
#pragma acc data present(ctl,met0,met1,atm,dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++)
    if (dt[ip] != 0) {

      double G, K, eta, lambda, p, r_p, rho, rho_p, T, v, v_p, cw[3];

      int ci[3];

      /* Convert units... */
      p = 100. * atm->p[ip];
      r_p = 1e-6 * atm->q[ctl->qnt_r][ip];
      rho_p = atm->q[ctl->qnt_rho][ip];

      /* Get temperature... */
      intpol_met_time_3d(met0, met0->t, met1, met1->t, atm->time[ip],
			 atm->p[ip], atm->lon[ip], atm->lat[ip], &T, ci, cw,
			 1);

      /* Density of dry air... */
      rho = p / (RA * T);

      /* Dynamic viscosity of air... */
      eta = 1.8325e-5 * (416.16 / (T + 120.)) * pow(T / 296.16, 1.5);

      /* Thermal velocity of an air molecule... */
      v = sqrt(8. * KB * T / (M_PI * 4.8096e-26));

      /* Mean free path of an air molecule... */
      lambda = 2. * eta / (rho * v);

      /* Knudsen number for air... */
      K = lambda / r_p;

      /* Cunningham slip-flow correction... */
      G = 1. + K * (1.249 + 0.42 * exp(-0.87 / K));

      /* Sedimentation (fall) velocity... */
      v_p = 2. * SQR(r_p) * (rho_p - rho) * G0 / (9. * eta) * G;

      /* Calculate pressure change... */
      atm->p[ip] += DZ2DP(v_p * dt[ip] / 1000., atm->p[ip]);
    }
}

/*****************************************************************************/

void module_oh_chem(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt) {

  /* Check quantity flags... */
  if (ctl->qnt_m < 0)
    ERRMSG("Module needs quantity mass!");

#ifdef _OPENACC
#pragma acc data present(ctl,atm,dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++)
    if (dt[ip] != 0) {

      double c, k, k0, ki, M, T, cw[3];

      int ci[3];

      /* Get temperature... */
      intpol_met_time_3d(met0, met0->t, met1, met1->t, atm->time[ip],
			 atm->p[ip], atm->lon[ip], atm->lat[ip], &T, ci, cw,
			 1);

      /* Calculate molecular density... */
      M = 7.243e21 * (atm->p[ip] / P0) / T;

      /* Calculate rate coefficient for X + OH + M -> XOH + M
         (JPL Publication 15-10) ... */
      k0 = ctl->oh_chem[0] *
	(ctl->oh_chem[1] > 0 ? pow(T / 300., -ctl->oh_chem[1]) : 1.);
      ki = ctl->oh_chem[2] *
	(ctl->oh_chem[3] > 0 ? pow(T / 300., -ctl->oh_chem[3]) : 1.);
      c = log10(k0 * M / ki);
      k = k0 * M / (1. + k0 * M / ki) * pow(0.6, 1. / (1. + c * c));

      /* Calculate exponential decay... */
      atm->q[ctl->qnt_m][ip] *=
	exp(-dt[ip] * k * clim_oh(atm->time[ip], atm->lat[ip], atm->p[ip]));
    }
}

/*****************************************************************************/

void module_wet_deposition(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double *dt) {

  /* Check quantity flags... */
  if (ctl->qnt_m < 0)
    ERRMSG("Module needs quantity mass!");

#ifdef _OPENACC
#pragma acc data present(ctl,atm,dt)
#pragma acc parallel loop independent gang vector
#else
#pragma omp parallel for default(shared)
#endif
  for (int ip = 0; ip < atm->np; ip++)
    if (dt[ip] != 0) {

      double H, Is, Si, T, cl, lambda, iwc, lwc, pc, cw[3];

      int inside, ci[3];

      /* Check whether particle is below cloud top... */
      intpol_met_time_2d(met0, met0->pc, met1, met1->pc, atm->time[ip],
			 atm->lon[ip], atm->lat[ip], &pc, ci, cw, 1);
      if (!isfinite(pc) || atm->p[ip] <= pc)
	continue;

      /* Check whether particle is inside or below cloud... */
      intpol_met_time_3d(met0, met0->lwc, met1, met1->lwc, atm->time[ip],
			 atm->p[ip], atm->lon[ip], atm->lat[ip], &lwc, ci, cw,
			 1);
      intpol_met_time_3d(met0, met0->iwc, met1, met1->iwc, atm->time[ip],
			 atm->p[ip], atm->lon[ip], atm->lat[ip], &iwc, ci, cw,
			 0);
      inside = (iwc > 0 || lwc > 0);

      /* Estimate precipitation rate (Pisso et al., 2019)... */
      intpol_met_time_2d(met0, met0->cl, met1, met1->cl, atm->time[ip],
			 atm->lon[ip], atm->lat[ip], &cl, ci, cw, 0);
      Is = pow(2. * cl, 1. / 0.36);
      if (Is < 0.01)
	continue;

      /* Calculate in-cloud scavenging for gases... */
      if (inside) {

	/* Get temperature... */
	intpol_met_time_3d(met0, met0->t, met1, met1->t, atm->time[ip],
			   atm->p[ip], atm->lon[ip], atm->lat[ip], &T, ci, cw,
			   0);

	/* Get Henry's constant (Sander, 2015)... */
	H = ctl->wet_depo[2] * 101.325
	  * exp(ctl->wet_depo[3] * (1. / T - 1. / 298.15));

	/* Get scavenging coefficient (Hertel et al., 1995)... */
	Si = 1. / ((1. - cl) / (H * RI / P0 * T) + cl);
	lambda = 6.2 * Si * Is / 3.6e6;
      }

      /* Calculate below-cloud scavenging for gases (Pisso et al., 2019)... */
      else
	lambda = ctl->wet_depo[0] * pow(Is, ctl->wet_depo[1]);

      /* Calculate exponential decay... */
      atm->q[ctl->qnt_m][ip] *= exp(-dt[ip] * lambda);
    }
}

/*****************************************************************************/

void write_output(
  const char *dirname,
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  double t) {

  char filename[2 * LEN];

  double r;

  int year, mon, day, hour, min, sec, updated = 0;

  /* Get time... */
  jsec2time(t, &year, &mon, &day, &hour, &min, &sec, &r);

  /* Write atmospheric data... */
  if (ctl->atm_basename[0] != '-' && fmod(t, ctl->atm_dt_out) == 0) {
    if (!updated) {
      RANGE_PUSH("W atm D2H", RED);
#ifdef _OPENACC
#pragma acc update host(atm[:1])
#endif
      RANGE_POP;
      updated = 1;
    }
    RANGE_PUSH("IO", YELLOW);
    sprintf(filename, "%s/%s_%04d_%02d_%02d_%02d_%02d.tab",
	    dirname, ctl->atm_basename, year, mon, day, hour, min);
    write_atm(filename, ctl, atm, t);
    RANGE_POP;
  }

  /* Write gridded data... */
  if (ctl->grid_basename[0] != '-' && fmod(t, ctl->grid_dt_out) == 0) {
    if (!updated) {
      RANGE_PUSH("W grd D2H", RED);
#ifdef _OPENACC
#pragma acc update host(atm[:1])
#endif
      RANGE_POP;
      updated = 1;
    }
    RANGE_PUSH("IO", YELLOW);
    sprintf(filename, "%s/%s_%04d_%02d_%02d_%02d_%02d.tab",
	    dirname, ctl->grid_basename, year, mon, day, hour, min);
    write_grid(filename, ctl, met0, met1, atm, t);
    RANGE_POP;
  }

  /* Write CSI data... */
  if (ctl->csi_basename[0] != '-') {
    if (!updated) {
      RANGE_PUSH("W csi D2H", RED);
#ifdef _OPENACC
#pragma acc update host(atm[:1])
#endif
      RANGE_POP;
      updated = 1;
    }
    RANGE_PUSH("IO", YELLOW);
    sprintf(filename, "%s/%s.tab", dirname, ctl->csi_basename);
    write_csi(filename, ctl, atm, t);
    RANGE_POP;
  }

  /* Write ensemble data... */
  if (ctl->ens_basename[0] != '-') {
    if (!updated) {
      RANGE_PUSH("W csi D2H", RED);
#ifdef _OPENACC
#pragma acc update host(atm[:1])
#endif
      RANGE_POP;
      updated = 1;
    }
    RANGE_PUSH("IO", YELLOW);
    sprintf(filename, "%s/%s.tab", dirname, ctl->ens_basename);
    write_ens(filename, ctl, atm, t);
    RANGE_POP;
  }

  /* Write profile data... */
  if (ctl->prof_basename[0] != '-') {
    if (!updated) {
      RANGE_PUSH("W prof D2H", RED);
#ifdef _OPENACC
#pragma acc update host(atm[:1])
#endif
      RANGE_POP;
      updated = 1;
    }
    RANGE_PUSH("IO", YELLOW);
    sprintf(filename, "%s/%s.tab", dirname, ctl->prof_basename);
    write_prof(filename, ctl, met0, met1, atm, t);
    RANGE_POP;
  }

  /* Write station data... */
  if (ctl->stat_basename[0] != '-') {
    if (!updated) {
      RANGE_PUSH("W st D2H", RED);
#ifdef _OPENACC
#pragma acc update host(atm[:1])
#endif
      RANGE_POP;
      updated = 1;
    }
    RANGE_PUSH("IO", YELLOW);
    sprintf(filename, "%s/%s.tab", dirname, ctl->stat_basename);
    write_station(filename, ctl, atm, t);
    RANGE_POP;
  }
}
