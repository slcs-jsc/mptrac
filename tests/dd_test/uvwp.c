/* Direct MPI regression for particle-specific turbulent velocity state. */
#include "mptrac.h"

static float expected_uvwp(
  const int id,
  const int component) {
  return (float) (10 * id + component + 1) + 0.25f;
}

int main(
  int argc,
  char *argv[]) {

  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2 || argc != 2)
    MPI_Abort(MPI_COMM_WORLD, 1);

  const int sort = strcmp(argv[1], "sort") == 0;
  const int push = strcmp(argv[1], "push") == 0;
  const int compact = strcmp(argv[1], "compact") == 0;
  if (!sort && !push && !compact)
    MPI_Abort(MPI_COMM_WORLD, 2);

#ifdef _OPENACC
  const int devices = acc_get_num_devices(acc_device_nvidia);
  if (devices < 1)
    MPI_Abort(MPI_COMM_WORLD, 3);
  acc_set_device_num(rank % devices, acc_device_nvidia);
#endif

  ctl_t *ctl;
  cache_t *cache;
  met_t *met;
  atm_t *atm;
  dd_t *dd;
  mptrac_alloc(&ctl, &cache, NULL, &met, NULL, &atm, NULL, &dd);

  ctl->dd = 1;
  ctl->dd_subdomains_zonal = 2;
  ctl->dd_subdomains_meridional = 1;
  ctl->nq = 3;
  ctl->qnt_idx = 0;
  ctl->qnt_current_subdomain = 1;
  ctl->qnt_target_subdomain = 2;
  ctl->dt_mod = 600;

  dd->nx_glob = 4;
  dd->ny_glob = 4;
  for (int i = 0; i < 4; i++) {
    dd->lon_glob[i] = 90. * i;
    dd->lat_glob[i] = -60. + 40. * i;
    met->lon[i] = dd->lon_glob[i];
    met->lat[i] = dd->lat_glob[i];
  }
  met->nx = 4;
  met->ny = 4;
  met->np = 2;
  met->p[0] = 100;
  met->p[1] = 500;

  atm->np = 4;
  for (int ip = 0; ip < atm->np; ip++) {
    const int id = 4 * rank + ip;
    atm->time[ip] = 1000. + id;
    atm->p[ip] = 300;
    atm->lon[ip] = rank == 0 ? 20. + 20. * ip : 200. + 20. * ip;
    atm->lat[ip] = 0;
    atm->q[0][ip] = id;
    cache->dt[ip] = 120;
    for (int j = 0; j < 3; j++)
      cache->uvwp[ip][j] = expected_uvwp(id, j);
  }

  mptrac_update_device(ctl, cache, NULL, &met, NULL, atm);
  dd_init(ctl, dd, atm);

  MPI_Aint lower_bound, extent, true_lower_bound, true_extent;
  MPI_Type_get_extent(dd->MPI_Particle, &lower_bound, &extent);
  MPI_Type_get_true_extent(dd->MPI_Particle, &true_lower_bound,
			   &true_extent);
  int type_size;
  MPI_Type_size(dd->MPI_Particle, &type_size);
  if (lower_bound != 0 || true_lower_bound != 0
      || extent != (MPI_Aint) sizeof(particle_t)
      || true_extent != (MPI_Aint) (offsetof(particle_t, uvwp)
				     + sizeof(((particle_t *) 0)->uvwp))
      || type_size != (int) ((4 + NQ) * sizeof(double)
			    + 3 * sizeof(float)))
    MPI_Abort(MPI_COMM_WORLD, 4);

  mptrac_update_host(NULL, NULL, NULL, NULL, NULL, atm);
  for (int ip = 0; ip < atm->np; ip++)
    if (ip % 2 == 1)
      atm->lon[ip] = rank == 0 ? 220. : 40.;
  mptrac_update_device(NULL, NULL, NULL, NULL, NULL, atm);
  dd_assign_subdomains(ctl, dd, atm, 0);

  int npart = 0, nkeep = 0, capacity = NP;
  particle_t *particles = calloc((size_t) capacity, sizeof(particle_t));
  if (particles == NULL)
    MPI_Abort(MPI_COMM_WORLD, 5);

  if (sort)
    dd_sort(ctl, met, cache, atm, dd, &npart);
  else if (push)
    dd_push(ctl, atm, cache, &npart);
  else
    dd_count_outgoing(ctl, atm, dd, &nkeep, &npart);
  if (npart != 2)
    MPI_Abort(MPI_COMM_WORLD, 6);

  if (compact)
    dd_pack_compact(ctl, cache, atm, dd, particles, nkeep, npart);
  else
    dd_atm2particles(ctl, cache, atm, particles, npart);

  dd_communicate_particles(ctl, dd, &particles, &npart, &capacity);
  if (npart != 2)
    MPI_Abort(MPI_COMM_WORLD, 7);
  dd_particles2atm(ctl, cache, particles, npart, atm);
  mptrac_update_host(NULL, cache, NULL, NULL, NULL, atm);

  int seen[8] = { 0 };
  for (int ip = 0; ip < atm->np; ip++) {
    if ((int) atm->q[ctl->qnt_current_subdomain][ip] < 0)
      continue;
    const int id = (int) atm->q[ctl->qnt_idx][ip];
    if (id < 0 || id >= 8 || seen[id]++)
      MPI_Abort(MPI_COMM_WORLD, 8);
    if ((id / 4 == rank && id % 2 != 0)
	|| (id / 4 != rank && id % 2 != 1))
      MPI_Abort(MPI_COMM_WORLD, 9);
    for (int j = 0; j < 3; j++) {
      const float expected = expected_uvwp(id, j);
      if (memcmp(&cache->uvwp[ip][j], &expected, sizeof(float)) != 0)
	MPI_Abort(MPI_COMM_WORLD, 10);
    }
  }
  int count = 0;
  for (int id = 0; id < 8; id++)
    count += seen[id];
  if (count != 4)
    MPI_Abort(MPI_COMM_WORLD, 11);

  free(particles);
  mptrac_free(ctl, cache, NULL, met, NULL, atm, NULL, dd);
  MPI_Finalize();
  return 0;
}
