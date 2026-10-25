/* Compare the complete local meteorological state produced by both DD readers.
 * Build with the same dimensions and options as the MPTRAC executable. */
#include "mptrac.h"

int main(int argc, char **argv) {

#ifndef DD
  (void) argc;
  (void) argv;
  ERRMSG("This comparison requires DD and MPI!");
#else
  MPI_Init(&argc, &argv);
  if (argc != 5)
    ERRMSG("Usage: met_reader_compare config.ctl wind.nc nx_domains ny_domains");

  ctl_t *ctl;
  met_t *old, *global;
  dd_t *old_dd, *global_dd;
  ALLOC(ctl, ctl_t, 1);
  ALLOC(old, met_t, 1);
  ALLOC(global, met_t, 1);
  ALLOC(old_dd, dd_t, 1);
  ALLOC(global_dd, dd_t, 1);

  /* The ordinary 2x2 DD control file supplies all meteorological settings.
     Override only the layout after its regular validation. */
  char *ctl_argv[] = { argv[0], argv[1] };
  mptrac_read_ctl(argv[1], 2, ctl_argv, ctl);
  ctl->dd = 1;
  ctl->dd_subdomains_zonal = atoi(argv[3]);
  ctl->dd_subdomains_meridional = atoi(argv[4]);
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (size != ctl->dd_subdomains_zonal * ctl->dd_subdomains_meridional)
    ERRMSG("Layout does not match MPI size!");

  ctl->met_nc_read_mode = 0;
  if (!read_met_nc(argv[2], ctl, old, old_dd))
    ERRMSG("Existing reader failed!");
  ctl->met_nc_read_mode = 1;
  if (!read_met_nc(argv[2], ctl, global, global_dd))
    ERRMSG("Global 2-D reader failed!");

  const unsigned char *a = (const unsigned char *) old;
  const unsigned char *b = (const unsigned char *) global;
  size_t differences = 0;
  for (size_t i = 0; i < sizeof(met_t); i++)
    if (a[i] != b[i]) {
      if (differences == 0)
	fprintf(stderr, "Rank %d: first meteorological byte mismatch at %zu "
		"(%u versus %u)\n", rank, i, a[i], b[i]);
      differences++;
    }
  printf("MET_READER_COMPARE rank=%d layout=%dx%d nx=%d ny=%d np=%d "
	 "different_bytes=%zu\n", rank, ctl->dd_subdomains_zonal,
	 ctl->dd_subdomains_meridional, old->nx, old->ny, old->np,
	 differences);
  free(global_dd);
  free(old_dd);
  free(global);
  free(old);
  free(ctl);
  MPI_Finalize();
  return differences ? EXIT_FAILURE : EXIT_SUCCESS;
#endif
}
