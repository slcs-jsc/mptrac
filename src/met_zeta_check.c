#include "libtrac.h"

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  ctl_t *ctl;
  ALLOC(ctl, ctl_t, 1);
  ctl->read_mode = 0;
  ctl->chunkszhint = 160000000;
  ctl->atm_type = 3;
  ctl->clams_met_data = 1;
  ctl->vert_coord_met = 1;
  ctl->vert_vel = 1;

  met_t *met;
  ALLOC(met, met_t, 1);
  //char filename[LEN]="/p/fastdata/slmet/slmet111/model_data/mptrac/input/pressure_0.3deg_v2/nc/2016/07/era5_2016_07_12_07.nc";
  int ncid = 0;

  /* Check arguments... */
  if (argc < 1)
    ERRMSG("Give parameters: <met>");

  /* Open netCDF file... */
  if (nc__open(argv[1], ctl->read_mode, &ctl->chunkszhint, &ncid) != NC_NOERR) {
    WARN("File not found!");
    return 0;
  }

  read_met_grid(argv[1], ncid, ctl, met);
  read_met_levels(ncid, ctl, met);


  for (int i = 0; i < met->nx; i++)
    for (int j = 0; j < met->ny; j++)
      for (int k = 1; k < met->np; k++)
	if ((met->zeta[i][j][k - 1] >=
	     met->zeta[i][j][k]) & (met->zeta[i][j][k - 1] >
				    0.0) & (met->zeta[i][j][k] > 0.0)) {
	  int l = 0;
	  printf("1:%f,%f\n", met->zeta[i][j][k - 1], met->zeta[i][j][k]);
	  while (met->zeta[i][j][k - 1] >= met->zeta[i][j][k + l])
	    l = l + 1;
	  float w =
	    (float) (met->p[k] - met->p[k - 1]) / (float) (met->p[k + l] -
							   met->p[k - 1]);
	  float d = (met->zeta[i][j][k + l] - met->zeta[i][j][k - 1]);
	  met->zeta[i][j][k] = d * w + met->zeta[i][j][k - 1];
	  printf("2:%f,%f\n", met->zeta[i][j][k - 1], met->zeta[i][j][k]);
	}


  int cnt = 0;
  for (int i = 0; i < met->nx; i++)
    for (int j = 0; j < met->ny; j++)
      for (int k = 1; k < met->np; k++)
	if ((met->zeta[i][j][k - 1] >
	     met->zeta[i][j][k]) & (met->zeta[i][j][k - 1] >
				    0.0) & (met->zeta[i][j][k] > 0.0)) {
	  cnt = cnt + 1;

	  //printf("%f,%f,%f,%f,%d,%f,%f\n",met->zeta[i][j][k],met->zeta[i][j][k-1],met->p[k],met->p[k-1],k,met->w[i][j][k],met->w[i][j][k-1]);
	  //break;
	}

  printf("%d\n", cnt);
  free(met);
  free(ctl);
}
