#include "libtrac.h"

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  ctl_t *ctl;
  ALLOC(ctl, ctl_t, 1);

  met_t *met0;
  ALLOC(met0, met_t, 1);
  met_t *met1;
  ALLOC(met1, met_t, 1);
  int ncid;

  float pressures0[EX][EY][EP];
  //ALLOC(pressures0,float,EX*EY*EP);
  float pressures1[EX][EY][EP];
  //ALLOC(pressures1,float,EX*EY*EP);

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <met0> <met1> <ctl>");

  /* Read control parameters... */
  read_ctl(argv[3], argc, argv, ctl);

  ctl->read_mode = 0;
  ctl->chunkszhint = 163840000;
  ctl->atm_type = 3;
  ctl->clams_met_data = 1;
  ctl->vert_coord_met = 1;
  ctl->vert_vel = 1;

  /* Open netCDF file... */
  if (nc__open(argv[1], ctl->read_mode, &ctl->chunkszhint, &ncid) != NC_NOERR) {
    WARN("File not found!");
    return 0;
  }

  /* Read meteo data on vertical levels... */
  read_met_grid(argv[1], ncid, ctl, met0);
  read_met_levels(ncid, ctl, met0);

  /* Open netCDF file... */
  if (nc__open(argv[2], ctl->read_mode, &ctl->chunkszhint, &ncid) != NC_NOERR) {
    WARN("File not found!");
    return 0;
  }

  /* Read meteo data on vertical levels... */
  read_met_grid(argv[2], ncid, ctl, met1);
  read_met_levels(ncid, ctl, met1);

  double lon_ap = 90;
  double lat_ap = 45;
  double time_ap = 521399800.00;
  //double zeta_ap = 400;
  double press = 0;
  double zeta_out = 0;


  for (int i = 0; i < met0->nx; i++)
    for (int j = 0; j < met0->ny; j++)
      for (int k = 0; k < met0->np; k++) {
	//printf("%d,%d,%d\n",i,j,k);
	pressures0[i][j][k] = (float) met0->p[k];
	pressures1[i][j][k] = (float) met1->p[k];
      }


  for (double zeta_ap = 10; zeta_ap < 2000; zeta_ap++) {
    //press = intpol_ap_ml2pl_time(met0,met1,lon_ap,lat_ap,zeta_ap,time_ap);
    INTPOL_INIT;
    double cw_apc[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    int ci_apc[3] = { 0, 0, 0 };

    intpol_met_time_3d_ap_coord(met0, pressures0, met1, pressures1, time_ap,
				zeta_ap, lon_ap, lat_ap, &press, ci_apc,
				cw_apc, 1);
    intpol_met_time_3d(met0, met0->zeta, met1, met1->zeta, time_ap, press,
		       lon_ap, lat_ap, &zeta_out, ci, cw, 1);

    printf("%f;%f;%f;%f;%f\n", zeta_ap, press, zeta_out, zeta_out - zeta_ap,
	   met1->zeta_dot[ci[0]][ci[1]][ci[2]]);
  }

  free(met0);
  free(met1);
  //free(pressures0);
  //free(pressures1);
  free(ctl);

}

	  //int ci[3] = {0, 0, 0};
	//double cw_apc[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; 
	//intpol_met_space_3d_ap_coord(met, met->zeta, press, lon_ap, lat_ap, &press,ci,cw_apc,1);
