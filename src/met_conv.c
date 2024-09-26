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
  
  Copyright (C) 2013-2023 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Convert file format of meteo data files.
*/

#include "mptrac.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  clim_t *clim;

  met_t *met;

  /* Start timers... */
  START_TIMERS;

  /* Check arguments... */
  if (argc < 6)
    ERRMSG("Give parameters: <ctl> <met_in> <met_in_type>"
	   " <met_out> <met_out_type>");

  /* Allocate... */
  ALLOC(clim, clim_t, 1);
  ALLOC(met, met_t, 1);

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);

  /* Read climatological data... */
  read_clim(&ctl, clim);

  /* Read meteo data... */
  ctl.met_type = atoi(argv[3]);
  if (!read_met(argv[2], &ctl, clim, met))
    ERRMSG("Cannot open file!");

  /* Write meteo data... */
  ctl.met_type = atoi(argv[5]);
  write_met(argv[4], &ctl, met);

  /* Free... */
  free(clim);
  free(met);

  /* Report timers... */
  PRINT_TIMERS;
  STOP_TIMERS;

  return EXIT_SUCCESS;
}
