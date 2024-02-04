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
  
  Copyright (C) 2013-2024 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Convert atmospheric data file to grid data file.
*/

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  atm_t *atm;

  met_t *met0 = NULL, *met1 = NULL;

  /* Allocate... */
  ALLOC(atm, atm_t, 1);

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <atm_in>");

  /* Read control parameters... */
  read_ctl(argv[1], argc, argv, &ctl);

  /* Check grid basename... */
  if (ctl.grid_basename[0] == '-')
    ERRMSG("You need to specify GRID_BASENAME!");

  /* Read atmospheric data... */
  if (!read_atm(argv[2], &ctl, atm))
    ERRMSG("Cannot open file!");

  /* Get time from filename... */
  int year, mon, day, hour, min, sec;
  double r, t = time_from_filename(argv[2], ctl.atm_type < 2 ? 20 : 19);
  jsec2time(t, &year, &mon, &day, &hour, &min, &sec, &r);

  /* Set output filename... */
  char filename[3 * LEN];
  sprintf(filename, "%s_%04d_%02d_%02d_%02d_%02d.%s",
	  ctl.grid_basename, year, mon, day, hour, min,
	  ctl.grid_type == 0 ? "tab" : "nc");

  /* Write grid data... */
  write_grid(filename, &ctl, met0, met1, atm, t);

  /* Free... */
  free(atm);

  return EXIT_SUCCESS;
}
