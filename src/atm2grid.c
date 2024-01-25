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
  
  Copyright (C) 2013-2019 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Convert atmosphric data file to grid data file.
  example:
  cd data
  for atmfile in $(ls atm_20*_00_00.tab); do
  ../../../src/atm2grid trac.ctl $atmfile GRID_BASENAME grid
  done
*/

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  atm_t *atm;

  char filename[3 * LEN];

  met_t *met0 = NULL, *met1 = NULL;
  
  /* Allocate... */
  ALLOC(atm, atm_t, 1);

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <atm_in>");

  read_ctl(argv[1], argc, argv, &ctl);

  /* Read atmospheric data... */
  if (!read_atm(argv[2], &ctl, atm))
    ERRMSG("Cannot open file!");

  /* Set time interval for output... */
  int year, mon, day, hour, min;
  double t;
  sscanf(argv[2], "atm_%d_%d_%d_%d_%d.tab", &year, &mon, &day, &hour, &min);
  time2jsec(year, mon, day, hour, min, 0, 0, &t);
  
  sprintf(filename, "%s_%04d_%02d_%02d_%02d_%02d.%s",
    ctl.grid_basename, year, mon, day, hour, min,
    ctl.grid_type == 0 ? "tab" : "nc");
    
  write_grid(filename, &ctl, met0, met1, atm, t);

  return EXIT_SUCCESS;

}