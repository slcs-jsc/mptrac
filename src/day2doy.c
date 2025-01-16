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
  
  Copyright (C) 2013-2025 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Convert date to day of year.
*/

#include "mptrac.h"

int main(
  int argc,
  char *argv[]) {

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <year> <mon> <day>");

  /* Read arguments... */
  const int year = atoi(argv[1]);
  const int mon = atoi(argv[2]);
  const int day = atoi(argv[3]);

  /* Convert... */
  int doy;
  day2doy(year, mon, day, &doy);
  printf("%d %d\n", year, doy);

  return EXIT_SUCCESS;
}
