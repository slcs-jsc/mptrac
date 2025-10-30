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
  Spectral analysis of meteorological data.
*/

#include "mptrac.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  clim_t *clim;

  met_t *met;

  dd_t *dd;

  FILE *out;

  static double cutImag[EX], cutReal[EX], lx[EX], A[EX], phi[EX];

  /* Allocate... */
  ALLOC(clim, clim_t, 1);
  ALLOC(met, met_t, 1);
  ALLOC(dd, dd_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <spec.tab> <met0>");

  /* Read control parameters... */
  mptrac_read_ctl(argv[1], argc, argv, &ctl);
  const double wavemax =
    (int) scan_ctl(argv[1], argc, argv, "SPEC_WAVEMAX", -1, "7", NULL);

  /* Read climatological data... */
  mptrac_read_clim(&ctl, clim);

  /* Read meteorological data... */
  if (!mptrac_read_met(argv[3], &ctl, clim, met, dd))
    ERRMSG("Cannot read meteo data!");

  /* Create output file... */
  LOG(1, "Write spectral data file: %s", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = altitude [km]\n"
	  "# $3 = longitude [deg]\n" "# $4 = latitude [deg]\n");
  for (int ix = 0; ix <= wavemax; ix++) {
    fprintf(out, "# $%d = wavelength (PW%d) [km]\n", 5 + 3 * ix, ix);
    fprintf(out, "# $%d = amplitude (PW%d) [K]\n", 6 + 3 * ix, ix);
    fprintf(out, "# $%d = phase (PW%d) [deg]\n", 7 + 3 * ix, ix);
  }

  /* Loop over pressure levels... */
  for (int ip = 0; ip < met->np; ip++) {

    /* Write output... */
    fprintf(out, "\n");

    /* Loop over latitudes... */
    for (int iy = 0; iy < met->ny; iy++) {

      /* Copy data... */
      for (int ix = 0; ix < met->nx; ix++) {
	cutReal[ix] = met->t[ix][iy][ip];
	cutImag[ix] = 0.0;
      }

      /* FFT... */
      fft_help(cutReal, cutImag, met->nx);

      /* 
         Get wavelength, amplitude, and phase:
         A(x) = A[0] + A[1] * cos(2 pi x / lx[1] + phi[1]) + A[2] * cos...
       */
      for (int ix = 0; ix < met->nx; ix++) {
	lx[ix] = DEG2DX(met->lon[met->nx - 1] - met->lon[0], met->lat[iy])
	  / ((ix < met->nx / 2) ? (double) ix : -(double) (met->nx - ix));
	A[ix] = (ix == 0 ? 1.0 : 2.0) / (met->nx)
	  * sqrt(gsl_pow_2(cutReal[ix]) + gsl_pow_2(cutImag[ix]));
	phi[ix] = RAD2DEG(atan2(cutImag[ix], cutReal[ix]));
      }

      /* Write data... */
      fprintf(out, "%.2f %g %g %g", met->time, Z(met->p[ip]), 0.0,
	      met->lat[iy]);
      for (int ix = 0; ix <= wavemax; ix++)
	fprintf(out, " %g %g %g", lx[ix], A[ix], phi[ix]);
      fprintf(out, "\n");
    }
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(clim);
  free(met);
  free(dd);

  return EXIT_SUCCESS;
}
