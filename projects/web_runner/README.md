# MPTRAC Web Runner

The MPTRAC Web Runner provides a browser-based interface for performing air parcel trajectory simulations with the [MPTRAC model](https://slcs-jsc.github.io/mptrac/).

It allows users to configure, execute, and visualize Lagrangian transport simulations using meteorological reanalysis and forecast data through an interactive web interface.

An instance of the MPTRAC Web Runner is currently running on the [JSC Cloud](https://www.fz-juelich.de/en/jsc/systems/scientific-clouds/jsc-cloud): [https://mptrac.jsc.fz-juelich.de/](https://mptrac.jsc.fz-juelich.de/)

## Features

### Simulation Setup and Configuration
- Web-based setup of initial positions, heights, and simulation duration (forward or backward in time).  
- Supports up to 10,000 air parcels per simulation.  
- Includes predefined example cases for events such as volcanic eruptions, wildfires, dust storms, and pollution plumes.  

### Meteorological Data
- Supports multiple datasets (ECMWF AIFS/IFS, ERA5, ERA-Interim, JRA-3Q, JRA-55, MERRA-2, NCEP, NCEP2).  
- Automatically manages dataset metadata, vertical coordinates, and time resolution.  
- Ensures simulation time spans do not exceed model constraints (up to 30 days).

### Model Execution
- Automatically generates required MPTRAC input files (`trac.ctl`, `dirlist.txt`, `atm_init.tab`).  
- Runs MPTRAC binaries (`atm_init`, `trac`) via subprocesses.  
- Controls concurrent execution using semaphores and performs automatic cleanup of temporary files.  

### Physics and Parameterization
- Adjustable diffusivity for PBL, troposphere, and stratosphere.  
- Configurable subgrid-scale winds, convective parametrization thresholds (CAPE/CIN), and PBL mixing.
- Randomized uniform or Gaussian perturbations for initial parcel distributions.

### Output and Visualization
- Configurable output frequency (15 min to 24 h).  
- Generates trajectory plots using Matplotlib and Cartopy with multiple projections and color mapping by altitude.  
- Interactive in-browser visualization with image playback controls.  
- Provides ZIP downloads of all simulation outputs and logs.

### Web Application and Error Handling
- Implemented in Flask with modular Jinja2 templates.  
- Responsive layout using W3.CSS and accessible dark-themed design.  
- Distinct pages for results, errors, and busy-server states.  
- Execution logs stored in `runs/logs/web_runner.log`.

## Project Structure

```
app.py              Flask backend and simulation logic
requirements.txt    Python dependencies
job.sh              Local test and run script
templates/          HTML templates for UI and results
static/             Static assets (images, styles)
runs/               Working directories, logs, and output archives
```

## Installation

1. Install and compile MPTRAC following the official installation instructions.  
2. Configure meteorological data paths in `app.py` by updating the `MET_OPTIONS` dictionary to match your system’s file locations.  
3. Navigate to `[mptrac_dir]/projects/web_runner/` and set up a Python environment.  
4. Run locally with `python app.py` or `./job.sh`.  
5. For production, deploy with Gunicorn behind Nginx.

## Usage

1. Open the web interface (`http://localhost:5000` or your server URL).  
2. Configure simulation parameters and submit a job.  
3. Download results as ZIP archives or explore trajectory plots interactively.

## Citation

If you use the MPTRAC Web Runner or the MPTRAC model in your research, please cite:

Hoffmann, L., Clemens, J., Griessbach, S., Haghighi Mood, K., Heng, Y., Khosrawi, F., Liu, M., Lu, Y.-S., Meyer, C., Wittwer, N. N., Wu, X., & Zou, L. (2025). MPTRAC: A high-performance Lagrangian transport model for atmospheric air parcel dispersion. Journal of Open Source Software, 10 (111), 8177. [https://doi.org/10.21105/joss.08177](https://doi.org/10.21105/joss.08177)

## License

The MPTRAC Web Runner is released under the GNU General Public License v3 (GPLv3), consistent with the [MPTRAC model](https://github.com/slcs-jsc/mptrac). You are free to use, modify, and redistribute the software under the same license terms.

## Contact

**Dr. Lars Hoffmann**

Jülich Supercomputing Centre (JSC), Forschungszentrum Jülich, Germany

Email: [l.hoffmann@fz-juelich.de](mailto:l.hoffmann@fz-juelich.de)
