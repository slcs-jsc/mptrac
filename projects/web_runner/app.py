from flask import Flask, render_template, request, send_from_directory
import argparse
import os, time, uuid, glob, shutil, textwrap, zipfile, subprocess, multiprocessing
from datetime import datetime, timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import logging
import json
from logging.handlers import RotatingFileHandler
from concurrent.futures import ProcessPoolExecutor, as_completed

# MPTRAC tool paths and local test-mode defaults.
ATM_INIT_CMD = '../../src/atm_init'
TRAC_CMD = '../../src/trac'
ATM_CONV_CMD = '../../src/atm_conv'
MET_ACCESS_TIMEOUT = 10
MET_ACCESS_PATH = '/mnt/slmet_mnt/met_data/'
TEST_MET_SOURCE = 'erai_6h'
TEST_METBASE = '../../tests/data/ei'
TEST_DT_MET = 86400
TEST_START_TIME = '2011-06-05 00:00'
TEST_STOP_TIME = '2011-06-08 00:00'
TEST_MZ = '30.0'
TEST_UZ = '60.0'

# Directories for working files, downloads, and logs.
RUNS_DIR, ZIPS_DIR, LOG_DIR = 'runs/working', 'runs/zips', 'runs/logs'
# Create the required runtime directories once.
for path in [RUNS_DIR, ZIPS_DIR, LOG_DIR]:
    os.makedirs(path, exist_ok=True)

# Rotating log file for the web application.
log_file = os.path.join(LOG_DIR, 'web_runner.log')
logger = logging.getLogger('web_runner')
logger.setLevel(logging.INFO)
handler = RotatingFileHandler(log_file, maxBytes=10*1024*1024, backupCount=5)
handler.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s in %(threadName)s: %(message)s'))
logger.addHandler(handler)

# Global constants for output variables and timeouts.

# Clean up old run directories and ZIP files.
def clean_old_runs(max_age_sec=3600):
    # Current time used as the cleanup reference.
    now = time.time()
    for directory in [RUNS_DIR, ZIPS_DIR]:
        logger.info(f"[CLEANUP] Scanning {directory} for files older than {max_age_sec} seconds.")
        for path in glob.glob(os.path.join(directory, "*")):
            if now - os.path.getmtime(path) > max_age_sec:
                try:
                    logger.info(f"[CLEANUP] Removed old file or directory: {path}")
                    if os.path.isdir(path):
                        shutil.rmtree(path)
                    else:
                        os.remove(path)
                except Exception as e:
                    logger.warning(f"[CLEANUP] Failed to remove {path}: {e}")

# Create a ZIP archive quickly without compression.
def fast_zip_no_compression(zip_path, source_dir):
    # Recursively package all files in the run directory.
    with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_STORED) as zipf:
        for root, _, files in os.walk(source_dir):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, source_dir)
                zipf.write(file_path, arcname)

# Initialize the Flask application and feature flags.
app = Flask(__name__)
app.config['TEST_MODE'] = False

# Available meteorological data sources and their metadata.
MET_OPTIONS = {
    'aifs_6h': dict(METBASE='/mnt/slmet_mnt/met_data/ecmwf/open_data/data/aifs-single_YYYY_MM_DD/aifs-single',
                    DT_MET=21600, MET_VERT_COORD=0, MET_PRESS_LEVEL_DEF=-1, MET_NAME='ECMWF AIFS'),
    'ifs_6h': dict(METBASE='/mnt/slmet_mnt/met_data/ecmwf/open_data/data/ifs_YYYY_MM_DD/ifs',
                   DT_MET=21600, MET_VERT_COORD=0, MET_PRESS_LEVEL_DEF=-1, MET_NAME='ECMWF IFS'),
    'gfs_3h': dict(METBASE='/mnt/slmet_mnt/met_data/ncep/gfs/data_YYYY_MM_DD/gfs',
                   DT_MET=10800, MET_VERT_COORD=0, MET_PRESS_LEVEL_DEF=-1, MET_NAME='NCEP GFS'),
    'era5low_6h': dict(METBASE='/mnt/slmet_mnt/met_data/ecmwf/era5.1/resolution_1x1/nc/YYYY/MM/era5',
                       DT_MET=21600, MET_VERT_COORD=1, MET_PRESS_LEVEL_DEF=6, MET_NAME='ERA5'),
    'erai_6h': dict(METBASE='/mnt/slmet_mnt/met_data/ecmwf/era_interim/pressure_0.75deg_v2/nc/YYYY/ei',
                    DT_MET=21600, MET_VERT_COORD=0, MET_PRESS_LEVEL_DEF=-1, MET_NAME='ERA-Interim'),
    'jra3q_6h': dict(METBASE='/mnt/slmet_mnt/met_data/jma/jra3q/nc/YYYY/jra3q',
                     DT_MET=21600, MET_VERT_COORD=2, MET_PRESS_LEVEL_DEF=6, MET_NAME='JRA-3Q'),
    'jra55_6h': dict(METBASE='/mnt/slmet_mnt/met_data/jma/jra55/nc/YYYY/jra55',
                     DT_MET=21600, MET_VERT_COORD=4, MET_PRESS_LEVEL_DEF=6, MET_NAME='JRA-55'),
    'merra2_3h': dict(METBASE='/mnt/slmet_mnt/met_data/nasa/merra-2/hybrid/YYYY/merra2',
                      DT_MET=10800, MET_VERT_COORD=1, MET_PRESS_LEVEL_DEF=6, MET_NAME='MERRA-2'),
    'ncep_6h': dict(METBASE='/mnt/slmet_mnt/met_data/ncep/reanalysis/nc/YYYY/ncep',
                    DT_MET=21600, MET_VERT_COORD=0, MET_PRESS_LEVEL_DEF=-1, MET_NAME='NCEP-NCAR Reanalysis 1'),
    'ncep2_6h': dict(METBASE='/mnt/slmet_mnt/met_data/ncep/reanalysis2/nc/YYYY/ncep2',
                     DT_MET=21600, MET_VERT_COORD=0, MET_PRESS_LEVEL_DEF=-1, MET_NAME='NCEP-DOE Reanalysis 2')
}

# Output quantities offered in the web interface.
MAX_OUTPUT_QUANTITIES = 15
DEFAULT_OUTPUT_QUANTITIES = ['zg', 'p', 't', 'theta', 'u', 'v', 'w', 'pv', 'h2o', 'o3', 'ps', 'pbl', 'pt', 'cc']
OUTPUT_QUANTITY_GROUPS = [
    {
        'title': 'Core State',
        'quantities': [
            {'name': 'zg', 'label': 'Geopotential height'},
            {'name': 'p', 'label': 'Pressure'},
            {'name': 't', 'label': 'Temperature'},
            {'name': 'theta', 'label': 'Potential temperature'},
            {'name': 'rho', 'label': 'Air density'},
        ],
    },
    {
        'title': 'Winds and Dynamics',
        'quantities': [
            {'name': 'u', 'label': 'Zonal wind'},
            {'name': 'v', 'label': 'Meridional wind'},
            {'name': 'w', 'label': 'Vertical velocity'},
            {'name': 'vh', 'label': 'Horizontal velocity'},
            {'name': 'vz', 'label': 'Vertical velocity in m/s'},
            {'name': 'pv', 'label': 'Potential vorticity'},
            {'name': 'zeta_d', 'label': 'Diagnosed zeta coordinate'},
        ],
    },
    {
        'title': 'Humidity and Thermodynamics',
        'quantities': [
            {'name': 'h2o', 'label': 'Water vapor'},
            {'name': 'o3', 'label': 'Ozone'},
            {'name': 'sh', 'label': 'Specific humidity'},
            {'name': 'rh', 'label': 'Relative humidity'},
            {'name': 'rhice', 'label': 'Relative humidity over ice'},
            {'name': 'tvirt', 'label': 'Virtual temperature'},
            {'name': 'tdew', 'label': 'Dew point temperature'},
            {'name': 'tice', 'label': 'Frost point temperature'},
            {'name': 'psat', 'label': 'Saturation pressure over water'},
            {'name': 'psice', 'label': 'Saturation pressure over ice'},
            {'name': 'pw', 'label': 'Partial water vapor pressure'},
        ],
    },
    {
        'title': 'Surface and Boundary Layer',
        'quantities': [
            {'name': 'ps', 'label': 'Surface pressure'},
            {'name': 'pbl', 'label': 'Planetary boundary layer'},
            {'name': 'ts', 'label': 'Surface temperature'},
            {'name': 'zs', 'label': 'Surface height'},
            {'name': 'us', 'label': 'Surface zonal wind'},
            {'name': 'vs', 'label': 'Surface meridional wind'},
            {'name': 'ess', 'label': 'Eastward turbulent surface stress'},
            {'name': 'nss', 'label': 'Northward turbulent surface stress'},
            {'name': 'shf', 'label': 'Surface sensible heat flux'},
            {'name': 'lsm', 'label': 'Land-sea mask'},
            {'name': 'sst', 'label': 'Sea surface temperature'},
        ],
    },
    {
        'title': 'Tropopause',
        'quantities': [
            {'name': 'pt', 'label': 'Tropopause pressure'},
            {'name': 'tt', 'label': 'Tropopause temperature'},
            {'name': 'zt', 'label': 'Tropopause geopotential height'},
            {'name': 'h2ot', 'label': 'Tropopause water vapor'},
        ],
    },
    {
        'title': 'Clouds and Convection',
        'quantities': [
            {'name': 'cc', 'label': 'Cloud cover'},
            {'name': 'cape', 'label': 'CAPE'},
            {'name': 'cin', 'label': 'CIN'},
            {'name': 'plcl', 'label': 'Lifted condensation level'},
            {'name': 'plfc', 'label': 'Level of free convection'},
            {'name': 'pel', 'label': 'Equilibrium level'},
            {'name': 'lwc', 'label': 'Cloud liquid water content'},
            {'name': 'rwc', 'label': 'Cloud rain water content'},
            {'name': 'iwc', 'label': 'Cloud ice water content'},
            {'name': 'swc', 'label': 'Cloud snow water content'},
            {'name': 'pct', 'label': 'Cloud top pressure'},
            {'name': 'pcb', 'label': 'Cloud bottom pressure'},
            {'name': 'cl', 'label': 'Total column cloud water'},
        ],
    },
]
ALLOWED_OUTPUT_QUANTITIES = {
    quantity['name']
    for group in OUTPUT_QUANTITY_GROUPS
    for quantity in group['quantities']
}

# Run external processes with timeout handling and logging.
def run_command(cmd, timeout, run_id):
    # Actual subprocess invocation.
    try:
        logger.info(f"[EXEC] [{run_id}] Executing command: {' '.join(cmd)}")
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=timeout)
        logger.info(f"[EXEC] [{run_id}] Command finished with exit code {result.returncode}.")
        return result.returncode, result.stdout
    except subprocess.TimeoutExpired as e:
        logger.warning(f"[EXEC] [{run_id}] Command timed out after {timeout} seconds.")
        return -1, (e.stdout or "") + f"\n❌ Command timed out after {timeout} seconds."
    except Exception as e:
        logger.exception(f"[EXEC] [{run_id}] Command failed before completion: {' '.join(cmd)}")
        return -1, f"❌ Command failed before completion: {' '.join(cmd)}\n{e}"

# Convert a UTC timestamp into MPTRAC seconds.
def seconds_since_2000(time_str):
    dt = datetime.strptime(time_str, "%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc)
    return (dt - datetime(2000, 1, 1, tzinfo=timezone.utc)).total_seconds()

# Parse a web-form timestamp string into a UTC datetime object.
def parse_form_datetime(time_str):
    return datetime.strptime(time_str, "%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc)

# Serialize form data, including multi-value selections.
def serialize_form_data(form):
    # Build a JSON-compatible dictionary.
    data = {}
    for key, values in form.lists():
        data[key] = values if len(values) != 1 else values[0]
    return data

# Probe a directory in a child process so hanging mounts do not block the request.
def _probe_directory_access(path, conn):
    # Actual filesystem access in the child process.
    try:
        with os.scandir(path) as entries:
            next(entries, None)
        conn.send((True, None))
    except Exception as exc:
        conn.send((False, str(exc)))
    finally:
        conn.close()

# Test directory reachability with a short timeout.
def check_directory_access(path, timeout_sec=MET_ACCESS_TIMEOUT):
    # Communication between the parent and child process.
    parent_conn, child_conn = multiprocessing.Pipe(duplex=False)
    process = multiprocessing.Process(
        target=_probe_directory_access,
        args=(path, child_conn),
        daemon=True
    )
    process.start()
    child_conn.close()
    process.join(timeout_sec)

    if process.is_alive():
        logger.warning(f"[METCHECK] Timed out while probing meteorological data path: {path}")
        process.terminate()
        process.join(1)
        parent_conn.close()
        return False, f"Timed out after {timeout_sec} seconds."

    if parent_conn.poll():
        result = parent_conn.recv()
    else:
        result = (False, "No response from filesystem probe.")
    parent_conn.close()
    return result

# Draw a quick-look map for one output time step.
def create_plot_file(
    lons, lats, heights, filepath,
    projection='cartesian', region='global',
    central_lon=0.0, central_lat=0.0,
    lon_min=-180, lat_min=-90, lon_max=180, lat_max=90,
    z_min=-999, z_max=-999,
    met_name='', time_str='',
    mlon=None, mlat=None
):
    # Select the map projection for the plot.
    proj_map = {
        'cartesian': ccrs.PlateCarree(),
        'orthographic': ccrs.Orthographic(central_longitude=central_lon, central_latitude=central_lat),
        'robinson': ccrs.Robinson(central_longitude=central_lon)
    }
    proj = proj_map.get(projection, ccrs.PlateCarree())
    fig, ax = plt.subplots(figsize=(15, 12), subplot_kw={'projection': proj})
    if region == 'global':
        ax.set_global()
    else:
        if lon_min < lon_max and lat_min < lat_max:
            ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(color='gray', resolution='10m')
    ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='gray')
    if z_max > z_min:
        sc = ax.scatter(lons, lats, c=heights, cmap='tab20c', s=10,
                        edgecolors='none', transform=ccrs.PlateCarree(),
                        vmin=z_min, vmax=z_max)
    else:
        sc = ax.scatter(lons, lats, c=heights, cmap='tab20c', s=10,
                        edgecolors='none', transform=ccrs.PlateCarree())
    if mlon is not None and mlat is not None:
        ax.plot(mlon, mlat, 'o', color='red', markersize=10, markeredgecolor='white',
                transform=ccrs.PlateCarree(), zorder=5, label='Initial location')
    gl = ax.gridlines(draw_labels=True, color='gray', linestyle='--', alpha=0.5)
    gl.xlabel_style = gl.ylabel_style = {'color': 'black'}
    for obj in [ax, fig]:
        obj.set_facecolor('white')
    cbar = plt.colorbar(sc, orientation='horizontal', pad=0.05, aspect=50)
    cbar.set_label('Log-pressure height [km]', color='black')
    cbar.ax.xaxis.set_tick_params(color='black')
    plt.setp(cbar.ax.get_xticklabels(), color='black')
    ax.tick_params(axis='both', colors='black')
    title = f"MPTRAC | {met_name} | {time_str}"
    ax.set_title(title, fontsize=16, color='black', pad=12)
    plt.savefig(filepath, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close()

# Load a trajectory file and generate a plot image.
def process_plot(
    run_id, file, work_dir,
    map_projection, plot_region,
    central_lon, central_lat,
    mlon, mlat, met_name,
    lon_min, lat_min,
    lon_max, lat_max,
    z_min, z_max
):
    # Derive the timestamp and plot filename from the output filename.
    filename = os.path.splitext(os.path.basename(file))[0]
    parts = filename.split('_')[1:6]
    dt = datetime.strptime("_".join(parts), "%Y_%m_%d_%H_%M")
    timestamp_iso = dt.strftime("%Y-%m-%dT%H:%MZ")
    plot_filename = f"plot_{timestamp_iso}.png"
    plot_file = os.path.join(work_dir, plot_filename)
    try:
        data = np.loadtxt(file)
        if data.ndim == 1:
            data = data[np.newaxis, :]
        if data.shape[1] < 4:
            raise ValueError("Insufficient data columns")
    except Exception as e:
        logger.warning(f"[PLOT] [{run_id}] Skipping plot for {file}: {e}")
        return None
    lons, lats, heights = data[:, 2], data[:, 3], data[:, 1]
    sort_idx = np.argsort(heights)
    create_plot_file(
        lons[sort_idx], lats[sort_idx], heights[sort_idx],
        filepath=plot_file,
        projection=map_projection,
        region=plot_region,
        central_lon=central_lon,
        central_lat=central_lat,
        lon_min=lon_min,
        lat_min=lat_min,
        lon_max=lon_max,
        lat_max=lat_max,
        z_min=z_min,
        z_max=z_max,
        met_name=met_name,
        time_str=timestamp_iso,
        mlon=mlon,
        mlat=mlat
    )
    return plot_filename

# Render and log validation errors consistently.
def validation_error(message, run_id):
    logger.error(f"[VALIDATE] [{run_id}] Validation failed: {message}")
    return render_template('error.html', stdout=f"❌ {message}"), 400

# Render the main form with defaults and UI metadata.
@app.route('/')
def index():
    # Switch between normal mode and local test mode.
    test_mode = app.config.get('TEST_MODE', False)
    return render_template(
        'index.html',
        test_mode=test_mode,
        default_met_source=TEST_MET_SOURCE if test_mode else 'era5low_6h',
        default_start_time=TEST_START_TIME if test_mode else '2000-01-01 00:00',
        default_stop_time=TEST_STOP_TIME if test_mode else '2000-01-05 00:00',
        default_mz=TEST_MZ if test_mode else '0.0',
        default_uz=TEST_UZ if test_mode else '0.0',
        max_output_quantities=MAX_OUTPUT_QUANTITIES,
        default_output_quantities=DEFAULT_OUTPUT_QUANTITIES,
        output_quantity_groups=OUTPUT_QUANTITY_GROUPS,
    )

# Download the ZIP archive of a completed run.
@app.route('/download/<run_id>')
def download(run_id):
    logger.info(f"[DOWNLOAD] [{run_id}] Download requested.")
    zip_path = os.path.join(ZIPS_DIR, f"mptrac_{run_id}.zip")
    if not os.path.exists(zip_path):
        return render_template('error.html', stdout="❌ ZIP file not found."), 404
    return send_from_directory(ZIPS_DIR, f"mptrac_{run_id}.zip", as_attachment=True)

# Download the saved setup JSON of a run.
@app.route('/download_setup/<run_id>')
def download_setup(run_id):
    setup_path = os.path.join(RUNS_DIR, run_id, "setup.json")
    if not os.path.exists(setup_path):
        return render_template('error.html', stdout="❌ Setup file not found."), 404

    return send_from_directory(
        os.path.join(RUNS_DIR, run_id),
        "setup.json",
        as_attachment=True,
        download_name=f"mptrac_setup_{run_id}.json"
    )

# Main run route with validation, model execution, and result rendering.
@app.route('/run', methods=['POST'])
def run():
    
    # Unique identifier for this simulation run.
    run_id = str(uuid.uuid4())
    logger.info(f"[RUN] [{run_id}] Run started from {request.remote_addr}.")
    
    # Clean up old runs before starting a new one.
    clean_old_runs()
    
    # Parse input form...
    try:
        f = request.form
        form_data = serialize_form_data(f)
        logger.info(f"[RUN] [{run_id}] Received input form: {form_data}")
        
        # Read and convert the form parameters.
        to_f = lambda x: float(f[x])
        met_source = f.get('met_source', 'erai_6h')
        if met_source not in MET_OPTIONS:
            return validation_error(f"Unknown meteorological source: {met_source}", run_id=run_id)

        met_name = MET_OPTIONS[met_source]['MET_NAME']
        METBASE = MET_OPTIONS[met_source]['METBASE']
        DT_MET = MET_OPTIONS[met_source]['DT_MET']

        # Special overrides for local test mode.
        if app.config.get('TEST_MODE', False):
            if met_source != TEST_MET_SOURCE:
                return validation_error(
                    "Local test mode only supports ERA-Interim test data. Please select ERA-Interim.",
                    run_id=run_id
                )
            METBASE = TEST_METBASE
            DT_MET = TEST_DT_MET

        MET_PRESS_LEVEL_DEF = MET_OPTIONS[met_source]['MET_PRESS_LEVEL_DEF']
        MET_VERT_COORD = MET_OPTIONS[met_source]['MET_VERT_COORD']
        start_dt = parse_form_datetime(f['start_time'])
        start_time, stop_time = map(seconds_since_2000, (f['start_time'], f['stop_time']))
        if met_source == 'aifs_6h':
            METBASE = f"/mnt/slmet_mnt/met_data/ecmwf/open_data/data/aifs-single_{start_dt.strftime('%Y_%m_%d')}/aifs-single"
        if met_source == 'ifs_6h':
            METBASE = f"/mnt/slmet_mnt/met_data/ecmwf/open_data/data/ifs_{start_dt.strftime('%Y_%m_%d')}/ifs"
        if met_source == 'gfs_3h':
            METBASE = f"/mnt/slmet_mnt/met_data/ncep/gfs/data_{start_dt.strftime('%Y_%m_%d')}/gfs"

        # Early reachability check of the meteo mount.
        if not app.config.get('TEST_MODE', False):
            met_access_path = MET_ACCESS_PATH
            met_access_ok, met_access_error = check_directory_access(met_access_path)
            if not met_access_ok:
                return validation_error(
                    f"Meteorological data directory is not reachable: {met_access_path} ({met_access_error})",
                    run_id=run_id
                )

        mlon, mlat, mz = to_f('mlon'), to_f('mlat'), to_f('mz')
        ulon, ulat, uz = to_f('ulon'), to_f('ulat'), to_f('uz')
        slon, slat, sz = to_f('slon'), to_f('slat'), to_f('sz')
        rep = int(f['rep'])
        turb = {k: to_f(k) for k in ['turb_dx_pbl','turb_dx_trop','turb_dx_strat',
                                     'turb_dz_pbl','turb_dz_trop','turb_dz_strat',
                                     'turb_mesox','turb_mesoz']}
        conv_cape, conv_cin = to_f('conv_cape'), to_f('conv_cin')
        conv_mix_pbl = int(f.get('conv_mix_pbl', 0))
        atm_dt_out = to_f('atm_dt_out')
        output_format = f.get('output_format', 'ascii')
        selected_quantities = list(dict.fromkeys(f.getlist('quantities')))
        plot_region = f.get('plot_region', 'global')
        map_projection = f.get('map_projection', 'cartesian')
        lon_min, lat_min, lon_max, lat_max = to_f('lon_min'), to_f('lat_min'), to_f('lon_max'), to_f('lat_max')
        z_min, z_max = to_f('z_min'), to_f('z_max')
        central_lon, central_lat = to_f('central_lon'), to_f('central_lat')
        
        # Semantic validation of the input values.
        if abs(start_time - stop_time) > 30*86400:
            return validation_error("Trajectory duration exceeds 30 days.", run_id=run_id)
        if not (-100 <= mz <= 100):
            return validation_error("Mean height must be between -100 and 100 km.", run_id=run_id)
        if not (-90 <= mlat <= 90):
            return validation_error("Mean latitude must be between -90° and 90°.", run_id=run_id)
        if not (-180 <= mlon <= 180):
            return validation_error("Mean longitude must be between -180° and 180°.", run_id=run_id)
        if not (1 <= rep <= 10000):
            return validation_error("Number of particles must be between 1 and 10,000.", run_id=run_id)
        if any(p < 0 for p in [ulon, ulat, uz, slon, slat, sz] + list(turb.values())):
            return validation_error("Turbulence parameters must be zero or positive.", run_id=run_id)
        if not (-180 <= central_lon <= 180):
            return validation_error("Central longitude must be between -180° and 180°.", run_id=run_id)
        if not (-90 <= central_lat <= 90):
            return validation_error("Central latitude must be between -90° and 90°.", run_id=run_id)
        if output_format not in ['ascii', 'netcdf']:
            return validation_error("Output format must be ASCII or netCDF.", run_id=run_id)
        if len(selected_quantities) > MAX_OUTPUT_QUANTITIES:
            return validation_error(
                f"Select at most {MAX_OUTPUT_QUANTITIES} output quantities.",
                run_id=run_id
            )
        invalid_quantities = [q for q in selected_quantities if q not in ALLOWED_OUTPUT_QUANTITIES]
        if invalid_quantities:
            return validation_error(
                "Unknown output quantities: " + ", ".join(invalid_quantities),
                run_id=run_id
            )
        
        # Set forward/backward trajectory flag...
        direction = -1 if stop_time < start_time else 1

    except Exception as e:
        return validation_error(str(e), run_id=run_id)
    
    # Create the run-specific working directory.
    work_dir = os.path.join(RUNS_DIR, run_id)
    os.makedirs(work_dir, exist_ok=True)

    # Save the complete form setup as JSON.
    setup_file = os.path.join(work_dir, "setup.json")
    with open(setup_file, "w") as fobj:
        json.dump(form_data, fobj, indent=2)
    
    # Create the MPTRAC control file.
    ctl_file, dirlist_file, init_file = map(lambda x: os.path.join(work_dir, x), ['trac.ctl', 'dirlist.txt', 'atm_init.tab'])
    atm_file = os.path.join(work_dir, 'atm')
    with open(dirlist_file, 'w') as f: f.write("./\n")

    # Build NQ and QNT_NAME entries dynamically.
    quantity_block = "\n".join(
        ['# Variables...', f'NQ = {len(selected_quantities)}']
        + [f'QNT_NAME[{idx}] = {name}' for idx, name in enumerate(selected_quantities)]
    )

    ctl_template = quantity_block + textwrap.dedent(f"""

    # atm_init...
    INIT_T0 = {start_time}
    INIT_T1 = {start_time}
    INIT_Z0 = {mz}
    INIT_Z1 = {mz}
    INIT_LON0 = {mlon}
    INIT_LON1 = {mlon}
    INIT_LAT0 = {mlat}
    INIT_LAT1 = {mlat}
    INIT_ULON = {ulon}
    INIT_ULAT = {ulat}
    INIT_UZ = {uz}
    INIT_SLON = {slon}
    INIT_SLAT = {slat}
    INIT_SZ = {sz}
    INIT_REP = {rep}

    # trac...
    DIRECTION = {direction}
    T_STOP = {stop_time}
    METBASE = {METBASE}
    DT_MET = {DT_MET}
    MET_PRESS_LEVEL_DEF = {MET_PRESS_LEVEL_DEF}
    MET_VERT_COORD = {MET_VERT_COORD}
    DIFFUSION = 1
    """ + '\n'.join(f"{k.upper()} = {v}" for k, v in turb.items()) + f"""
    CONV_CAPE = {conv_cape}
    CONV_CIN = {conv_cin}
    CONV_MIX_PBL = {conv_mix_pbl}
    ATM_BASENAME = {atm_file}
    ATM_DT_OUT = {atm_dt_out}
    """)
    if met_source in ['era5low_6h']:
        ctl_template += "MET_CLAMS = 1\n"
    elif met_source in ['jra55_6h']:
        ctl_template += """
MET_NLEV = 61
MET_LEV_HYAM[0] = 0.00000000000000E+00
MET_LEV_HYAM[1] = 0.00000000000000E+00
MET_LEV_HYAM[2] = 0.00000000000000E+00
MET_LEV_HYAM[3] = 0.00000000000000E+00
MET_LEV_HYAM[4] = 0.00000000000000E+00
MET_LEV_HYAM[5] = 0.00000000000000E+00
MET_LEV_HYAM[6] = 0.00000000000000E+00
MET_LEV_HYAM[7] = 0.00000000000000E+00
MET_LEV_HYAM[8] = 1.33051011276943E+02
MET_LEV_HYAM[9] = 3.64904148871589E+02
MET_LEV_HYAM[10] = 6.34602716447362E+02
MET_LEV_HYAM[11] = 9.59797167291774E+02
MET_LEV_HYAM[12] = 1.34768004165515E+03
MET_LEV_HYAM[13] = 1.79090739595110E+03
MET_LEV_HYAM[14] = 2.29484168994850E+03
MET_LEV_HYAM[15] = 2.84748477771176E+03
MET_LEV_HYAM[16] = 3.46887148811864E+03
MET_LEV_HYAM[17] = 4.16295646296916E+03
MET_LEV_HYAM[18] = 4.89188083250491E+03
MET_LEV_HYAM[19] = 5.67182423980408E+03
MET_LEV_HYAM[20] = 6.47671299638532E+03
MET_LEV_HYAM[21] = 7.29746989472049E+03
MET_LEV_HYAM[22] = 8.12215979124915E+03
MET_LEV_HYAM[23] = 8.91408220106234E+03
MET_LEV_HYAM[24] = 9.65618191050164E+03
MET_LEV_HYAM[25] = 1.03294361777746E+04
MET_LEV_HYAM[26] = 1.09126384442387E+04
MET_LEV_HYAM[27] = 1.13696478308432E+04
MET_LEV_HYAM[28] = 1.16953715974700E+04
MET_LEV_HYAM[29] = 1.18612530873948E+04
MET_LEV_HYAM[30] = 1.18554343163493E+04
MET_LEV_HYAM[31] = 1.16633553655803E+04
MET_LEV_HYAM[32] = 1.12854040644942E+04
MET_LEV_HYAM[33] = 1.07299494055679E+04
MET_LEV_HYAM[34] = 1.00146150535107E+04
MET_LEV_HYAM[35] = 9.16724703583310E+03
MET_LEV_HYAM[36] = 8.22624490770442E+03
MET_LEV_HYAM[37] = 7.20156898029828E+03
MET_LEV_HYAM[38] = 6.08867300853392E+03
MET_LEV_HYAM[39] = 4.95000000000000E+03
MET_LEV_HYAM[40] = 4.00000000000000E+03
MET_LEV_HYAM[41] = 3.23000000000000E+03
MET_LEV_HYAM[42] = 2.61000000000000E+03
MET_LEV_HYAM[43] = 2.10500000000000E+03
MET_LEV_HYAM[44] = 1.70000000000000E+03
MET_LEV_HYAM[45] = 1.37000000000000E+03
MET_LEV_HYAM[46] = 1.10500000000000E+03
MET_LEV_HYAM[47] = 8.93000000000000E+02
MET_LEV_HYAM[48] = 7.20000000000000E+02
MET_LEV_HYAM[49] = 5.81000000000000E+02
MET_LEV_HYAM[50] = 4.69000000000000E+02
MET_LEV_HYAM[51] = 3.77000000000000E+02
MET_LEV_HYAM[52] = 3.01000000000000E+02
MET_LEV_HYAM[53] = 2.37000000000000E+02
MET_LEV_HYAM[54] = 1.82000000000000E+02
MET_LEV_HYAM[55] = 1.36000000000000E+02
MET_LEV_HYAM[56] = 9.70000000000000E+01
MET_LEV_HYAM[57] = 6.50000000000000E+01
MET_LEV_HYAM[58] = 3.90000000000000E+01
MET_LEV_HYAM[59] = 2.00000000000000E+01
MET_LEV_HYAM[60] = 0.00000000000000E+00
MET_LEV_HYBM[0] = 1.00000000000000E+00
MET_LEV_HYBM[1] = 9.97000000000000E-01
MET_LEV_HYBM[2] = 9.94000000000000E-01
MET_LEV_HYBM[3] = 9.89000000000000E-01
MET_LEV_HYBM[4] = 9.82000000000000E-01
MET_LEV_HYBM[5] = 9.72000000000000E-01
MET_LEV_HYBM[6] = 9.60000000000000E-01
MET_LEV_HYBM[7] = 9.46000000000000E-01
MET_LEV_HYBM[8] = 9.26669489887231E-01
MET_LEV_HYBM[9] = 9.04350958511284E-01
MET_LEV_HYBM[10] = 8.79653972835526E-01
MET_LEV_HYBM[11] = 8.51402028327082E-01
MET_LEV_HYBM[12] = 8.19523199583449E-01
MET_LEV_HYBM[13] = 7.85090926040489E-01
MET_LEV_HYBM[14] = 7.48051583100515E-01
MET_LEV_HYBM[15] = 7.09525152222882E-01
MET_LEV_HYBM[16] = 6.68311285118814E-01
MET_LEV_HYBM[17] = 6.24370435370308E-01
MET_LEV_HYBM[18] = 5.80081191674951E-01
MET_LEV_HYBM[19] = 5.34281757601959E-01
MET_LEV_HYBM[20] = 4.88232870036147E-01
MET_LEV_HYBM[21] = 4.42025301052795E-01
MET_LEV_HYBM[22] = 3.95778402087509E-01
MET_LEV_HYBM[23] = 3.50859177989377E-01
MET_LEV_HYBM[24] = 3.07438180894984E-01
MET_LEV_HYBM[25] = 2.65705638222254E-01
MET_LEV_HYBM[26] = 2.25873615557613E-01
MET_LEV_HYBM[27] = 1.89303521691568E-01
MET_LEV_HYBM[28] = 1.55046284025300E-01
MET_LEV_HYBM[29] = 1.24387469126052E-01
MET_LEV_HYBM[30] = 9.64456568365075E-02
MET_LEV_HYBM[31] = 7.23664463441966E-02
MET_LEV_HYBM[32] = 5.21459593550578E-02
MET_LEV_HYBM[33] = 3.57005059443214E-02
MET_LEV_HYBM[34] = 2.28538494648935E-02
MET_LEV_HYBM[35] = 1.33275296416689E-02
MET_LEV_HYBM[36] = 6.73755092295582E-03
MET_LEV_HYBM[37] = 2.48431019701722E-03
MET_LEV_HYBM[38] = 1.13269914660783E-04
MET_LEV_HYBM[39] = 0.00000000000000E+00
MET_LEV_HYBM[40] = 0.00000000000000E+00
MET_LEV_HYBM[41] = 0.00000000000000E+00
MET_LEV_HYBM[42] = 0.00000000000000E+00
MET_LEV_HYBM[43] = 0.00000000000000E+00
MET_LEV_HYBM[44] = 0.00000000000000E+00
MET_LEV_HYBM[45] = 0.00000000000000E+00
MET_LEV_HYBM[46] = 0.00000000000000E+00
MET_LEV_HYBM[47] = 0.00000000000000E+00
MET_LEV_HYBM[48] = 0.00000000000000E+00
MET_LEV_HYBM[49] = 0.00000000000000E+00
MET_LEV_HYBM[50] = 0.00000000000000E+00
MET_LEV_HYBM[51] = 0.00000000000000E+00
MET_LEV_HYBM[52] = 0.00000000000000E+00
MET_LEV_HYBM[53] = 0.00000000000000E+00
MET_LEV_HYBM[54] = 0.00000000000000E+00
MET_LEV_HYBM[55] = 0.00000000000000E+00
MET_LEV_HYBM[56] = 0.00000000000000E+00
MET_LEV_HYBM[57] = 0.00000000000000E+00
MET_LEV_HYBM[58] = 0.00000000000000E+00
MET_LEV_HYBM[59] = 0.00000000000000E+00
MET_LEV_HYBM[60] = 0.00000000000000E+00
"""
    elif met_source in ['merra2_3h']:
        ctl_template += "MET_NC_SCALE = 0\n"
    elif met_source in ['ncep2_6h']:
        ctl_template += "MET_RELHUM = 1\n"
    
    # Write the final control file into the working directory.
    with open(ctl_file, 'w') as f: f.write(ctl_template)
    logger.info(f"[RUN] [{run_id}] Control file written at {ctl_file}.")
    
    # Launch the MPTRAC programs with timeouts.
    atm_init_code, atm_init_output = run_command([ATM_INIT_CMD, ctl_file, init_file], timeout=120, run_id=run_id)
    trac_code, trac_output = run_command([TRAC_CMD, dirlist_file, ctl_file, init_file], timeout=600, run_id=run_id)
    
    # Summary logging for success or failure.
    if atm_init_code != 0 or trac_code != 0:
        logger.error(f"[RUN] [{run_id}] Run failed: atm_init_code={atm_init_code}, trac_code={trac_code}")
    else:
        logger.info(f"[RUN] [{run_id}] Run succeeded.")

    # Combine tool output for the result page.
    combined_output = f"=== atm_init Output ===\n{atm_init_output}\n\n=== trac Output ===\n{trac_output}"
    
    # Follow-up workflow after a successful trac run.
    if trac_code == 0:

        # Initialize the list of generated plot files.
        plot_filenames = []
        
        # Collect all trajectory output files.
        files = sorted(glob.glob(os.path.join(work_dir, 'atm_19*.tab')) +
                       glob.glob(os.path.join(work_dir, 'atm_20*.tab')))

        # Generate the quick-look plots in parallel.
        logger.info(f"[PLOT] [{run_id}] {len(files)} plot files to process.")
        with ProcessPoolExecutor(max_workers = min(8, os.cpu_count())) as executor:
            futures = [
                executor.submit(process_plot, run_id, file, work_dir,
                                map_projection, plot_region, central_lon, central_lat,
                                mlon, mlat, met_name, lon_min, lat_min, lon_max, lat_max, z_min, z_max)
                for file in files
            ]
            for future in as_completed(futures):
                result = future.result()
                if result:
                    plot_filenames.append(result)
                    
        # Sort the plot filenames by time.
        plot_filenames = sorted(
            [f for f in plot_filenames if f],
            key=lambda name: datetime.strptime(name.split("_")[1][:-4], "%Y-%m-%dT%H:%MZ")
        )
        logger.info(f"[PLOT] [{run_id}] Plots finished.")

        # Optional netCDF conversion after plotting.
        atm_conv_output = ""
        if output_format == 'netcdf':
            conv_logs = []
            for ascii_file in files:
                nc_file = os.path.splitext(ascii_file)[0] + '.nc'
                conv_code, conv_output = run_command(
                    [ATM_CONV_CMD, ctl_file, ascii_file, '0', nc_file, '2'],
                    timeout=120,
                    run_id=run_id
                )
                conv_logs.append(
                    f"=== atm_conv Output | {os.path.basename(ascii_file)} ===\n{conv_output}"
                )
                if conv_code != 0:
                    combined_output += "\n\n" + "\n\n".join(conv_logs)
                    logger.error(
                        f"[CONVERT] [{run_id}] netCDF conversion failed for {ascii_file}: exit code {conv_code}"
                    )
                    return render_template('error.html', stdout=combined_output), 500
                try:
                    os.remove(ascii_file)
                except Exception as e:
                    combined_output += "\n\n" + "\n\n".join(conv_logs)
                    logger.error(f"[CONVERT] [{run_id}] Failed to remove {ascii_file}: {e}")
                    return render_template('error.html', stdout=combined_output), 500
            atm_conv_output = "\n\n" + "\n\n".join(conv_logs)
            combined_output += atm_conv_output
            logger.info(f"[CONVERT] [{run_id}] Converted {len(files)} atmosphere files to netCDF.")
        
        # Package all run results into a ZIP archive.
        zip_path = os.path.join(ZIPS_DIR, f"mptrac_{run_id}.zip")
        logger.info(f"[ZIP] [{run_id}] Creating zip file: {zip_path}")
        fast_zip_no_compression(zip_path, work_dir)
        logger.info(f"[ZIP] [{run_id}] Zip file created successfully: {zip_path}")
        
        # Return the result page to the browser.
        return render_template(
            'result.html',
            run_id=run_id,
            stdout=combined_output,
            plot_filenames=plot_filenames
        )
    
    # Error page for a failed model run.
    return render_template('error.html', stdout=combined_output), 500

# Render the terms and conditions page.
@app.route('/terms')
def terms():
    return render_template('terms.html')

# Serve generated plot images from the run directory.
@app.route('/runs/working/<run_id>/<filename>')
def serve_plot_image(run_id, filename):
    return send_from_directory(os.path.join(RUNS_DIR, run_id), filename)

# Show the log file for browser-based debugging.
@app.route('/logs')
def show_logs():
    # Read the current log file for the browser view.
    try:
        with open(log_file, 'r') as f:
            log_contents = f.read()
        return f"<pre>{log_contents}</pre>"
    except Exception as e:
        logger.error(f"[LOG] Error reading log file: {e}")
        return render_template('error.html', stdout=f"❌ Error reading log file: {e}"), 500

# Parse command-line options and start the Flask server.
if __name__ == '__main__':
    # Define the optional command-line arguments.
    parser = argparse.ArgumentParser(description='Run the MPTRAC web runner.')
    parser.add_argument(
        '--test',
        action='store_true',
        help='Use the local ERA-Interim test dataset from ../../tests/data/ei.'
    )
    args = parser.parse_args()

    app.config['TEST_MODE'] = args.test
    if app.config['TEST_MODE']:
        logger.info('[CONFIG] Running in local test mode with ../../tests/data/ei.')

    app.run(debug=False)
