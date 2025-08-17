from flask import Flask, render_template, request, send_from_directory
import os, time, uuid, glob, shutil, textwrap, zipfile, subprocess, threading
from datetime import datetime, timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import logging
from logging.handlers import RotatingFileHandler
from concurrent.futures import ProcessPoolExecutor, as_completed

# Config...
ATM_INIT_CMD = '../../src/atm_init'
TRAC_CMD = '../../src/trac'

# Create directories...
RUNS_DIR, ZIPS_DIR, LOG_DIR = 'runs/working', 'runs/zips', 'runs/logs'
for path in [RUNS_DIR, ZIPS_DIR, LOG_DIR]:
    os.makedirs(path, exist_ok=True)

# Semaphore...
process_semaphore = threading.Semaphore(3)

# Logging setup...
log_file = os.path.join(LOG_DIR, 'web_runner.log')
logger = logging.getLogger('web_runner')
logger.setLevel(logging.INFO)
handler = RotatingFileHandler(log_file, maxBytes=10*1024*1024, backupCount=5)
handler.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s in %(threadName)s: %(message)s'))
logger.addHandler(handler)

# Clean up working directory...
def clean_old_runs(max_age_sec=3600):
    now = time.time()
    for directory in [RUNS_DIR, ZIPS_DIR]:
        logger.info(f"[CLEANUP] Scanning {directory} for files older than {max_age_sec} seconds.")
        for path in glob.glob(os.path.join(directory, "*")):
            if now - os.path.getmtime(path) > max_age_sec:
                try:
                    logger.info(f"[CLEANUP] Removed old file or directory: {path}")
                    shutil.rmtree(path) if os.path.isdir(path) else os.remove(path)
                except Exception as e:
                    logger.warning(f"[CLEANUP] Failed to remove {path}: {e}")

# Periodic cleaning...
def schedule_periodic_cleanup(interval=3600):
    threading.Thread(target=lambda: [clean_old_runs() or time.sleep(interval) for _ in iter(int, 1)],
                     daemon=True).start()

# Create zip file...
def fast_zip_no_compression(zip_path, source_dir):
    with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_STORED) as zipf:
        for root, _, files in os.walk(source_dir):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, source_dir)
                zipf.write(file_path, arcname)

# Flask app...
app = Flask(__name__)

# Meteo options...
MET_OPTIONS = {
    'aifs_6h': dict(METBASE='/mnt/slmet-mnt/met_data/ecmwf/open_data/data/aifs-single_YYYY_MM_DD/aifs-single',
                    DT_MET=21600, MET_VERT_COORD=0, MET_PRESS_LEVEL_DEF=-1, MET_NAME='ECMWF AIFS'),
    'ifs_6h': dict(METBASE='/mnt/slmet-mnt/met_data/ecmwf/open_data/data/ifs_YYYY_MM_DD/ifs',
                   DT_MET=21600, MET_VERT_COORD=0, MET_PRESS_LEVEL_DEF=-1, MET_NAME='ECMWF IFS'),
    'gfs_3h': dict(METBASE='/mnt/slmet-mnt/met_data/ncep/gfs/data_YYYY_MM_DD/gfs',
                   DT_MET=10800, MET_VERT_COORD=0, MET_PRESS_LEVEL_DEF=-1, MET_NAME='NCEP GFS'),
    'era5low_6h': dict(METBASE='/mnt/slmet-mnt/met_data/ecmwf/era5.1/resolution_1x1/nc/YYYY/MM/era5',
                       DT_MET=21600, MET_VERT_COORD=1, MET_PRESS_LEVEL_DEF=6, MET_NAME='ERA5'),
    'erai_6h': dict(METBASE='/mnt/slmet-mnt/met_data/ecmwf/era_interim/pressure_0.75deg_v2/nc/YYYY/ei',
                    DT_MET=21600, MET_VERT_COORD=0, MET_PRESS_LEVEL_DEF=-1, MET_NAME='ERA-Interim'),
    'jra3q_6h': dict(METBASE='/mnt/slmet-mnt/met_data/jma/jra3q/nc/YYYY/jra3q',
                     DT_MET=21600, MET_VERT_COORD=2, MET_PRESS_LEVEL_DEF=6, MET_NAME='JRA-3Q'),
    'jra55_6h': dict(METBASE='/mnt/slmet-mnt/met_data/jma/jra55/nc/YYYY/jra55',
                     DT_MET=21600, MET_VERT_COORD=4, MET_PRESS_LEVEL_DEF=6, MET_NAME='JRA-55'),
    'merra2_3h': dict(METBASE='/mnt/slmet-mnt/met_data/nasa/merra-2/hybrid/YYYY/merra2',
                      DT_MET=10800, MET_VERT_COORD=1, MET_PRESS_LEVEL_DEF=6, MET_NAME='MERRA-2'),
    'ncep_6h': dict(METBASE='/mnt/slmet-mnt/met_data/ncep/reanalysis/nc/YYYY/ncep',
                    DT_MET=21600, MET_VERT_COORD=0, MET_PRESS_LEVEL_DEF=-1, MET_NAME='NCEP-NCAR Reanalysis 1'),
    'ncep2_6h': dict(METBASE='/mnt/slmet-mnt/met_data/ncep/reanalysis2/nc/YYYY/ncep2',
                     DT_MET=21600, MET_VERT_COORD=0, MET_PRESS_LEVEL_DEF=-1, MET_NAME='NCEP-DOE Reanalysis 2')
}

# Execute shell command and check result...
def run_command(cmd, timeout, run_id):
    if not process_semaphore.acquire(timeout=5):
        logger.error(f"[EXEC] [{run_id}] Server is too busy.")
        return -999, "❌ Server is too busy."
    try:
        logger.info(f"[EXEC] [{run_id}] Executing command: {' '.join(cmd)}")
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=timeout)
        logger.info(f"[EXEC] [{run_id}] Command finished with exit code {result.returncode}.")
        return result.returncode, result.stdout
    except subprocess.TimeoutExpired as e:
        logger.warning(f"[EXEC] [{run_id}] Command timed out after {timeout} seconds.")
        return -1, (e.stdout or '') + f"\n❌ Command timed out after {timeout} seconds."
    finally:
        process_semaphore.release()

# Calculate seconds since 2000-01-01, 00:00 UTC...
def seconds_since_2000(time_str):
    dt = datetime.strptime(time_str, "%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc)
    return (dt - datetime(2000, 1, 1, tzinfo=timezone.utc)).total_seconds()

# Create map plot of air parcel data...
def create_plot_file(
    lons, lats, heights, filepath,
    projection='cartesian', region='global',
    central_lon=0.0, central_lat=0.0,
    lon_min=-180, lat_min=-90, lon_max=180, lat_max=90,
    z_min=-999, z_max=-999,
    met_name='', time_str='',
    mlon=None, mlat=None
):
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

# Parallel processing of plots...
def process_plot(
    run_id, file, work_dir,
    map_projection, plot_region,
    central_lon, central_lat,
    mlon, mlat, met_name,
    lon_min, lat_min,
    lon_max, lat_max,
    z_min, z_max
):
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

# Exit with error message...
def validation_error(message, run_id):
    logger.error(f"[VALIDATE] [{run_id}] Validation failed: {message}")
    return render_template('error.html', stdout=f"❌ {message}"), 400

@app.route('/')
def index(): return render_template('index.html')

@app.route('/download/<run_id>')
def download(run_id):
    logger.info(f"[DOWNLOAD] [{run_id}] Download requested.")
    return send_from_directory(ZIPS_DIR, f"{run_id}.zip", as_attachment=True)

@app.route('/run', methods=['POST'])
def run():
    
    # Get run ID...
    run_id = str(uuid.uuid4())
    logger.info(f"[RUN] [{run_id}] Run started from {request.remote_addr}.")
    
    # Parse input form...
    try:
        f = request.form
        logger.info(f"[RUN] [{run_id}] Received input form: {dict(f)}")
        
        # Parse parameters...
        to_f = lambda x: float(f[x])
        met_source = f.get('met_source', 'erai_6h')
        met_name = MET_OPTIONS[met_source]['MET_NAME']
        METBASE = MET_OPTIONS[met_source]['METBASE']
        DT_MET = MET_OPTIONS[met_source]['DT_MET']
        # METBASE = '../../tests/data/ei'
        # DT_MET = 86400
        MET_PRESS_LEVEL_DEF = MET_OPTIONS[met_source]['MET_PRESS_LEVEL_DEF']
        MET_VERT_COORD = MET_OPTIONS[met_source]['MET_VERT_COORD']
        start_dt = datetime.strptime(f['start_time'], "%Y-%m-%d %H:%M")
        start_time, stop_time = map(seconds_since_2000, (f['start_time'], f['stop_time']))
        if met_source == 'aifs_6h':
            METBASE = f"/mnt/slmet-mnt/met_data/ecmwf/open_data/data/aifs-single_{start_dt.strftime('%Y_%m_%d')}/aifs-single"
        if met_source == 'ifs_6h':
            METBASE = f"/mnt/slmet-mnt/met_data/ecmwf/open_data/data/ifs_{start_dt.strftime('%Y_%m_%d')}/ifs"
        if met_source == 'gfs_3h':
            METBASE = f"/mnt/slmet-mnt/met_data/ncep/gfs/data_{start_dt.strftime('%Y_%m_%d')}/gfs"
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
        plot_region = f.get('plot_region', 'global')
        map_projection = f.get('map_projection', 'cartesian')
        lon_min, lat_min, lon_max, lat_max = to_f('lon_min'), to_f('lat_min'), to_f('lon_max'), to_f('lat_max')
        z_min, z_max = to_f('z_min'), to_f('z_max')
        central_lon, central_lat = to_f('central_lon'), to_f('central_lat')
        
        # Check values...
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
        
        # Set forward/backward trajectory flag...
        direction = -1 if stop_time < start_time else 1

    except Exception as e:
        return validation_error(str(e), run_id=run_id)
    
    # Create working directory...
    work_dir = os.path.join(RUNS_DIR, run_id)
    os.makedirs(work_dir, exist_ok=True)
    
    # Create control parameter file...
    ctl_file, dirlist_file, init_file = map(lambda x: os.path.join(work_dir, x), ['trac.ctl', 'dirlist.txt', 'atm_init.tab'])
    atm_file = os.path.join(work_dir, 'atm')
    with open(dirlist_file, 'w') as f: f.write("./\n")

    ctl_template = textwrap.dedent(f"""\
    # Variables...
    NQ = 14
    QNT_NAME[0] = zg
    QNT_NAME[1] = p
    QNT_NAME[2] = t
    QNT_NAME[3] = theta
    QNT_NAME[4] = u
    QNT_NAME[5] = v
    QNT_NAME[6] = w
    QNT_NAME[7] = pv
    QNT_NAME[8] = h2o
    QNT_NAME[9] = o3
    QNT_NAME[10] = cc
    QNT_NAME[11] = ps
    QNT_NAME[12] = pbl
    QNT_NAME[13] = pt
    
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
    
    with open(ctl_file, 'w') as f: f.write(ctl_template)
    logger.info(f"[RUN] [{run_id}] Control file written at {ctl_file}.")
    
    # Run atm_init and trac...
    atm_init_code, atm_init_output = run_command([ATM_INIT_CMD, ctl_file, init_file], timeout=120, run_id=run_id)
    trac_code, trac_output = run_command([TRAC_CMD, dirlist_file, ctl_file, init_file], timeout=600, run_id=run_id)
    
    # Catch errors...
    if trac_code == -999:
        logger.error(f"[RUN] [{run_id}] Server busy.")
        return render_template('server_busy.html'), 503
    
    # Logging...
    if atm_init_code != 0 or trac_code != 0:
        logger.error(f"[RUN] [{run_id}] Run failed: atm_init_code={atm_init_code}, trac_code={trac_code}")
    else:
        logger.info(f"[RUN] [{run_id}] Run succeeded.")

    # Combine log message...
    combined_output = f"=== atm_init Output ===\n{atm_init_output}\n\n=== trac Output ===\n{trac_output}"
    
    # Run successfull...
    if trac_code == 0:

        # Init...
        plot_filenames = []
        
        # Get list of filenames...
        files = sorted(glob.glob(os.path.join(work_dir, 'atm_19*.tab')) +
                       glob.glob(os.path.join(work_dir, 'atm_20*.tab')))

        # Parallel creation of plots...
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
                    
        # Sort filenames...
        plot_filenames = sorted(
            [f for f in plot_filenames if f],
            key=lambda name: datetime.strptime(name.split("_")[1][:-4], "%Y-%m-%dT%H:%MZ")
        )
        logger.info(f"[PLOT] [{run_id}] Plots finished.")
        
        # Create zip file...
        zip_path = os.path.join(ZIPS_DIR, f"{run_id}.zip")
        logger.info(f"[ZIP] [{run_id}] Creating zip file: {zip_path}")
        fast_zip_no_compression(zip_path, work_dir)
        logger.info(f"[ZIP] [{run_id}] Zip file created successfully: {zip_path}")
        
        # Show results...
        return render_template(
            'result.html',
            run_id=run_id,
            stdout=combined_output,
            plot_filenames=plot_filenames
        )
    
    # Show error...
    return render_template('error.html', stdout=combined_output), 500

@app.route('/terms')
def terms():
    return render_template('terms.html')

@app.route('/runs/working/<run_id>/<filename>')
def serve_plot_image(run_id, filename):
    return send_from_directory(os.path.join(RUNS_DIR, run_id), filename)

@app.route('/logs')
def show_logs():
    try:
        with open(log_file, 'r') as f:
            log_contents = f.read()
        return f"<pre>{log_contents}</pre>"
    except Exception as e:
        logger.error(f"[LOG] Error reading log file: {e}")
        return render_template('error.html', stdout=f"❌ Error reading log file: {e}"), 500

if __name__ == '__main__':
    clean_old_runs()
    schedule_periodic_cleanup()
    app.run(debug=False)
