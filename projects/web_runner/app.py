from flask import Flask, render_template, request, send_from_directory
import subprocess, threading, os, uuid, shutil, time, io, base64, glob, textwrap
import numpy as np
from datetime import datetime, timezone
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Config...
ATM_INIT_CMD = '../../src/atm_init'
TRAC_CMD = '../../src/trac'
RUNS_DIR, ZIPS_DIR = 'runs/working', 'runs/zips'
os.makedirs(RUNS_DIR, exist_ok=True)
os.makedirs(ZIPS_DIR, exist_ok=True)
process_semaphore = threading.Semaphore(3)
app = Flask(__name__)

# Meteo options...
MET_OPTIONS = {
    'era5low_6h': {
        'METBASE': '/mnt/slmet-mnt/met_data/ecmwf/era5.1/resolution_1x1/nc/YYYY/MM/era5',
        'DT_MET': 21600.0,
        'MET_PRESS_LEVEL_DEF': 6
    },
    'erai_6h': {
        'METBASE': '/mnt/slmet-mnt/met_data/ecmwf/era_interim/pressure_0.75deg_v2/nc/YYYY/ei',
        'DT_MET': 21600.0,
        'MET_PRESS_LEVEL_DEF': -1
    },
    'merra2_3h': {
        'METBASE': '/mnt/slmet-mnt/met_data/nasa/merra-2/hybrid/YYYY/merra2',
        'DT_MET': 10800.0,
        'MET_PRESS_LEVEL_DEF': 6
    },
    'merra2_6h': {
        'METBASE': '/mnt/slmet-mnt/met_data/nasa/merra-2/hybrid/YYYY/merra2',
        'DT_MET': 21600.0,
        'MET_PRESS_LEVEL_DEF': 6
    },
    'ncep_6h': {
        'METBASE': '/mnt/slmet-mnt/met_data/ncep/reanalysis/nc/YYYY/ncep',
        'DT_MET': 21600.0,
        'MET_PRESS_LEVEL_DEF': -1
    },
    'ncep2_6h': {
        'METBASE': '/mnt/slmet-mnt/met_data/ncep/reanalysis2/nc/YYYY/ncep2',
        'DT_MET': 21600.0,
        'MET_PRESS_LEVEL_DEF': -1
    }
}

# Clean up working directory...
def delayed_cleanup(path, delay=600):
    threading.Thread(target=lambda: (time.sleep(delay), shutil.rmtree(path, ignore_errors=True)), daemon=True).start()

# Execute shell command and check result...
def run_command(cmd, timeout=300):
    if not process_semaphore.acquire(timeout=5):
        return -999, "❌ Server is too busy."
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=timeout)
        return result.returncode, result.stdout
    except subprocess.TimeoutExpired as e:
        return -1, (e.stdout or '') + f"\n❌ Command timed out after {timeout} seconds."
    finally:
        process_semaphore.release()

# Calculate seconds since 2000-01-01, 00:00 UTC...
def seconds_since_2000(time_str):
    dt = datetime.strptime(time_str, "%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc)
    return (dt - datetime(2000, 1, 1, tzinfo=timezone.utc)).total_seconds()

# Create map plot of air parcel data...
def create_plot(lons, lats, heights, projection='cartesian', region='global', central_lon=0.0, central_lat=0.0):
    proj_map = {
        'cartesian': ccrs.PlateCarree(),
        'orthographic': ccrs.Orthographic(central_longitude=central_lon, central_latitude=central_lat),
        'robinson': ccrs.Robinson(central_longitude=central_lon)
    }
    proj = proj_map.get(projection, ccrs.PlateCarree())
    fig, ax = plt.subplots(figsize=(15, 12), subplot_kw={'projection': proj})
    if region == 'global':
        ax.set_global()
    ax.coastlines(color='gray', resolution='10m')
    ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='gray')
    sc = ax.scatter(lons, lats, c=heights, cmap='tab20c', s=10, edgecolors='none', transform=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, color='gray', linestyle='--', alpha=0.5)
    gl.xlabel_style = gl.ylabel_style = {'color': 'black'}
    for obj in [ax, fig]:
        obj.set_facecolor('white')
    cbar = plt.colorbar(sc, orientation='horizontal', pad=0.05, aspect=50)
    cbar.set_label('Log-pressure height [km]', color='black')
    cbar.ax.xaxis.set_tick_params(color='black')
    plt.setp(cbar.ax.get_xticklabels(), color='black')
    ax.tick_params(axis='both', colors='black')
    ax.set_title('MPTRAC | Air Parcel Trajectories', fontsize=16, color='black')
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close()
    return base64.b64encode(buf.getvalue()).decode()

# Exit with error message...
def validation_error(message):
    return render_template('error.html', stdout=f"❌ {message}"), 400

@app.route('/')
def index(): return render_template('index.html')

@app.route('/download/<run_id>')
def download(run_id):
    return send_from_directory(ZIPS_DIR, f"{run_id}.zip", as_attachment=True)

@app.route('/run', methods=['POST'])
def run():
    try:
        # Get control parameters from input form...
        f = request.form
        to_f = lambda x: float(f[x])
        met_source = f.get('met_source', 'erai_6h')
        METBASE = MET_OPTIONS[met_source]['METBASE']
        DT_MET = MET_OPTIONS[met_source]['DT_MET']
        #METBASE = '../../tests/data/ei'
        #DT_MET = 86400
        MET_PRESS_LEVEL_DEF = MET_OPTIONS[met_source]['MET_PRESS_LEVEL_DEF']
        start_time, stop_time = map(seconds_since_2000, (f['start_time'], f['stop_time']))
        z0, z1, dz = to_f('z0'), to_f('z1'), to_f('dz')
        lat0, lat1, dlat = to_f('lat0'), to_f('lat1'), to_f('dlat')
        lon0, lon1, dlon = to_f('lon0'), to_f('lon1'), to_f('dlon')
        rep = int(f['rep'])
        ulon, ulat, uz = to_f('ulon'), to_f('ulat'), to_f('uz')
        slon, slat, sz = to_f('slon'), to_f('slat'), to_f('sz')
        turb = {k: to_f(k) for k in ['turb_dx_pbl','turb_dx_trop','turb_dx_strat','turb_dz_pbl','turb_dz_trop','turb_dz_strat','turb_mesox','turb_mesoz']}
        atm_dt_out = to_f('atm_dt_out')
        plot_region = f.get('plot_region', 'global')
        map_projection = f.get('map_projection', 'cartesian')
        central_lon = to_f('central_lon')
        central_lat = to_f('central_lat')

        # Check values...
        if abs(start_time - stop_time) > 30*86400: return validation_error("Duration exceeds 30 days.")
        if not (-100 <= z0 <= 100 and -100 <= z1 <= 100): return validation_error("Height range invalid.")
        if not (0 < dz <= 100): return validation_error("dz must be > 0 and ≤ 100 km.")
        if not (-90 <= lat0 <= 90 and -90 <= lat1 <= 90): return validation_error("Latitude range invalid.")
        if not (0 < dlat <= 180): return validation_error("Latitude sampling invalid.")
        if not (-180 <= lon0 <= 180 and -180 <= lon1 <= 180): return validation_error("Longitude range invalid.")
        if not (0 < dlon <= 360): return validation_error("Longitude sampling invalid.")
        if not (1 <= rep <= 10000): return validation_error("Repetitions must be 1–10000.")
        if any(p < 0 for p in [ulon, ulat, uz, slon, slat, sz] + list(turb.values())): return validation_error("Negative values not allowed.")
        if atm_dt_out <= 0: return validation_error("Output frequency must be > 0.")
        if not (-180 <= central_lon <= 180): return validation_error("Central longitude is out of range.")
        if not (-90 <= central_lat <= 90): return validation_error("Central latitude is out of range.")

        # Check number of air parcels...
        nz, nlat, nlon = map(lambda x: int(round(x[0]/x[1])) + 1, [(z1-z0, dz), (lat1-lat0, dlat), (lon1-lon0, dlon)])
        if nz * nlat * nlon * rep > 10000: return validation_error("Too many trajectories.")

        # Set forward/backward trajectory flag...
        direction = -1 if stop_time < start_time else 1

    except Exception as e:
        return validation_error(str(e))

    # Set run ID...
    run_id = str(uuid.uuid4())

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
    INIT_Z0 = {z0}
    INIT_Z1 = {z1}
    INIT_DZ = {dz}
    INIT_LON0 = {lon0}
    INIT_LON1 = {lon1}
    INIT_DLON = {dlon}
    INIT_LAT0 = {lat0}
    INIT_LAT1 = {lat1}
    INIT_DLAT = {dlat}
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
    DIFFUSION = 1
    """ + '\n'.join(f"{k.upper()} = {v}" for k, v in turb.items()) + f"""
    ATM_BASENAME = {atm_file}
    ATM_DT_OUT = {atm_dt_out}
    """)
    if met_source in ['era5low_6h']:
        ctl_template += "MET_CLAMS = 1\nMET_VERT_COORD = 1\n"
    if met_source in ['merra2_3h', 'merra2_6h']:
        ctl_template += "MET_NC_SCALE = 0\nMET_VERT_COORD = 1\n"
    if met_source in ['ncep2_6h']:
        ctl_template += "MET_RELHUM = 1\n"
    
    with open(ctl_file, 'w') as f: f.write(ctl_template)

    # Run atm_init...
    atm_init_code, atm_init_output = run_command([ATM_INIT_CMD, ctl_file, init_file], timeout=120)

    # Run trac...
    trac_code, trac_output = run_command([TRAC_CMD, dirlist_file, ctl_file, init_file], timeout=600)

    # Catch errors...
    if trac_code == -999:
        return render_template('server_busy.html'), 503

    # Combine log message...
    combined_output = f"=== atm_init Output ===\n{atm_init_output}\n\n=== trac Output ===\n{trac_output}"

    # Run successfull...
    if trac_code == 0:

        # Create zip file...
        zip_path = os.path.join(ZIPS_DIR, f"{run_id}.zip")
        shutil.make_archive(zip_path.replace('.zip', ''), 'zip', work_dir)

        # Gather data...
        lons, lats, heights = [], [], []
        files = sorted(glob.glob(os.path.join(work_dir, 'atm_19*.tab')) +
                       glob.glob(os.path.join(work_dir, 'atm_20*.tab')))
        for file in files:
            data = np.loadtxt(file)
            if data.ndim == 1 and data.shape[0] != 4:
                continue
            data = data[np.newaxis, :] if data.ndim == 1 else data
            lons.extend(data[:, 2])
            lats.extend(data[:, 3])
            heights.extend(data[:, 1])
            
        # Sort by height...
        combined = np.array([lons, lats, heights]).T
        sorted_combined = combined[np.argsort(combined[:, 2])]
        lons_sorted = sorted_combined[:, 0]
        lats_sorted = sorted_combined[:, 1]
        heights_sorted = sorted_combined[:, 2]

        # Plot...
        plot_url = create_plot(
            lons_sorted, lats_sorted, heights_sorted,
            projection=map_projection,
            region=plot_region,
            central_lon=central_lon,
            central_lat=central_lat
        )

        # Clean up...
        delayed_cleanup(work_dir)
        return render_template('result.html', run_id=run_id, stdout=combined_output, plot_url=plot_url)

    # Clean up...
    delayed_cleanup(work_dir)
    return render_template('error.html', stdout=combined_output), 500

@app.route('/terms')
def terms():
    return render_template('terms.html')

if __name__ == '__main__':
    app.run(debug=True)
