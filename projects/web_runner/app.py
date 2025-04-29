from flask import Flask, render_template, request, send_from_directory, jsonify
import subprocess
import threading
import os
import uuid
import shutil
import time
from datetime import datetime, timezone
import io
import base64
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import textwrap

# === Configuration ===
ATM_INIT_CMD = os.environ.get('ATM_INIT_CMD', '/home/lars/wrk/mptrac/src/atm_init')
TRAC_CMD = os.environ.get('TRAC_CMD', '/home/lars/wrk/mptrac/src/trac')
METBASE = os.environ.get('METBASE', '/home/lars/wrk/mptrac/tests/data/ei')

RUNS_DIR = 'runs/working'
ZIPS_DIR = 'runs/zips'
os.makedirs(RUNS_DIR, exist_ok=True)
os.makedirs(ZIPS_DIR, exist_ok=True)

# Maximum number of simulations running at once:
process_semaphore = threading.Semaphore(3)

# Flask app:
app = Flask(__name__)


def delayed_cleanup(path, delay=600):
    def cleanup():
        time.sleep(delay)
        if os.path.exists(path):
            shutil.rmtree(path)
    threading.Thread(target=cleanup, daemon=True).start()


def run_command(cmd, timeout=300):
    """Runs a command with concurrency limit and timeout."""
    acquired = process_semaphore.acquire(timeout=5)
    if not acquired:
        return -999, "❌ Server is too busy."

    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=timeout
        )
        return result.returncode, result.stdout
    except subprocess.TimeoutExpired as e:
        return -1, (e.stdout or '') + f"\n❌ Command timed out after {timeout} seconds."
    finally:
        process_semaphore.release()


def seconds_since_2000(time_str):
    """Converts a UTC time string 'YYYY-MM-DD HH:MM:SS' to seconds since January 1, 2000."""
    dt = datetime.strptime(time_str, "%Y-%m-%d %H:%M:%S")
    dt = dt.replace(tzinfo=timezone.utc)
    ref = datetime(2000, 1, 1, 0, 0, 0, tzinfo=timezone.utc)
    return (dt - ref).total_seconds()


def create_plot(lons, lats, heights):
    """Create a scatter plot."""
    fig = plt.figure(figsize=(15, 12))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()

    sc = ax.scatter(lons, lats, c=heights, cmap='viridis', s=10, edgecolors='none', transform=ccrs.PlateCarree())

    ax.coastlines(color='white')
    gl = ax.gridlines(draw_labels=True, color='gray', linestyle='--', alpha=0.5)
    gl.xlabel_style = {'color': 'white'}
    gl.ylabel_style = {'color': 'white'}

    ax.set_facecolor('#222')
    fig.patch.set_facecolor('#111')

    cbar = plt.colorbar(sc, orientation='horizontal', pad=0.05, aspect=50)
    cbar.set_label('Height (km)', color='white')
    cbar.ax.xaxis.set_tick_params(color='white')
    plt.setp(cbar.ax.get_xticklabels(), color='white')

    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.tick_params(axis='both', colors='white')
    ax.set_title('MPTRAC | Air Parcel Trajectories', fontsize=16, color='white')

    img = io.BytesIO()
    plt.savefig(img, format='png', bbox_inches='tight')
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode()
    plt.close()

    return plot_url


def validation_error(message):
    return render_template('error.html', stdout=f"❌ {message}"), 400


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/terms')
def terms():
    return render_template('terms.html')


@app.route('/run', methods=['POST'])
def run():
    try:
        start_time = seconds_since_2000(request.form['start_time'])
        stop_time  = seconds_since_2000(request.form['stop_time'])

        z0 = float(request.form['z0'])
        z1 = float(request.form['z1'])
        dz = float(request.form['dz'])

        lat0 = float(request.form['lat0'])
        lat1 = float(request.form['lat1'])
        dlat = float(request.form['dlat'])

        lon0 = float(request.form['lon0'])
        lon1 = float(request.form['lon1'])
        dlon = float(request.form['dlon'])

        ulon = float(request.form['ulon'])
        ulat = float(request.form['ulat'])
        uz   = float(request.form['uz'])

        slon = float(request.form['slon'])
        slat = float(request.form['slat'])
        sz   = float(request.form['sz'])

        rep = int(request.form['rep'])

        turb_dx_pbl   = float(request.form['turb_dx_pbl'])
        turb_dx_trop  = float(request.form['turb_dx_trop'])
        turb_dx_strat = float(request.form['turb_dx_strat'])
        turb_dz_pbl   = float(request.form['turb_dz_pbl'])
        turb_dz_trop  = float(request.form['turb_dz_trop'])
        turb_dz_strat = float(request.form['turb_dz_strat'])
        turb_mesox    = float(request.form['turb_mesox'])
        turb_mesoz    = float(request.form['turb_mesoz'])

        atm_dt_out = float(request.form.get('atm_dt_out', 3600))

        if not (-100 <= z0 <= 100 and -100 <= z1 <= 100):
            return validation_error("Invalid height range. Must be between -100 and 100 km.")
        if not (-90 <= lat0 <= 90 and -90 <= lat1 <= 90):
            return validation_error("Invalid latitude range. Must be between -90 and 90°.")
        if not (-180 <= lon0 <= 180 and -180 <= lon1 <= 180):
            return validation_error("Invalid longitude range. Must be between -180 and 180°.")
        if not (1 <= rep <= 10000):
            return validation_error("Invalid repetition value. Must be between 1 and 10000.")
        if atm_dt_out <= 0:
            return validation_error("Output frequency must be > 0.")

    except Exception as e:
        return validation_error(str(e))

    run_id = str(uuid.uuid4())
    work_dir = os.path.join(RUNS_DIR, run_id)
    os.makedirs(work_dir, exist_ok=True)

    ctl_file = os.path.join(work_dir, 'trac.ctl')
    dirlist_file = os.path.join(work_dir, 'dirlist.txt')
    init_file = os.path.join(work_dir, 'atm_init.tab')
    atm_file = os.path.join(work_dir, 'atm')

    with open(dirlist_file, 'w') as f:
        f.write("./\n")

    config_content = textwrap.dedent(f"""\
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
    T_STOP = {stop_time}
    METBASE = {METBASE}
    DT_MET = 86400.0
    DIFFUSION = 1
    TURB_DX_PBL = {turb_dx_pbl}
    TURB_DX_TROP = {turb_dx_trop}
    TURB_DX_STRAT = {turb_dx_strat}
    TURB_DZ_PBL = {turb_dz_pbl}
    TURB_DZ_TROP = {turb_dz_trop}
    TURB_DZ_STRAT = {turb_dz_strat}
    TURB_MESOX = {turb_mesox}
    TURB_MESOZ = {turb_mesoz}
    ATM_BASENAME = {atm_file}
    ATM_DT_OUT = {atm_dt_out}
    """)

    with open(ctl_file, 'w') as f:
        f.write(config_content)

    atm_init_cmd = [ATM_INIT_CMD, ctl_file, init_file]
    atm_init_code, atm_init_output = run_command(atm_init_cmd, timeout=120)

    trac_cmd = [TRAC_CMD, dirlist_file, ctl_file, init_file]
    trac_code, trac_output = run_command(trac_cmd, timeout=600)

    if trac_code == -999:
        return render_template('server_busy.html'), 503

    if trac_code == 0:
        zip_filename = f"{run_id}.zip"
        zip_path = os.path.join(ZIPS_DIR, zip_filename)
        shutil.make_archive(zip_path.replace('.zip', ''), 'zip', work_dir)

        files = sorted(glob.glob(os.path.join(work_dir, 'atm_20*.tab')))
        lons, lats, heights = [], [], []

        for file in files:
            data = np.loadtxt(file)
            lons.extend(data[:, 2])
            lats.extend(data[:, 3])
            heights.extend(data[:, 1])

        plot_url = create_plot(np.array(lons), np.array(lats), np.array(heights))

        delayed_cleanup(work_dir)

        return render_template('result.html', run_id=run_id, stdout=trac_output, plot_url=plot_url)

    else:
        delayed_cleanup(work_dir)
        return render_template('error.html', stdout=trac_output), 500


@app.route('/download/<run_id>')
def download(run_id):
    """Serves the generated simulation ZIP file for download."""
    zip_filename = f"{run_id}.zip"
    return send_from_directory(ZIPS_DIR, zip_filename, as_attachment=True)


if __name__ == '__main__':
    app.run(debug=True)
