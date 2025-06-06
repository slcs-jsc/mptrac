from flask import Flask, render_template, request, send_from_directory
import subprocess, threading, os, uuid, shutil, time, io, base64, glob, textwrap
import numpy as np
from datetime import datetime, timezone
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# === Config ===
ATM_INIT_CMD = os.getenv('ATM_INIT_CMD', '/home/lars/wrk/mptrac/src/atm_init')
TRAC_CMD = os.getenv('TRAC_CMD', '/home/lars/wrk/mptrac/src/trac')
METBASE = os.getenv('METBASE', '/home/lars/wrk/mptrac/tests/data/ei')
RUNS_DIR, ZIPS_DIR = 'runs/working', 'runs/zips'
os.makedirs(RUNS_DIR, exist_ok=True)
os.makedirs(ZIPS_DIR, exist_ok=True)
process_semaphore = threading.Semaphore(3)
app = Flask(__name__)

def delayed_cleanup(path, delay=600):
    threading.Thread(target=lambda: (time.sleep(delay), shutil.rmtree(path, ignore_errors=True)), daemon=True).start()

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

def seconds_since_2000(time_str):
    dt = datetime.strptime(time_str, "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)
    return (dt - datetime(2000, 1, 1, tzinfo=timezone.utc)).total_seconds()

def create_plot(lons, lats, heights):
    fig = plt.figure(figsize=(15, 12))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()
    sc = ax.scatter(lons, lats, c=heights, cmap='viridis', s=10, edgecolors='none', transform=ccrs.PlateCarree())
    ax.coastlines(color='white')
    gl = ax.gridlines(draw_labels=True, color='gray', linestyle='--', alpha=0.5)
    gl.xlabel_style = gl.ylabel_style = {'color': 'white'}
    ax.set_facecolor('#222')
    fig.patch.set_facecolor('#111')
    cbar = plt.colorbar(sc, orientation='horizontal', pad=0.05, aspect=50)
    cbar.set_label('Height (km)', color='white')
    cbar.ax.xaxis.set_tick_params(color='white')
    plt.setp(cbar.ax.get_xticklabels(), color='white')
    ax.tick_params(axis='both', colors='white')
    ax.set_title('MPTRAC | Air Parcel Trajectories', fontsize=16, color='white')
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight')
    plt.close()
    return base64.b64encode(buf.getvalue()).decode()

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
        f = request.form
        to_f = lambda x: float(f[x])
        start_time, stop_time = map(seconds_since_2000, (f['start_time'], f['stop_time']))
        z0, z1, dz = to_f('z0'), to_f('z1'), to_f('dz')
        lat0, lat1, dlat = to_f('lat0'), to_f('lat1'), to_f('dlat')
        lon0, lon1, dlon = to_f('lon0'), to_f('lon1'), to_f('dlon')
        rep = int(f['rep'])
        ulon, ulat, uz = to_f('ulon'), to_f('ulat'), to_f('uz')
        slon, slat, sz = to_f('slon'), to_f('slat'), to_f('sz')
        turb = {k: to_f(k) for k in ['turb_dx_pbl','turb_dx_trop','turb_dx_strat','turb_dz_pbl','turb_dz_trop','turb_dz_strat','turb_mesox','turb_mesoz']}
        atm_dt_out = to_f('atm_dt_out') if 'atm_dt_out' in f else 3600

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

        nz, nlat, nlon = map(lambda x: int(round(x[0]/x[1])) + 1, [(z1-z0, dz), (lat1-lat0, dlat), (lon1-lon0, dlon)])
        if nz * nlat * nlon * rep > 10000: return validation_error("Too many trajectories.")

        direction = -1 if stop_time < start_time else 1

    except Exception as e:
        return validation_error(str(e))

    run_id = str(uuid.uuid4())
    work_dir = os.path.join(RUNS_DIR, run_id)
    os.makedirs(work_dir, exist_ok=True)

    ctl_file, dirlist_file, init_file = map(lambda x: os.path.join(work_dir, x), ['trac.ctl', 'dirlist.txt', 'atm_init.tab'])
    atm_file = os.path.join(work_dir, 'atm')
    with open(dirlist_file, 'w') as f: f.write("./\n")

    ctl_template = textwrap.dedent(f"""\
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
    DT_MET = 86400.0
    DIFFUSION = 1
    """ + '\n'.join(f"{k.upper()} = {v}" for k, v in turb.items()) + f"""
    ATM_BASENAME = {atm_file}
    ATM_DT_OUT = {atm_dt_out}
    """)

    with open(ctl_file, 'w') as f: f.write(ctl_template)

    atm_init_code, atm_init_output = run_command([ATM_INIT_CMD, ctl_file, init_file], timeout=120)
    trac_code, trac_output = run_command([TRAC_CMD, dirlist_file, ctl_file, init_file], timeout=600)

    if trac_code == -999:
        return render_template('server_busy.html'), 503

    combined_output = f"=== atm_init Output ===\n{atm_init_output}\n\n=== trac Output ===\n{trac_output}"
    if trac_code == 0:
        zip_path = os.path.join(ZIPS_DIR, f"{run_id}.zip")
        shutil.make_archive(zip_path.replace('.zip', ''), 'zip', work_dir)

        lons, lats, heights = [], [], []
        for file in sorted(glob.glob(os.path.join(work_dir, 'atm_20*.tab'))):
            data = np.loadtxt(file)
            if data.ndim == 1 and data.shape[0] != 4: continue
            data = data[np.newaxis, :] if data.ndim == 1 else data
            lons.extend(data[:, 2])
            lats.extend(data[:, 3])
            heights.extend(data[:, 1])

        plot_url = create_plot(np.array(lons), np.array(lats), np.array(heights))
        delayed_cleanup(work_dir)
        return render_template('result.html', run_id=run_id, stdout=combined_output, plot_url=plot_url)
    
    delayed_cleanup(work_dir)
    return render_template('error.html', stdout=combined_output), 500

@app.route('/terms')
def terms():
    return render_template('terms.html')

if __name__ == '__main__':
    app.run(debug=True)
