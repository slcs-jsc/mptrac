from flask import Flask, render_template, request, send_from_directory, jsonify
import subprocess
import threading
import os
import uuid
import shutil
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

# Maximum number of simulations running at once:
process_semaphore = threading.Semaphore(3)

# Flask app:
app = Flask(__name__)


def run_command(cmd, timeout=300):
    """
    Runs a command with concurrency limit and timeout.
    Combines stdout and stderr into a single output stream.
    """
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
        return -1, (e.stdout or '') + "\n❌ Command timed out after {} seconds.".format(timeout)
    finally:
        process_semaphore.release()

        
def seconds_since_2000(time_str):
    """
    Converts a UTC time string 'YYYY-MM-DD HH:MM:SS'
    to seconds since January 1, 2000, 00:00 UTC.
    """
    dt = datetime.strptime(time_str, "%Y-%m-%d %H:%M:%S")
    dt = dt.replace(tzinfo=timezone.utc)
    ref = datetime(2000, 1, 1, 0, 0, 0, tzinfo=timezone.utc)
    delta = dt - ref
    return delta.total_seconds()


def create_plot(lons, lats, heights):
    """
    Create a scatter plot with longitude and latitude, color-coded by height.
    Returns the plot as a base64 string to be embedded in an HTML page.
    """

    fig = plt.figure(figsize=(15, 12))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()

    # Scatter plot
    sc = ax.scatter(
        lons, lats, c=heights, cmap='viridis',
        s=10, edgecolors='none', transform=ccrs.PlateCarree()
    )

    # Add coastlines
    ax.coastlines(color='white')

    # Add a grid
    gl = ax.gridlines(draw_labels=True, color='gray', linestyle='--', alpha=0.5)
    
    # Customize gridline labels
    gl.xlabel_style = {'color': 'white'}
    gl.ylabel_style = {'color': 'white'}

    # Set background colors
    ax.set_facecolor('#222')
    fig.patch.set_facecolor('#111')

    # Colorbar
    cbar = plt.colorbar(sc, orientation='horizontal', pad=0.05, aspect=50)
    cbar.set_label('Height (km)', color='white')
    cbar.ax.xaxis.set_tick_params(color='white')
    plt.setp(cbar.ax.get_xticklabels(), color='white')

    # Set axis tick colors
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.tick_params(axis='both', colors='white')

    # Title
    ax.set_title('MPTRAC | Air Parcel Trajectories', fontsize=16, color='white')

    # Save to base64
    img = io.BytesIO()
    plt.savefig(img, format='png', bbox_inches='tight')
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode()
    plt.close()

    return plot_url


@app.route('/')
def index():
    return render_template('index.html')


@app.route("/terms")
def terms():
    return render_template("terms.html")


@app.route('/run', methods=['POST'])
def run():
    
    # Collect form inputs...
    try:
        # Time
        start_time = seconds_since_2000(request.form['start_time'])
        stop_time  = seconds_since_2000(request.form['stop_time'])
        
        # Heights (grid)
        z0 = float(request.form['z0'])
        z1 = float(request.form['z1'])
        dz = float(request.form['dz'])
        
        # Latitude (grid)
        lat0 = float(request.form['lat0'])
        lat1 = float(request.form['lat1'])
        dlat = float(request.form['dlat'])
        
        # Longitude (grid)
        lon0 = float(request.form['lon0'])
        lon1 = float(request.form['lon1'])
        dlon = float(request.form['dlon'])
        
        # Uniform random ranges
        ulon = float(request.form['ulon'])
        ulat = float(request.form['ulat'])
        uz   = float(request.form['uz'])
        
        # Gaussian random standard deviations
        slon = float(request.form['slon'])
        slat = float(request.form['slat'])
        sz   = float(request.form['sz'])
        
        # Repetitions
        rep = int(request.form['rep'])
        
        # Diffusivity parameters
        turb_dx_pbl   = float(request.form['turb_dx_pbl'])
        turb_dx_trop  = float(request.form['turb_dx_trop'])
        turb_dx_strat = float(request.form['turb_dx_strat'])
        turb_dz_pbl   = float(request.form['turb_dz_pbl'])
        turb_dz_trop  = float(request.form['turb_dz_trop'])
        turb_dz_strat = float(request.form['turb_dz_strat'])
        turb_mesox    = float(request.form['turb_mesox'])
        turb_mesoz    = float(request.form['turb_mesoz'])
        
        # Output frequency
        atm_dt_out = float(request.form.get('atm_dt_out', 3600))
        
        # Validation
        if not (-100 <= z0 <= 100 and -100 <= z1 <= 100):
            return render_template('error.html', stdout="❌ Invalid height range. Must be between -100 and 100 km."), 400
        if not (-90 <= lat0 <= 90 and -90 <= lat1 <= 90):
            return render_template('error.html', stdout="❌ Invalid latitude range. Must be between -90 and 90°."), 400
        if not (-180 <= lon0 <= 180 and -180 <= lon1 <= 180):
            return render_template('error.html', stdout="❌ Invalid longitude range. Must be between -180 and 180°."), 400
        if not (1 <= rep <= 10000):
            return render_template('error.html', stdout="❌ Invalid repetition value. Must be between 1 and 10000."), 400
        if atm_dt_out <= 0:
            return render_template('error.html', stdout="❌ Output frequency must be > 0."), 400
        
    except Exception as e:
        return render_template('error.html', stdout=f"❌ Invalid input: {str(e)}"), 400
    
    # Create working directory...
    run_id = str(uuid.uuid4())
    work_dir = os.path.join("runs", run_id)
    os.makedirs(work_dir, exist_ok=True)

    try:
        # Set filenames...
        ctl_file = os.path.join(work_dir, f'trac.ctl')
        dirlist_file = os.path.join(work_dir, f'dirlist.txt')
        init_file = os.path.join(work_dir, f'atm_init.tab')
        atm_file = os.path.join(work_dir, f'atm')

        # Generate dirlist...
        with open(dirlist_file, 'w') as f:
            f.write(f"./\n")
        
        # Generate control file...
        config_content = f"""
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
        METBASE = /home/lars/wrk/mptrac/tests/data/ei
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
        """
        with open(ctl_file, 'w') as f:
            f.write(config_content)

        # Run atm_init...
        atm_init_cmd = ['/home/lars/wrk/mptrac/src/atm_init', ctl_file, init_file]
        atm_init_code, atm_init_output = run_command(atm_init_cmd, timeout=120)

        # Run trac...
        trac_cmd = ['/home/lars/wrk/mptrac/src/trac', dirlist_file, ctl_file, init_file]
        trac_code, trac_output = run_command(trac_cmd, timeout=600)
        
        # Server busy...
        if trac_code == -999:
            return render_template('server_busy.html'), 503

        # Success...
        elif trac_code == 0:

            # Create ZIP file...
            zip_filename = f"{run_id}.zip"
            zip_path = os.path.join("runs", zip_filename)
            shutil.make_archive(zip_path.replace(".zip", ""), 'zip', work_dir)
            
            # Find all `.tab` files (ATM files from your simulation)
            files = sorted(glob.glob(os.path.join(work_dir, 'atm_20*.tab')))
            
            # Initialize lists
            lons = []
            lats = []
            heights = []
            
            # Read the data from all files
            for file in files:
                data = np.loadtxt(file)
                lons.extend(data[:, 2])
                lats.extend(data[:, 3])
                heights.extend(data[:, 1])
            
            lons = np.array(lons)
            lats = np.array(lats)
            heights = np.array(heights)
            
            # Create plot...
            plot_url = create_plot(lons, lats, heights)
             
            # Show results...
            return render_template('result.html', run_id=run_id, stdout=trac_output, plot_url=plot_url)

        # Failure...
        else:
            return render_template('error.html', stdout=trac_output), 500
        
    finally:
        
        # Clean working directory...
        if os.path.exists(work_dir):
            shutil.rmtree(work_dir)

            
@app.route('/download/<run_id>')
def download(run_id):
    """Serves the generated simulation ZIP file for download."""
    zip_filename = f"{run_id}.zip"
    return send_from_directory("runs", zip_filename, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)
