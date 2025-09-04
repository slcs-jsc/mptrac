#! /p/software/juwels/stages/2025/software/Python/3.12.3-GCCcore-13.3.0/bin/python
from glob import glob
from tqdm import tqdm
import matplotlib.pylab as plt
import numpy as np
from glob import glob
import xarray as xa
import matplotlib.pylab as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Hack to get the ltm_tools...
import sys
import os
path_ltm_tools="../../../"
sys.path.append(os.path.abspath(path_ltm_tools))

# Load ltm_tools...
from ltm_tools.io_manager import read_atm_pos, read_met, write_atm_pos

# Set-up...
atm_prefix = "atm"

# Define pathways to data...
setup_ext_1 = "data"
data_paths_files_data0 = f"./{setup_ext_1}/data.0/atm_2022_06_*_*_00.tab"
data_paths_1 = f"./{setup_ext_1}/data*/"
data_fig_path = "./plots/"

# Plot different time steps...	
init = True
init_legend = True
for ind, data_path_1 in enumerate(np.sort(glob(data_paths_files_data0))):
	date = data_path_1[-21-len(atm_prefix):]
	print(f"[INFO] Plotting {data_paths_1 + date} ...")
	print([f for f in glob(data_paths_1 + date)])
	atm1 = xa.concat( [read_atm_pos(f) for f in glob(data_paths_1 + date)], dim="index")
	atm1 = atm1.where( atm1["subdomain"] != -1, drop=True)
	
	atm1["idx"] = atm1["idx"]
	
	atm1["lon"][atm1["lon"]<0] += 360

	idx1 = np.unique(atm1["idx"])

	dublicates = idx1[np.unique(atm1["idx"], return_counts=True)[1]>1]

	for idx in dublicates:
		mask = ( atm1["idx"] == idx )
		
	atm1 = atm1.sortby("idx")
	
	if init:
		subdomains = np.array(atm1["subdomain"])
		init = False

	# Plot domains...
	plt.figure()
	ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
	ax.add_feature(cfeature.LAND, edgecolor='gray')	
	
	# Set the longitude and latitude window
	lon_min, lon_max = -180, 180  # Example values, adjust as needed
	lat_min, lat_max = -90, 90  # Example values, adjust as needed
	ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree(central_longitude=180))
	
	ndomains_meridional = 3
	ndomains_zonal = 3
	vlines = np.linspace(-180,180, ndomains_zonal+1)
	hlines = np.linspace(-90,90, ndomains_meridional+1) 
	for vline in vlines:
	    ax.axvline(vline, c="k", lw=0.5)
	for hline in hlines:
	    ax.axhline(hline, c="k", lw=0.5)
	    	    
	# Plot halos...
	vlines = np.linspace(-180,180, ndomains_zonal+1)+2
	hlines = np.linspace(-90,90, ndomains_meridional+1)+2
	for vline in vlines:
	    ax.axvline(vline, c="k", lw=0.5)
	for hline in hlines[1:-1]:
	    ax.axhline(hline, c="k", lw=0.5)
	    
	vlines = np.linspace(-180,180, ndomains_zonal+1)-2
	hlines = np.linspace(-90,90, ndomains_meridional+1)-2
	for vline in vlines:
	    ax.axvline(vline, c="k", lw=0.5)
	for hline in hlines[1:-1]:
	    ax.axhline(hline, c="k", lw=0.5)

	vcentres = 0.5*(vlines[1:]+vlines[:-1])
	hcentres = 0.5*(hlines[1:]+hlines[:-1])
	dhcentre = 0.5*(hcentres[1]-hcentres[0])
	dvcentre = 0.5*(vcentres[1]-vcentres[0])
	
	ndm = 0
	for vline in vcentres[::]:
		for hline in hcentres[::-1]:
			ax.text(vline, hline, f"D{ndm}")
			ndm+=1
	
	mask = (subdomains >=0)
	atm1["lon"][mask]+=180
	
	sdd=5.0
	sbase=1.0 #0.5
	thr=1e-16
	#if init_legend:
	ax.scatter( atm1["lon"][mask], atm1["lat"][mask], c=subdomains[mask], marker=".", cmap="hsv",vmin=0, vmax=np.max(subdomains), label="DD", s=sdd)
	
	# Adding labels
	ax.set_title("3x3 Domain Decomposition", fontsize=14)
	ax.text(0.5, -0.10, 'Longitude [°]', va='bottom', ha='center', transform=ax.transAxes)
	ax.text(-0.07, 0.5, 'Latitude [°]', va='center', ha='center', rotation='vertical', transform=ax.transAxes)

	ax.set_xlim(-180,180)
	ax.set_ylim(-90,90)
	#plt.legend(loc='upper center', bbox_to_anchor=(0.2, -0.01), ncols=2, fontsize=9)
	
	plt.savefig(f"./{data_fig_path}/{data_path_1[-24:-4]}.png", dpi=300, bbox_inches='tight')	
		

        
