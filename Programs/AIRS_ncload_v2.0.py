#!/usr/bin/env python3
"""
---------------------
AIRS_ncload_v2.0.py
---------------------
Load and plot AIRS data from MATLAB format
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Additional Flexibility and uses the Lambert conformal conic---
----------projection.
---v1.2-----------------------------------------------------------------
---v1.3---Fixed for final GRL submission--------------------------------
---v2.0---Final Published Version---------------------------------------
----------[CURRENT]-This_is_the_current_version_of_this_file------------
------------------------------------------------------------------------
========================================================================
Reads .mat files of AIRS data and plots along scan track over the Andes
========================================================================
"""

# Imports
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import matplotlib.gridspec as gridspec
import matplotlib.text as text
import os
import coda
import errno
from datetime import timedelta, datetime
import time
from scipy.interpolate import griddata
from scipy.io import savemat
from scipy.signal import savgol_filter, butter
import scipy.ndimage as ndimage
from itertools import groupby
from mpl_toolkits.basemap import Basemap, pyproj
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import random
from scipy.io import savemat, loadmat
import numpy.ma as ma
import matplotlib

# Change current working directory to parent directory
os.chdir('..')

# Import from functions file
import sys
from AEGW_functions import *
from cmapstore import *

# cmaps
rdbucmap = LinearSegmentedColormap('RDBUcmap', segmentdata=fetchcmap('RdBunew'), N=265)
whitehatchescmap_r = LinearSegmentedColormap('Whitehatchescmap', segmentdata=fetchcmap('whitehatches'), N=265)
whitehatchescmap = whitehatchescmap_r.reversed()

# matplotlib settings
matplotlib.rcParams['axes.unicode_minus'] = False

########################################################################
####### The settings in this section can be altered by the user ########
########################################################################
# Please note, the plots created by this program use MATLAB files      #
# created by Corwin Wright, saved as two separate files for the before #
# and after orbital sections, as a concatenation of AIRS granules      #
# 54 + 55 + 56 and 189 + 190 + 191. Below, the original netCDF files   #
# are also loaded, in case the user prefers to use these instead.      #
########################################################################
# Find directory and read AIRS netCDF data
strdirectory = './AIRS_Data/'
filename = 'airs_2019_207_055.nc'
infile = strdirectory + str(filename) # Specifies data file
print("\n===========================================================")
print('AIRS netCDF file path:', infile, '\n')
data = nc.Dataset(infile)

# Use AIRS data instead to produce plots. Toggle the before and after  #
# overpasses below:
filename = 'airs_054055056.mat' # Before
# ~ filename = 'airs_189190191.mat' # After

# Select desired map projection. Options: 'lcc', 'cyl'
proj = 'lcc'

# Select desired colorbar orientation. Options: 'horiz', 'verti'
horiz_verti = 'verti'
########################################################################
########################################################################

"""=============================================================="""
"""======================Download Variables======================"""
"""=============================================================="""
# Longitude
data_lon = data.variables['l1_lon'][:]
# Latitude
data_lat = data.variables['l1_lat'][:]
# Altitude
data_alt = data.variables['ret_z'][:]
# Temperature
data_temp = data.variables['ret_temp'][:]
data.close()

# Calculating perturbations using a Savitzky-Golay filter
data_temp2 = savgol_filter(data_temp, 55, 2, axis = 1)
data_temp = data_temp - data_temp2
data_temp = ndimage.gaussian_filter(data_temp, sigma=0.75, order=0)
# data_temp may be plotted under the user's own plotting regime

"""=============LOAD AIRS INSTEAD FROM MATLAB FILES=============="""
infile = strdirectory + str(filename) # Specifies data file
matlat1 = loadmat(infile)['Lat'] # Latitudes
matlon1 = loadmat(infile)['Lon'] # Longitudes
matTp1 = loadmat(infile)['Tp'] # Temperature Perturbations
matlat1, matlon1, matTp1 = ma.asarray(matlat1), ma.asarray(matlon1), ma.asarray(matTp1)

"""=========================================================="""
"""=======================Plotting==========================="""
"""=========================================================="""
print("Current working directory:", os.getcwd())
os.chdir('Plots') # Change directory to Plots folder

# Initialise figure
print('\nInitialise plot...\n')
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)

# Plotting map
if proj == 'cyl':
	bmlowlat, bmupperlat, bmleftlon, bmrightlon = -70, -40, -100, -50
	map = Basemap(projection='cyl',llcrnrlat=bmlowlat,urcrnrlat=bmupperlat,\
				llcrnrlon=bmleftlon,urcrnrlon=bmrightlon,resolution='i', ax=ax1)
	# ~ map.fillcontinents(color='#ffdd99', lake_color='#cceeff')
	# ~ map.drawmapboundary(linewidth=0.75, fill_color='#cceeff')
	map.drawcoastlines(linewidth=0.25, color='black', zorder=3)
	map.drawmeridians([-90,-80,-70, -60, -50], linewidth=0.3)
	map.drawparallels([-70, -60, -50], linewidth=0.3)

	# Plotting
	cs = plt.contourf(matlon1, matlat1, matTp1, cmap=rdbucmap, zorder=2, vmin = -20, vmax = 20, levels=np.linspace(-20,20,17), extend='both')
	
	# Fix axes
	ax1.set_xticks([-100,-90,-80, -70, -60, -50])
	ax1.set_yticks([-70, -60, -50, -40])
	ax1.set_xlabel('Longitude / deg')
	ax1.set_ylabel('Latitude / deg')
	ax1.set_aspect('auto') # Stretch map to fill subplot
	for axis in ['top','bottom','left','right']: # Set axes thickness
		ax1.spines[axis].set_linewidth(0.75)

if proj == 'lcc':
	# Create plot border mask
	border_lons = np.arange(-120, 0.01, 0.01)
	border_lats = np.arange(-80, -29.99, 0.01)
	border_lons, border_lats = np.meshgrid(border_lons, border_lats)
	box = np.where(border_lats<-70.05, 0, (np.where(border_lats>-39.95, 0,
		(np.where(border_lons>-49.95, 0, (np.where(border_lons<-90.05, 0, 1)))))))
	
	# Plot basemap
	map = Basemap(width=3450000,height=3560000,
				rsphere=(6378137.00,6356752.3142),\
				resolution='l',area_thresh=250.,projection='lcc',\
				lat_1=-45.,lat_2=-65,lat_0=-55,lon_0=-70.)
	projstring = map.proj4string
	
	# Transforming data latitudes and longitudes to LCC
	pcyl = pyproj.Proj(init="epsg:4326")
	plcc = pyproj.Proj(projstring)
	data_lon, data_lat = pyproj.transform(pcyl, plcc, data_lon, data_lat)
	matlon1, matlat1 = pyproj.transform(pcyl, plcc, matlon1.transpose(), matlat1.transpose())
	
	# Transforming border mask latitudes and longitudes to LCC
	border_lons, border_lats = pyproj.transform(pcyl, plcc, border_lons, border_lats)

	# Plotting
	plt.plot(matlon1[:50], matlat1[:50])
	cs = map.contourf(matlon1, matlat1, matTp1.transpose(), cmap=rdbucmap, zorder=2, vmin = -20, vmax = 20, levels=np.linspace(-20,20,17), extend='both')
	cs2 = map.contourf(border_lons, border_lats, box, cmap = whitehatchescmap_r, zorder=4)
		
	# Map cosmetics
	for axis in ['top','bottom','left','right']: # Set axes thickness
		ax1.spines[axis].set_linewidth(0)
	map.drawcoastlines(linewidth=0.25, color='black', zorder=3)
	map.drawmeridians(np.linspace(-160, 90, 26), linewidth=0.3)
	map.drawparallels(np.linspace(-70, -20, 6), linewidth=0.3)
	# Map edges
	map.drawmeridians([-90,-50], linewidth=1, zorder=3, dashes=(None,None))
	map.drawparallels([-40,-70], linewidth=1, zorder=3, dashes=(None,None))
	
	# Labelling meridians
	meridians_x = np.arange(-90, -40, 20)
	meridians_y = np.full(len(meridians_x), -70)
	meridian_labels_x, meridian_labels_y = pyproj.transform(pcyl, plcc, meridians_x, meridians_y)
	for i in range(len(meridians_x)):
		if horiz_verti == 'horiz':
			plt.annotate(str(abs(meridians_x[i]))+"$\!^\circ\!$"+"W", xy=(meridian_labels_x[i], meridian_labels_y[i]), xytext=(meridian_labels_x[i]-250000, meridian_labels_y[i]-230000), zorder = 5, fontsize=12)
		elif horiz_verti == 'verti':
			plt.annotate(str(abs(meridians_x[i]))+"$\!^\circ\!$"+"W", xy=(meridian_labels_x[i], meridian_labels_y[i]), xytext=(meridian_labels_x[i]-250000, meridian_labels_y[i]-230000), zorder = 5, fontsize=14)
		else:
			print('Error: Invalid colorbar orientation')
			sys.exit(1)
		
	# Labelling parallels
	parallels_y = np.arange(-70, -30, 10)
	parallels_x = np.full(len(parallels_y), -90)
	parallel_labels_x, parallel_labels_y = pyproj.transform(pcyl, plcc, parallels_x, parallels_y)
	for i in range(len(parallels_y)):
		if horiz_verti == 'horiz':
			plt.annotate(str(abs(parallels_y[i]))+"$\!^\circ\!$"+"S", xy=(parallel_labels_x[i], parallel_labels_y[i]), xytext=(parallel_labels_x[i]-525000, parallel_labels_y[i]-40000), zorder = 5, fontsize=12)
		elif horiz_verti == 'verti':
			plt.annotate(str(abs(parallels_y[i]))+"$\!^\circ\!$"+"S", xy=(parallel_labels_x[i], parallel_labels_y[i]), xytext=(parallel_labels_x[i]-525000, parallel_labels_y[i]-40000), zorder = 5, fontsize=14)
		else:
			print('Error: Invalid colorbar orientation')
			sys.exit(1)
	
else:
	print('Error: Invalid map projection string')
	sys.exit(1)

# Title
if filename == 'airs_054055056.mat':
	plt.title('AIRS granules 54 to 56 for 2019-07-26 at 30 km')
elif filename == 'airs_189190191.mat':
	plt.title('AIRS granules 189 to 191 for 2019-07-26 at 30 km')

# Add horizontal colorbar to figure
if horiz_verti == 'horiz':
	fig.subplots_adjust(bottom=0.3, right=0.88, left=0.12)
	cbar_ax = fig.add_axes([0.27, 0.15, 0.46, 0.025])
	fig.colorbar(cs, cmap=rdbucmap, orientation='horizontal',
		label='Temperature Perturbation / K', cax=cbar_ax, ticks=np.linspace(-20,20,9))
	cbar_ax.tick_params(labelsize=14)
	cbar_ax.xaxis.label.set_size(14)

# Add vertical colorbar to figure
elif horiz_verti == 'verti':
	fig.subplots_adjust(bottom=0.1, right=0.75, left=0.125)
	cbar_ax = fig.add_axes([0.8, 0.05, 0.02, 0.9])
	fig.colorbar(cs, cmap=rdbucmap, orientation='vertical',
		label='Temperature Perturbation / K', cax=cbar_ax, ticks=np.linspace(-20,20,9))
	cbar_ax.tick_params(labelsize=14)
	cbar_ax.yaxis.label.set_size(14)

# Saving figure
if filename == 'airs_054055056.mat':
	pngsavename = 'airs_054055056.png'
elif filename == 'airs_189190191.mat':
	pngsavename = 'airs_189190191.png'
	
plt.savefig(pngsavename, dpi=600)
print("\n", pngsavename, "saved in: ", os.getcwd())

plt.close()
