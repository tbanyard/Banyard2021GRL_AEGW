#!/usr/bin/env python3
"""
---------------------
profilesplot_v3.py
---------------------
Plotting code for profiles plot in figure 2 GRL.
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v2.0---Updated for GRL submission: Fixed AIRS profiles and axes------
---v3.0---Final Published Version---------------------------------------
----------[CURRENT]-This_is_the_current_version_of_this_file------------
------------------------------------------------------------------------
========================================================================
Reads data files for CORAL, AIRS, AEOLUS and ERA5 profiles figure.
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
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib
from scipy.io import savemat, loadmat

# Change current working directory to parent directory
os.chdir('..')

# Import from functions file
import sys
from AEGW_functions import *
from cmapstore import *

# matplotlib
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['legend.frameon'] = 'True'
matplotlib.rcParams["legend.facecolor"] = 'inherit'
matplotlib.rcParams["legend.edgecolor"] = 'inherit'

"""=============================ERA5==============================="""
# Load any text files into arrays
ERA5_p_levels = np.loadtxt("ERA5_pressure_levels.txt")*100
ERA5_altitudes = np.loadtxt("ERA5_altitudes.txt")
ERA5_levels = np.loadtxt("ERA5_pressure_levels_labels.txt")

# ERA5 directory
# Data is from time points corresponding to the selected profiles
ERA5_dir = './ERA5_special/'
ncfile = ERA5_dir + 'era5_2019d207_all.nc'
print('ERA5 netCDF file:')
print(ncfile, '\n')
data = ERA5_dataload(ncfile)

# Pass ERA5 data into their respective variables
ERA5_data_lon, ERA5_data_lat, ERA5_data_lev, \
ERA5_data_temp, ERA5_data_u, ERA5_data_v, \
ERA5_data_time = data[0], data[1], data[2], data[3], \
data[4], data[5], data[6]

# Selecting profiles
for i in range(3):
	if i == 0: # 10am
		i_t = 2
	elif i == 1: # 6am
		i_t = 1
	elif i == 2: # 7pm
		i_t = 3
	
	x_temp = ERA5_data_temp[i_t,:,36,74]
	x_u = ERA5_data_u[i_t,:,36,74]
	x_v = ERA5_data_v[i_t,:,36,74]
	phi = 104.65894740762619 # Satellite bearing for this profile
	x_hlosproj = -x_u * np.sin(phi*(np.pi/180)) - x_v * np.cos(phi*(np.pi/180))
	y_era5 = np.arange(80000, 0, -500)
	
	# Temperature
	z_era5 = griddata(ERA5_altitudes, x_temp, y_era5)
	fixnanswithmean(z_era5)
	z_era5_T = np.copy(z_era5)
	
	# Wind
	z_era5 = griddata(ERA5_altitudes, x_hlosproj, y_era5)
	fixnanswithmean(z_era5)
	z_era5_V = np.copy(z_era5)

	# Savitzky-Golay band pass filter of 7-25 profiles
	zTa = savgol_filter(z_era5_T, 25, 2) # + S-G filter Temp
	zTb = savgol_filter(z_era5_T, 7, 2) # - S-G filter Temp
	zVa = savgol_filter(z_era5_V, 25, 2) # + S-G filter Wind
	zVb = savgol_filter(z_era5_V, 7, 2) # - S-G filter Wind
	
	# Getting Perturbations
	if i == 0:
		z_era5_10am = zTb - zTa
		z_era5_10am_hlos = zVb - zVa
	elif i == 1:
		z_era5_6am = zTb - zTa
		z_era5_6am_hlos = zVb - zVa
	elif i == 2:
		z_era5_7pm = zTb - zTa
		z_era5_7pm_hlos = zVb - zVa

"""=============================AEOLUS==============================="""
strdirectory = './AE_Data/'
filename = 'AE_2B_2019-07-26_092323.nc'
infile = strdirectory + str(filename) # Specifies data file
data = nc.Dataset(infile)

# Longitude
data_lon = data.variables['lon'][:]
# Latitude
data_lat = data.variables['lat'][:]
# Altitude
data_alt = data.variables['alt'][:]
# Horizontal Line of Sight Wind speed
data_HLOS = data.variables['Rayleigh_HLOS_wind_speed'][:]
# Zonal Wind Projection
data_u_wind = data.variables['Zonal_wind_projection'][:]
# Meridional Wind Projection
data_v_wind = data.variables['Meridional_wind_projection'][:]
# Line-of-sight azimuth
data_azimuth = data.variables['LOS_azimuth'][:]
# Rayleigh_Grouping
RG = data.variables['RG'][:]
# Time
rayleigh_times = data.variables['time'][:]
# Both QC Flags	
data_QC = data.variables['QC_Flag_Both'][:]	
# Close data file
data.close()

"""=============================================================="""
"""=====Test to see if orbit is sufficiently within Andes box===="""
"""=============================================================="""
minlat, maxlat, minlon, maxlon = -80, -40, 280, 320
mnopib = 150 # minimum_number_of_profiles_in_box
# Print full arrays without truncation
np.set_printoptions(threshold=sys.maxsize)
# ~ print(np.where(data_lat<-80, 0, (np.where(data_lat>-40, 0, 1))))
# Find where the satellite is within the Andes box
box = np.where(data_lat<minlat, 0, (np.where(data_lat>maxlat, 0,
	(np.where(data_lon>maxlon, 0, (np.where(data_lon<minlon, 0, 1)))))))
diffs = np.diff(box) # Array of diffs for box
# Grouping the differences between the elements in diffs
# (2nd derivative)
grouped_diffs = [(k, sum(1 for i in g)) for k,g in groupby(diffs)]
# Returns:[(0, 6206),(1, 1),...,(0, 1748),(-1, 1),...,(0, 8617)]
# ~ print(box)
# ~ print(grouped_diffs)
# Finding the start and end elements of the desired section
itrn = 0
start_elmnt = 0
end_elmnt = 0
complete_boxes = 0
# Iterate through grouped_diffs, adding up continually
for u in grouped_diffs:
	if u[0] != 0: # Bypass 1's and -1's
		itrn += u[1]
	elif u[0] == 0:
		# Is this section the Andes box?
		if data_lat[itrn] < maxlat and data_lat[itrn] > minlat and \
		box[itrn] == 1:
			# Are there enough profiles in the box?
			if u[1] > mnopib:
				if start_elmnt == 0:
					# First profile in box
					start_elmnt = itrn
					# Last profile in box
					end_elmnt = itrn + u[1]
				itrn += u[1]
			else:
				itrn += u[1]
		else:
			itrn += u[1]

	if end_elmnt == 0:
		continue
	
	# Amend arrays
	data_lon_new = data_lon[start_elmnt:end_elmnt+1]
	data_lat_new = data_lat[start_elmnt:end_elmnt+1]
	data_alt_new = data_alt[start_elmnt:end_elmnt+1]
	data_HLOS_new = data_HLOS[start_elmnt:end_elmnt+1]
	rayleigh_times_new = rayleigh_times[start_elmnt:end_elmnt+1]
	data_u_wind_new = data_u_wind[start_elmnt:end_elmnt+1]
	data_v_wind_new = data_v_wind[start_elmnt:end_elmnt+1]
	data_azimuth_new = data_azimuth[start_elmnt:end_elmnt+1]
	data_QC_new = data_QC[start_elmnt:end_elmnt+1]

# Amend arrays again to select profile
data_lon_new = data_lon[5833:6852]
data_lat_new = data_lat[5833:6852]
data_alt_new = data_alt[5833:6852]
data_HLOS_new = data_HLOS[5833:6852]
rayleigh_times_new = rayleigh_times[5833:6852]
data_u_wind_new = data_u_wind[5833:6852]
data_v_wind_new = data_v_wind[5833:6852]
data_azimuth_new = data_azimuth[5833:6852]
data_QC_new = data_QC[5833:6852]

# Initialise meshgrids for x, y and z
# Choose vertical resolution
vert_res = 500 # Given in metres
maxheight = 24000
levnum = (maxheight / vert_res) + 1
alts = np.linspace(0,maxheight, levnum)

# Running new interpolation routine
print('Running Aeolus interpolation routine...\n')
z_new = np.zeros((len(alts),len(RG))) # NumPy Arrays
locs = np.zeros((3, len(RG))) # For topography
points = np.empty(0) # For AE/ERA5
values = np.empty(0) # For AE/ERA5
lons = np.empty(0) # For topography
lats = np.empty(0) # For topography
times = np.empty(0) # For topography
azimuths = np.empty(0) # For ERA5 projection
lastgroupstarttime = 0
RG_start = 0
for RG_elmnt in range(len(RG)):
	for t in range(len(rayleigh_times_new)):
		# Find all elements inside this sandwich and add to z
		# and z_itrn:
		if rayleigh_times_new[t] < RG[RG_elmnt] and \
		rayleigh_times_new[t] >= lastgroupstarttime:
			if RG_start == 0:
				RG_start = RG_elmnt
			# Cap wind speeds to 250 m/s
			if np.abs(data_HLOS_new[t]) < 25000:
				if data_QC_new[t] == 1:
					points = np.append(points, data_alt_new[t])
					values = np.append(values, data_HLOS_new[t])
					azimuths = np.append(azimuths, data_azimuth_new[t])
			# Lons, Lats and Times for later topography
			lons = np.append(lons, data_lon_new[t])
			lats = np.append(lats, data_lat_new[t])
			times = np.append(times, rayleigh_times_new[t])
			# Rest RG_elmnt
			RG_end = RG_elmnt
	lastgroupstarttime = RG[RG_elmnt]
	z_new[:, RG_elmnt] = griddatainterpolation(points, values, alts)
	# locs contains the avg loc info for each RG profile [lons,lats,times]
	locs[0, RG_elmnt], locs[1, RG_elmnt], locs[2, RG_elmnt] = \
		np.mean(lons), np.mean(lats), np.mean(times)
	# Code to select out one profile and save data from it. E.g. Prof 189
	if RG_elmnt == 189:
		prof_values = values
		prof_locs = locs[:,189]
		prof_points = points
		prof_azimuths = azimuths
	points = np.empty(0)
	values = np.empty(0)
	lons = np.empty(0)
	lats = np.empty(0)
	times = np.empty(0)
	azimuths = np.empty(0)

# Extracting interpolated profile 189
prof_values_interp = z_new[:, 189]

"""1D Profile plot"""
# Plotting Aeolus profile
fig = plt.figure()
gridspec.GridSpec(5,5) # Initialise gridspec

# Getting Perturbations
z_aeolus = np.copy(prof_values_interp)
fixnanswithmean(z_aeolus)
z25 = savgol_filter(z_aeolus, 25, 2) # + S-G filter
z7 = savgol_filter(z_aeolus, 7, 2) # - S-G filter
z_aeolus = z7 - z25

"""=========================CORAL lidar=============================="""
strdirectory = './CORAL_Data/'
filename = '20190725-2110_T15Z900.nc'
infile = strdirectory + str(filename) # Specifies data file
data = nc.Dataset(infile)

# Longitude
station_lon = data.variables['station_longitude'][:]
# Latitude
station_lat = data.variables['station_latitude'][:]
# Altitude
data_alt = data.variables['altitude'][:]
# Temperature
data_temp = data.variables['temperature'][:]

# Converted time_offset
data_time_offset = nc.num2date(data.variables['time_offset'][:],\
calendar = 'standard', units = data.variables['time_offset'].units)

# Time units
time_units = 'milliseconds since ' + data_time_offset[0].strftime("%Y-%m-%d %H:%M:%S")

# Converted time
data_time = nc.num2date(data.variables['time'][:],\
calendar = 'standard', units = time_units)

data.close()

# -------------Selecting profile for 10am and finding perturbations------------ #
z_coral = data_temp[145]
z_coral[z_coral == 0] = 'nan'
z_temp = np.copy(z_coral)
fixnanswithmean(z_coral)
z2_coral = savgol_filter(z_coral, 223, 2, axis = 0) # + Vertical S-G filter (15km)
coral_profile_10am = z_temp - z2_coral
coral_alt = data_alt
coral_profile_10am[601:] = np.nan

# -------------Selecting profile for 6am and finding perturbations------------ #
z_coral = data_temp[91]
z_coral[z_coral == 0] = 'nan'
z_temp = np.copy(z_coral)
fixnanswithmean(z_coral)
z2_coral = savgol_filter(z_coral, 223, 2, axis = 0) # + Vertical S-G filter (15km)
coral_profile_6am = z_temp - z2_coral
coral_alt = data_alt

"""==============LOAD AIRS DATA FROM MATLAB FILES============"""
strdirectory = './AIRS_Data/'
filename = 'airs_20190726_055.mat'
infile = strdirectory + str(filename) # Specifies data file
matz1 = loadmat(infile)['z']
matTp1 = loadmat(infile)['Tp']
matTp1[0:6] = np.nan
filename = 'airs_20190726_190.mat'
infile = strdirectory + str(filename) # Specifies data file
matz2 = loadmat(infile)['z']
matTp2 = loadmat(infile)['Tp']
matTp2[0:6] = np.nan

"""=========================================================="""
"""=======================Plotting==========================="""
"""=========================================================="""
os.chdir('Plots')
# Main plot
lw_all = 1.1
# Gridding background
fig = plt.figure(figsize=(10,5))
ax7 = fig.add_subplot(111, label="4", frame_on=True, zorder=0)
ax7.set_yticks(np.arange(5, 60, 5))
ax7.set_ylim(0, 60)
ax7.yaxis.grid(which='both', color='gray', linestyle='--', linewidth=0.4)
for axis in ['top', 'bottom','left','right']: # Set axes thickness
	ax7.spines[axis].set_visible(False)
ax7.tick_params(axis=u'both', which=u'both',length=0)
ax7.set_xticklabels([])
ax7.set_yticklabels([])
# Axis range
axrange = np.linspace(-15,15,7)

# Middle Panel
ax1 = fig.add_subplot(132, label="2", frame_on=True)
# Aeolus Perturbations
plt.plot(z_aeolus[2:]/100, alts[2:]/1000, color = 'black', linestyle = '-', label='Aeolus wind', linewidth = 1.3, zorder=1)
ax2 = ax1.twiny()
# CORAL Perturbations
plt.plot(coral_profile_10am, coral_alt/1000, color = 'orange', linestyle = '--', linewidth = lw_all, label='CORAL temperature', clip_on=False, zorder = 1)
# ERA5 HLOS Perturbations
plt.plot(z_era5_10am, y_era5/1000, color = '#5f98c6', linestyle = '--', linewidth = lw_all, zorder = 1)
plt.plot(z_era5_10am_hlos, y_era5/1000, color = '#5f98c6', linestyle = 'solid', linewidth = lw_all, zorder = 1)
ax1.xaxis.tick_top()
ax1.xaxis.set_label_position('top') 
ax2.xaxis.tick_bottom()
ax2.xaxis.set_label_position('bottom') 

# Cosmetics ax1
ax1.set_yticks([])
ax1.set_yticklabels([])
ax1.set_xticks(axrange)
ax1.set_xlim(-15,15)
ax1.set_ylim(0, 60)
ax1.yaxis.set_major_locator(plt.MaxNLocator(13))
ax1.yaxis.set_minor_locator(plt.MaxNLocator(13))
ax1.set_xlabel('Wind Perturbation / ms$^{-1}$', fontsize = 12)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.tick_params(axis=u'y', which=u'both',length=0)

# Cosmetics ax2
ax2.tick_params(axis=u'y', which=u'both',length=0)
ax2.set_xticks(axrange)
ax2.set_xlim(-15,15)
ax2.set_xlabel('Temperature Perturbation / K', fontsize = 12)
ax2.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)

# Before (AIRS, CORAL, ERA5): i.e. Left Panel
ax3 = fig.add_subplot(131, label="1", frame_on=True)
# AIRS Perturbations
plt.plot(matTp1/2, matz1, color = '#cc193a', linestyle = '--', linewidth = lw_all, label=r"AIRS temperature $ \times$ 0.5", zorder = 1)
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
h3, l3 = ax3.get_legend_handles_labels()
# ERA5 HLOS Temperature Perturbations
plt.plot(z_era5_6am, y_era5/1000, color = '#5f98c6', linestyle = '--', linewidth = lw_all, zorder = 1)
# CORAL Perturbations
plt.plot(coral_profile_6am, coral_alt/1000, color = 'orange', linestyle = '--', linewidth = lw_all, zorder = 1)
ax4 = ax3.twiny()
# ERA5 HLOS Wind Perturbations
plt.plot(z_era5_6am_hlos, y_era5/1000, color = '#5f98c6', linestyle = 'solid', linewidth = lw_all, zorder = 1)
ax4.tick_params(axis=u'y', which=u'both',length=0)
ax4.set_xticks(axrange)
ax4.set_xlim(-15,15)
ax4.spines['right'].set_visible(False)
ax4.spines['left'].set_visible(False)

# Cosmetics ax3
ax3.set_xticks(axrange)
ax3.set_xlim(-15,15)
ax3.set_ylabel('Altitude / km', fontsize = 12)
ax3.set_yticks(np.arange(0, 60, 5))
ax3.set_ylim(0, 60)
ax3.yaxis.set_major_locator(plt.MaxNLocator(13))
ax3.yaxis.set_minor_locator(plt.MaxNLocator(13))
ax3.spines['right'].set_visible(False)

# After (AIRS, ERA5): i.e. Right Panel
ax5 = fig.add_subplot(133, label="3", frame_on=True)
# ERA5 HLOS Temperature Perturbations
plt.plot(z_era5_7pm, y_era5/1000, color = '#5f98c6', linestyle = '--', linewidth = lw_all, label = 'ERA5 temperature', zorder = 1)
# AIRS Perturbations
plt.plot(matTp2/2, matz2, color = '#cc193a', linestyle = '--', linewidth = lw_all, zorder = 1)
ax6 = ax5.twiny()
# ERA5 HLOS Wind Perturbations
plt.plot(z_era5_7pm_hlos, y_era5/1000, color = '#5f98c6', linestyle = 'solid', linewidth = lw_all, label = 'ERA5 wind', zorder = 1)
ax6.tick_params(axis=u'y', which=u'both',length=0)
ax6.set_xticks(axrange)
ax6.set_xlim(-15,15)
ax6.spines['right'].set_visible(False)
ax6.spines['left'].set_visible(False)

# Cosmetics ax5
ax5.set_xticks(axrange)
ax5.set_xlim(-15,15)
ax5.set_yticks(np.arange(0, 60, 5))
ax5.set_ylim(0, 60)
ax5.set_ylabel('Altitude / km', fontsize = 12)
ax5.yaxis.set_label_position("right")
ax5.yaxis.set_major_locator(plt.MaxNLocator(13))
ax5.yaxis.set_minor_locator(plt.MaxNLocator(13))
ax5.yaxis.tick_right()
ax5.spines['left'].set_visible(False)

# All grids
for ax in [ax1, ax3, ax5]:
	ax.yaxis.grid(which='both', color='gray', linestyle='--', linewidth=0.4)
	for num in [-10, 0, 10]:
		ax.axvline(num, color='gray', lw=0.4, linestyle='--', zorder = 0)
for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
	ax.tick_params(axis='both', labelsize = 12)

# Add legend to figure
fig.subplots_adjust(bottom=0.2, right=0.88, left=0.12)
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
h3, l3 = ax3.get_legend_handles_labels()
h5, l5 = ax5.get_legend_handles_labels()
h6, l6 = ax6.get_legend_handles_labels()
leg = plt.legend(h1+h6+h5+h2+h3, l1+l6+l5+l2+l3, loc=3, bbox_to_anchor=(-2.68, -0.28, 8, 1.2), fontsize = 10.7, ncol=8, columnspacing=0.8, handletextpad=0.4).set_zorder(10)

# Get time and location of profile
profile_time = coda.time_to_utcstring(RG[189]) # End of RG group
profile_time = coda.time_to_utcstring(prof_locs[2])
profile_time = datetime.strptime(profile_time,
			'%Y-%m-%d %H:%M:%S.%f')
plottitle = str(profile_time) + "\n" + str(prof_locs[0]) + " longitude\n" + str(prof_locs[1]) + " latitude"
plt.savefig('Profiles_Plot.png', dpi=600)
print("\nProfiles_Plot.png saved in: ", os.getcwd())
