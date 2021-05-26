#!/usr/bin/env python3
"""
---------------------
CORAL_ncload_v2.0.py
---------------------
CORAL data load from netCDF format
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Timeseries using qbocmap and raw data-------------------------
---v1.2---Improving plot------------------------------------------------
---v1.3---Time domain flexibility and improved appearance for fig2 GRL--
---v2.0---Final Published Version---------------------------------------
----------[CURRENT]-This_is_the_current_version_of_this_file------------
------------------------------------------------------------------------
========================================================================
Reads .nc files of data from the CORAL lidar in Argentina, and plots it.
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
import random
import matplotlib

# Change current working directory to parent directory
os.chdir('..')

# Import from functions file
import sys
from AEGW_functions import *
from cmapstore import *

# cmaps
rdbucmap = LinearSegmentedColormap('RDBUcmap', segmentdata=fetchcmap('RdBunew'), N=265)

# matplotlib settings
matplotlib.rcParams['axes.unicode_minus'] = False

########################################################################
####### The settings in this section can be altered by the user ########
########################################################################
strdirectory = './CORAL_Data/'
filename = '20190725-2110_T15Z900.nc'
infile = strdirectory + str(filename) # Specifies data file
print("\n===========================================================")
print('CORAL netCDF file path:', infile, '\n')
data = nc.Dataset(infile)
########################################################################
########################################################################

"""=============================================================="""
"""======================Download Variables======================"""
"""=============================================================="""
# Longitude
station_lon = data.variables['station_longitude'][:]
# Latitude
station_lat = data.variables['station_latitude'][:]
# Altitude
data_alt = data.variables['altitude'][:]
# Temperature
data_temp = data.variables['temperature'][:]

# Time offset units
print("Time offset units:", data.variables['time_offset'].units)

# Converted time_offset
data_time_offset = nc.num2date(data.variables['time_offset'][:],\
calendar = 'standard', units = data.variables['time_offset'].units)

# Time units
time_units = 'milliseconds since ' + data_time_offset[0].strftime("%Y-%m-%d %H:%M:%S")
print("Time units:", time_units)

# Converted time
data_time = nc.num2date(data.variables['time'][:],\
calendar = 'standard', units = time_units)
# ~ print(data_time)
data.close()

"""=========================================================="""
"""=======================Plotting==========================="""
"""=========================================================="""
print("\nCurrent working directory:", os.getcwd())
os.chdir('Plots')

# Initialise figure
print('\nInitialise raw temperatures plot...\n')
fig = plt.figure(figsize=(10,5))
ax1 = fig.add_subplot(1, 1, 1)

# Fixing time dimension and setting up domain
data_time = np.array([datetime.strftime(date,
	'%Y-%m-%d %H:%M:%S.%f') for date in data_time])
data_time = np.array([datetime.strptime(date,
	'%Y-%m-%d %H:%M:%S.%f') for date in data_time])
data_alt = data_alt[100:1001]
x, y = np.meshgrid(data_time, data_alt)

# Recasting x and y arrays
x_lims = [np.ndarray.flatten(x)[0], np.ndarray.flatten(x)[-1]]
x_lims = dates.date2num(x_lims)
# Y limits
y_lims = [np.ndarray.flatten(y)[0]/1000, np.ndarray.flatten(y)[-1]/1000]

# Set z
z = np.flip(np.transpose(data_temp[:,100:1001]), axis = 0)
z[z == 0] = 'nan'

# Plotting
cs = plt.imshow(z, aspect='auto', cmap='magma', extent=[x_lims[0],
				x_lims[1], y_lims[0], y_lims[1]], vmin=150, vmax=300,
				interpolation='none')

# Date axis				
ax1.xaxis_date() # Initialises date axis
date_form = dates.DateFormatter('%H:%M') # Sets date format
ax1.xaxis.set_major_formatter(date_form)
# ~ plt.gca().invert_yaxis() # Invert axis for imshow

# Cosmetics
plt.title('CORAL lidar timeseries for 2019-07-25 - 2019-07-26', fontdict = {'fontsize': 15})
ax1.set_ylabel('Altitude / km')
ax1.set_yticks(np.linspace(10,100,10))
ax1.yaxis.set_ticks_position('both')
ax1.tick_params(axis='y', which='both', labelleft='on', labelright='on')
ax1.set_xlabel('Time')
ax1.xaxis.grid(True, which='minor', color='gray', linestyle='--', linewidth=0.4)
ax1.grid(color='gray', linestyle = '--', linewidth = 0.4, axis='x',
			which='both')
ax1.grid(color='gray', linestyle = '--', linewidth = 0.4, axis='y',
			which='both')
			
# Ensure the number of date ticks is sensible
ax1.xaxis.set_minor_locator(dates.MinuteLocator(byminute=0))

# Label fontsize
ax1.tick_params(labelsize=13)
ax1.yaxis.label.set_size(13)
ax1.xaxis.label.set_size(13)

# Add colorbar to figure
fig.subplots_adjust(bottom=0.3, right=0.88, left=0.12)
cbar_ax = fig.add_axes([0.12, 0.15, 0.76, 0.025])
fig.colorbar(cs, cmap='magma', orientation='horizontal',
	label='Temperature / K', cax=cbar_ax)
cbar_ax.tick_params(labelsize=13)
cbar_ax.xaxis.label.set_size(13)
			
# Saving figure
plt.savefig('CORAL_timeseries.png',dpi=600)
print("CORAL_timeseries.png saved in: ", os.getcwd())

plt.close()

"""Perturbations"""
# Initialise figure
print('\nInitialise temperature perturbations plot...')
fig = plt.figure(figsize=(10,5))
ax1 = fig.add_subplot(1, 1, 1)

# Set z
z_temp = np.copy(z)
fixnanswithmean(z)
z2 = savgol_filter(z, 223, 2, axis = 0) # + Vertical S-G filter (15km)
z = z_temp - z2

# Plotting
cs = plt.imshow(z, aspect='auto', cmap=rdbucmap, extent=[x_lims[0],
				x_lims[1], y_lims[0], y_lims[1]], vmin=-20, vmax=20,
				interpolation='none')

# Date axis				
ax1.xaxis_date() # Initialises date axis
date_form = dates.DateFormatter('%H:%M') # Sets date format
ax1.xaxis.set_major_formatter(date_form)
# ~ plt.gca().invert_yaxis() # Invert axis for imshow

# Cosmetics
plt.title('CORAL lidar timeseries for 2019-07-26', fontdict = {'fontsize': 15})
ax1.set_ylabel('Altitude / km')
ax1.set_yticks(np.linspace(10,100,10))
ax1.yaxis.set_ticks_position('both')
ax1.tick_params(axis='y', which='both', labelleft='on', labelright='on')
ax1.set_xlabel('Time')
ax1.xaxis.grid(True, which='minor')
ax1.grid(color='gray', linestyle = 'dotted', linewidth = 0.25, axis='x',
			which='both')
ax1.grid(color='gray', linestyle = 'dotted', linewidth = 0.25, axis='y',
			which='both')
			
# Ensure the number of date ticks is sensible
ax1.xaxis.set_minor_locator(dates.MinuteLocator(byminute=0))
ax1.set_xlim([datetime(2019, 7, 26, 3, 0, 0), datetime(2019, 7, 26, 19, 0, 0)])

# Label fontsize
ax1.tick_params(labelsize=13)
ax1.yaxis.label.set_size(13)
ax1.xaxis.label.set_size(13)

# Add colorbar to figure
fig.subplots_adjust(bottom=0.3, right=0.88, left=0.12)
cbar_ax = fig.add_axes([0.12, 0.15, 0.76, 0.025])
fig.colorbar(cs, cmap=rdbucmap, orientation='horizontal',
	label='Temperature Perturbation / K', cax=cbar_ax, extend='both', boundaries=np.linspace(-22.5,22.5,19),
	ticks=np.linspace(-20,20,9))
cbar_ax.tick_params(labelsize=13)
cbar_ax.xaxis.label.set_size(13)
			
# Saving figure
plt.savefig('CORAL_timeseries_perturbations.png',dpi=600)
print("\nCORAL_timeseries_perturbations.png saved in: ", os.getcwd())
