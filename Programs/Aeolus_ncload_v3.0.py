#!/usr/bin/env python3
"""
---------------------
Aeolus_ncload_v3.0.py
---------------------
Load and plot Aeolus data from netCDF format as a profile timeseries
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Looping through dataset and plotting all orbits that travel---
----------through the 40-80 box.----------------------------------------
---v1.2---Various pdates------------------------------------------------
---v1.3---Deals with multiple orbits------------------------------------
----------simple_contourf-----------------------------------------------
---v1.4---pcolormesh/imshow---------------------------------------------
---v1.5---Additional subplot with lat/lon grid showing satellite track--
---v1.6---More Updates: Applying S-G Filter and 2D Boxcar on data-------
---v1.7---Experimenting with 500m bins and interpolation, different-----
----------scheme for dealing with NaNs etc.-----------------------------
---v1.8---ERA5 interpolated onto Aeolus track---------------------------
---v1.9---Testing New NC_FullQC files, smoothing ERA5 first (in u and v)
---v1.10--Consolidating code and tidying up-----------------------------
---v2.0---Initial File for GRL paper------------------------------------
---v2.1---New colormaps, topography, more ERA5 fields etc.--------------
---v3.0---Final Published Version---------------------------------------
----------[CURRENT]-This_is_the_current_version_of_this_file------------
------------------------------------------------------------------------
========================================================================
Reads .nc files converted from DBL files from the Aeolus database and 
produces plots of key parameters from the datasets
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
from scipy.signal import savgol_filter
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

# Import colormaps
rdbucmap = LinearSegmentedColormap('RDBUcmap',
		   segmentdata=fetchcmap('RdBunew'), N=265)
grayhatchescmap = LinearSegmentedColormap('Grayhatchescmap',
				  segmentdata=fetchcmap('grayhatches'), N=265)

# matplotlib settings
matplotlib.rcParams['axes.unicode_minus'] = False

########################################################################
####### The settings in this section can be altered by the user ########
########################################################################
# Find directory and read Aeolus netCDF data
strdirectory = './AE_Data/'
# ERA5 directory
ERA5_dir = './ERA5_Data/'

# Choose Aeolus or ERA5. Options: 'ae', 'era5'
ae_era5 = 'ae'

# Choose Raw HLOS or Perturbations. Options: 'hlos', 'perts'
hlos_perts = 'perts'

# Choose pcolor or imshow
pc_or_im = 'im'
im_interp = 'none' # See documentation on imshow interpolations

# Maximum height of timeseries plot
maxh_km = 24 # Given in km
# Choose vertical resolution
vert_res = 500 # Given in metres

# Set oldinterp to True to use the old interpolation routine (nearest
# neighbour)
oldinterp = False

# Choice of region (Comment in the desired region) - Deprecated Tool
region = 'andes'
if region == 'andes':
	minlat, maxlat, minlon, maxlon = -80, -40, 280, 320
	bmlowlat, bmupperlat, bmleftlon, bmrightlon = -80, -40, -80, -40

# Choose ERA5 Overlay. Options: 'none', 'mslp', 'cloud', 'mslp+cloud'
era5_type = 'mslp+cloud'
########################################################################
########################################################################

# Program Timing (Time taken to get through entire program following
# toggles)	
pstartTime = datetime.now()

# Iterate through files in directory
directory = os.fsencode(strdirectory)
for file in os.listdir(directory):
	
	print("\n===========================================================")
	print("Beginning working directory:", os.getcwd())
	# Program Timing (Time taken to get through one file)	
	fstartTime = datetime.now()
	
	# Setting filename for dataload
	filename = os.fsdecode(file)
	print('Aeolus Data filename:', str(filename))
	infile = strdirectory + str(filename) # Specifies data file
	print('Aeolus netCDF file path:', infile, '\n')
	data = nc.Dataset(infile)
	
	"""=============================================================="""
	"""======================Download Variables======================"""
	"""=============================================================="""
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
	mnopib = 150 # minimum_number_of_profiles_in_box
	printfullarrays()
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

		"""=========================================================="""
		"""=====================ERA5 Interpolation==================="""
		"""=========================================================="""
		# =========================Load ERA5 Data=======================
		# Read data for the current date using Aeolus data filename
		YYYY = str(filename)[6:10]
		MM = str(filename)[11:13]
		DD = str(filename)[14:16]
		doy = yyyymmdd_to_doy(YYYY, MM, DD)
		ncfile = ERA5_dir + str(YYYY) + '/era5_' + str(YYYY) + 'd' \
		+ doy + '.nc'
		print('ERA5 netCDF file path:', ncfile)
		data = ERA5_dataload(ncfile)
		
		# Create special variables for ERA5 interpolation
		data_lon_new_ae = np.copy(data_lon_new)
		data_lat_new_ae = np.copy(data_lat_new)
		data_alt_new_ae = np.copy(data_alt_new)
		
		# Pass ERA5 data into their respective variables
		ERA5_data_lon, ERA5_data_lat, ERA5_data_lev, \
		ERA5_data_temp, ERA5_data_u, ERA5_data_v, \
		ERA5_data_time = data[0], data[1], data[2], data[3], \
		data[4], data[5], data[6]
		
		# Read data for the next date in case orbit crosses days
		nextday_dt = yyyymmdd_nextday(YYYY, MM, DD)
		YYYY2, MM2, DD2 = nextday_dt.year, nextday_dt.month, \
		nextday_dt.day
		doy2 = yyyymmdd_to_doy(YYYY2, MM2, DD2)
		# Specify ERA5 data file
		ncfile = ERA5_dir + str(YYYY2) + '/era5_' + str(YYYY2) + \
		'd' + doy2 + '.nc'
		print('Following ERA5 netCDF file path:', ncfile, '\n')
		data = ERA5_dataload(ncfile)
		
		# Pass ERA5 data for the next date into their respective
		# variables
		ERA5_data_temp, ERA5_data_u, ERA5_data_v, ERA5_data_time = \
		np.append(ERA5_data_temp, data[3], axis=0), \
		np.append(ERA5_data_u, data[4], axis=0), \
		np.append(ERA5_data_v, data[5], axis=0), \
		np.append(ERA5_data_time, data[6], axis=0)
		
		# Load any text files into arrays
		ERA5_p_levels = np.loadtxt("ERA5_pressure_levels.txt")*100
		ERA5_altitudes = np.loadtxt("ERA5_altitudes.txt")
		ERA5_levels = np.loadtxt("ERA5_pressure_levels_labels.txt")
		# ==============================================================
				
		# =======================Run Interpolation======================
		print('Running interpolation of ERA5 data...')
		# Create a datetime format version of rayleigh_times_new
		ERA5_u_interpolated = np.empty(0)
		ERA5_v_interpolated = np.empty(0)
		ERA5_interpolated = np.empty(0)
		no_result = 0
		rayleigh_times_new_dt = \
		coda.time_to_utcstring(rayleigh_times_new[:])
		rayleigh_times_new_dt = np.array([datetime.strptime(date,
			'%Y-%m-%d %H:%M:%S.%f') for date in rayleigh_times_new_dt])
		for t in range(len(rayleigh_times_new_dt)):
			# Finding the correct time array elements  
			rounded_down = round_down_3_hours(rayleigh_times_new_dt[t])
			rounded_up = round_up_3_hours(rayleigh_times_new_dt[t])
			down_elmt = np.where(ERA5_data_time == rounded_down)
			up_elmt = np.where(ERA5_data_time == rounded_up)
			
			for each in ['u', 'v', 'hlos']:
				# If 'hlos' then we have already got a u and v value so:
				if each == 'hlos':
					# Simple calculation for HLOS from u, v, and azimuth
					hlos = (-ERA5_u_interpolated[t] * np.sin(data_azimuth_new[t] * (np.pi/180))) - \
						(ERA5_v_interpolated[t] * np.cos(data_azimuth_new[t] * (np.pi/180)))
					ERA5_interpolated = np.append(ERA5_interpolated, hlos)
					# Now skip the rest of the loop
					continue
					
				# Toggle between ERA5 sandwich times at 3hr intervals
				# and interpolate between them: Down=Before|Up=After
				for t_E in [down_elmt, up_elmt]:
					t_E = t_E[0][0] # t_E = ERA5 time element
					if t_E == down_elmt[0][0]:
						time_toggle = 0
					elif t_E == up_elmt[0][0]:
						time_toggle = 1
					else:
						print("Error")
						
					# Call run_interpolation to execute 3D spatial
					# interpolation of ERA5 data onto Aeolus track
					if each == 'u':
						ERA5_var = ERA5_data_u
					elif each == 'v':
						ERA5_var = ERA5_data_v
					result = run_interpolation(data_lat_new_ae,
					data_lon_new_ae, data_alt_new_ae, 
					ERA5_data_lat, ERA5_data_lon, ERA5_altitudes, ERA5_var, t, t_E)
					# Raise flag if there is no result
					if result == None:
						no_result = 1
						continue
					# Set up and down times
					if time_toggle == 0:
						down = result
						# ~ print("Down = ", down[0])			
					elif time_toggle == 1:
						up = result
						# ~ print("Up = ", up[0])
					if no_result == 1:
						break
				
				# Break out of loop if there is no result		
				if no_result == 1:
					if each == 'u':
						ERA5_u_interpolated = np.append(ERA5_u_interpolated, 99999)
					elif each == 'v':
						ERA5_v_interpolated = np.append(ERA5_v_interpolated, 99999)
					no_result = 0
					continue
				
				# Linear interpolation in time using linearinterp2kp
				# function
				x_0 = toTimestamp(ERA5_data_time[down_elmt[0][0]])
				x_1 = toTimestamp(ERA5_data_time[up_elmt[0][0]])
				x = toTimestamp(rayleigh_times_new_dt[t])
				y_0 = down[0]
				y_1 = up[0]
				y = linearinterp2kp(x, x_0, x_1, y_0, y_1)
				
				if each == 'u':
					ERA5_u_interpolated = np.append(ERA5_u_interpolated, y)
				elif each == 'v':
					ERA5_v_interpolated = np.append(ERA5_v_interpolated, y)
					
		"""=========================================================="""
		"""=================Creating arrays for plot================="""
		"""=========================================================="""

		# Initialise meshgrids for x, y and z
		maxheight = maxh_km * 1000
		levnum = (maxheight / vert_res) + 1
		alts = np.linspace(0, maxheight, levnum)
		
		# Interpolation onto plotting grid for both ERA5 and observations
		data_obs = np.copy(data_HLOS_new) # Store observations for later use
		data_mod = np.copy(ERA5_interpolated)*100 # Converts ERA5 data into cm/s
		for each in ['mod', 'obs']:
			if each == 'mod':
				data_HLOS_new = data_mod
				print("\nInterpolating ERA5 onto plotting grid...\n")
			elif each == 'obs':
				data_HLOS_new = data_obs
				print("\nInterpolating Aeolus observations onto plotting grid...\n")
			
			# Initialising arrays for HLOS data
			z = np.zeros((len(alts),len(RG))) # NumPy Arrays
			z_itrn = np.zeros((len(alts),len(RG)))
			
			# ============OLD INTERPOLATION ROUTINE==============
			if oldinterp == True:
				# Placing wind values into bins of height 1km and width
				# 1 rayleigh group
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
							# Find the nearest altitude level
							val = find_nearest(alts, data_alt_new[t])
							alt_elmnt = np.where(alts == val)[0][0]
							# Cap wind speeds to 250 m/s
							if np.abs(data_HLOS_new[t]) < 25000:
								z[alt_elmnt][RG_elmnt] += data_HLOS_new[t]
								z_itrn[alt_elmnt][RG_elmnt] += 1
							RG_end = RG_elmnt
					lastgroupstarttime = RG[RG_elmnt]
			# ===================================================			
			else:
			# ============NEW INTERPOLATION ROUTINE==============
				z_new = np.zeros((len(alts),len(RG))) # NumPy Arrays
				locs = np.zeros((3, len(RG))) # For topography
				points = np.empty(0) # For AE/ERA5
				values = np.empty(0) # For AE/ERA5
				lons = np.empty(0) # For topography
				lats = np.empty(0) # For topography
				times = np.empty(0) # For topography
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
								# Apply QC for Aeolus data only
								if ae_era5 == 'ae':
									if data_QC_new[t] == 1:
										points = np.append(points, data_alt_new[t])
										values = np.append(values, data_HLOS_new[t])
								if ae_era5 == 'era5':
									if data_QC_new[t] == 1 or data_QC_new[t] == 0:
										points = np.append(points, data_alt_new[t])
										values = np.append(values, data_HLOS_new[t])
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
					points = np.empty(0)
					values = np.empty(0)
					lons = np.empty(0)
					lats = np.empty(0)
					times = np.empty(0)
			# ===================================================
								
			# Find the mean for each bin (Toggled on for old interpolation)
			if oldinterp == True:
				z /= 100 * z_itrn # Factor of 100 for conversion from cm/s - m/s
			else:
				# Using new interpolation routine
				z = np.copy(z_new) / 100
					
			# Amend RG array
			RG_new = RG[RG_start:RG_end+1]
			z = z[:, RG_start:RG_end+1]
			locs = locs[:, RG_start:RG_end+1]
			
			# Dealing with empty profiles at beginning and end
			for RG_i in [0, -1]:
				if not 0 not in z[:,RG_i]:
					z[:,RG_i] = np.NaN
			
			if each == 'mod':			
				data_mod = np.copy(z)
			elif each == 'obs':
				data_obs = np.copy(z)
		
		# Fixing time dimension to match coda format
		date_time = coda.time_to_utcstring(RG_new[:])
		date_time = np.array([datetime.strptime(date,
			'%Y-%m-%d %H:%M:%S.%f') for date in date_time])
		rayleigh_times_new = \
		coda.time_to_utcstring(rayleigh_times_new[:])
		rayleigh_times_new = np.array([datetime.strptime(date,
			'%Y-%m-%d %H:%M:%S.%f') for date in rayleigh_times_new])
		x, y = np.meshgrid(date_time, alts)
		
		# Applying NaN mask
		isnanarray = np.isnan(data_obs)
		hatcharray = np.copy(data_obs)
		for horidx in range(len(isnanarray)):
			for veridx in range(len(isnanarray[horidx])):
				if isnanarray[horidx][veridx] == True:
					hatcharray[horidx][veridx] = 1
				else:
					hatcharray[horidx][veridx] = 0
					
		# Adding topography time-series
		topodir = './elev.0.25-deg.nc'
		topo = nc.Dataset(topodir)	
		# Longitude
		topo_lon = topo.variables['lon'][:]
		# Latitude
		topo_lat = topo.variables['lat'][:]
		# Altitude
		topo_alt = topo.variables['data'][0]
		
		# Creating x_topo for topography plot
		x_topo = coda.time_to_utcstring(locs[2][:])
		x_topo = np.array([datetime.strptime(date,
			'%Y-%m-%d %H:%M:%S.%f') for date in x_topo])
		x_topo = dates.date2num(x_topo)
				
		# Creating y-topo for topography plot using topography_interpolation func
		y_topo = np.zeros(0)
		for profile in range(len(locs[0])):
			y_topo = np.append(y_topo,
				topography_interpolation(locs[0][profile], locs[1][profile], 
					topo_lon, topo_lat, topo_alt))
		
		# ================Topography and ERA5 data for Map==============
		# Creating topography for map
		x_map_topo = topo_lon[1119:1281]-360
		y_map_topo = topo_lat[519:681]
		z_map_topo = topo_alt[519:681, 1119:1281]
		x_map_topo, y_map_topo = np.meshgrid(x_map_topo, y_map_topo)
		
		# Fetching ERA5 file again for ERA5 map
		if era5_type in ['mslp', 'cloud', 'mslp+cloud']:
			ERA5_dir_andes = './ERA5_andessfcvars/'
			midpoint_time = coda.time_to_utcstring(RG_new[int(np.ceil(len(RG_new)/2))])
			myear, mmonth, mday, mtime = midpoint_time[0:4], midpoint_time[5:7], \
				midpoint_time[8:10], midpoint_time[11:16]
			print("Fetching ERA5 netCDF file for map for", myear, "-", mmonth, "-", mday, mtime)
			mdoy = yyyymmdd_to_doy(myear, mmonth, mday)
			ncfile = ERA5_dir_andes + myear + '/era5_' + myear + 'd' \
			+ mdoy + '_asv.nc'
			print('ERA5 netCDF file:', ncfile)

		# ERA5 Map
		import xarray as xr		
		ds = xr.open_dataset(ncfile)
		
		# ERA5_Original Code
		if era5_type == 'original':
			era5_data = {
				'data_t': ds.data_vars['t'][:],
				'data_u': ds.data_vars['u'][:],
				'data_v': ds.data_vars['v'][:]}	
			andes_slice = era5_data['data_u'].isel(level=slice(19,20)).isel(time=slice(0, 1)).isel(latitude=slice(86, 115)).isel(longitude=slice(66,95))
			z_era5_map = andes_slice.values[0][0]	
			x_era5_map = np.arange(-81, -37.5, 1.5)
			y_era5_map = np.arange(-81, -37.5, 1.5)
			y_era5_map = y_era5_map[::-1]
			x_era5_map, y_era5_map = np.meshgrid(x_era5_map, y_era5_map)
		
		# ERA5 Andes Surface variables (andessfvars)
		if era5_type in ['mslp', 'cloud', 'mslp+cloud']:
			# Hour beginning time for andessfvars slice
			hh = int(mtime[0:2])
			era5_data = {
				'data_anor': ds.data_vars['anor'][:], # Angle of sub-gridscale orography
				'data_isor': ds.data_vars['isor'][:], # Anisotropy of sub-gridscale orography
				'data_lgws': ds.data_vars['lgws'][:], # Eastward GW surface stress
				'data_gwd': ds.data_vars['gwd'][:], # GW dissipation
				'data_msl': ds.data_vars['msl'][:], # MSLP
				'data_mgws': ds.data_vars['mgws'][:], # Northward GW surface stress
				'data_z': ds.data_vars['z'][:], # Geopotential
				'data_slor': ds.data_vars['slor'][:], # Slope of sub-gridscale orography
				'data_sdor': ds.data_vars['sdor'][:], # Stddev of orography
				'data_tcc': ds.data_vars['tcc'][:]} # Cloud area fraction
			
			# ERA5 Map Timestamp
			tstamp = "ERA5 Timestamp: " + str(ds.coords['time'][hh].values)[11:16]
			print(tstamp, '\n')
			
			# Simulated IR Cloud Fraction
			andes_cloud = era5_data['data_tcc'].isel(time=slice(hh,hh+1)).isel(latitude=slice(160, 321)).isel(longitude=slice(40, 201))
			z_era5_cloud = andes_cloud.values[0]
			x_era5_cloud = np.arange(-80, -39.75, 0.25)
			y_era5_cloud = np.arange(-80, -39.75, 0.25)
			y_era5_cloud = y_era5_cloud[::-1]
			x_era5_cloud, y_era5_cloud = np.meshgrid(x_era5_cloud, y_era5_cloud)
			
			# MSLP contour plot
			andes_mslp = era5_data['data_msl'].isel(time=slice(hh,hh+1)).isel(latitude=slice(160, 321)).isel(longitude=slice(40, 201))
			z_era5_mslp = andes_mslp.values[0]
			x_era5_mslp = np.arange(-80, -39.75, 0.25)
			y_era5_mslp = np.arange(-80, -39.75, 0.25)
			y_era5_mslp = y_era5_mslp[::-1]
			x_era5_mslp, y_era5_mslp = np.meshgrid(x_era5_mslp, y_era5_mslp)
			
			andes_slice = era5_data['data_anor'].isel(time=slice(hh,hh+1)).isel(latitude=slice(160, 321)).isel(longitude=slice(40, 201))
			z_era5_map = andes_slice.values[0]
			x_era5_map = np.arange(-80, -39.75, 0.25)
			y_era5_map = np.arange(-80, -39.75, 0.25)
			y_era5_map = y_era5_map[::-1]
			x_era5_map, y_era5_map = np.meshgrid(x_era5_map, y_era5_map)
			
		# ==============================================================
		
		# Adding track axis on top in km
		km_counter = np.empty(0)
		dist = 0
		hvlon1 = np.deg2rad(locs[0][0]-360)
		hvlat1 = np.deg2rad(locs[1][0])
		for pr in range(len(locs[0])):
			hvlon2 = np.deg2rad(locs[0][pr]-360)
			hvlat2 = np.deg2rad(locs[1][pr])
			dist += haversine([hvlat1, hvlat2], [hvlon1, hvlon2])
			km_counter = np.append(km_counter, dist)
			hvlon1, hvlat1 = np.copy(hvlon2), np.copy(hvlat2)
		
		# Preparing arrays for plotting
		if ae_era5 == 'era5':
			z = data_mod
		fixnanswithmean(z) # Uses my own function 'fixnanswithmean'	
		
		if hlos_perts == 'perts':
			# ====================== S-G Filter Schemes ====================== #
			"""# Test using Gaussian noise (Enable by deleting quotation marks)
			for i in range(len(z)):
				for j in range(len(z[0])):
					z[i][j] = random.gauss(0, 10)"""
			
			z1 = np.copy(z)
			# Calculate appropriate upper and lower bounds to band-pass
			sg_upper = int(np.ceil(15 * (1000/vert_res)))
			sg_upper += (1-(sg_upper % 2))
			sg_lower = int(np.floor(5 * (1000/vert_res)))
			sg_lower += (1-(sg_lower % 2))
			try:
				z25 = savgol_filter(z, 25, 2, axis = 1) # Horizontal + S-G filter
				z7 = savgol_filter(z, 7, 2, axis = 1) # Horizontal - S-G filter
			except:
				# Return to programs directory
				print("Error: Savitzky-Golay filtering failed")
				# Finish program
				start_elmnt = 0
				end_elmnt = 0
				complete_boxes += 1
				continue
			z = z7 - z25
			# Boxcar smooth:
			# ~ z = ndimage.uniform_filter(z, size=(2,2), mode = 'reflect')
			# ================================================================= #
		
		# Setting limits for colorbar
		if hlos_perts == 'hlos':
			vminval = np.mean(z) - 2*np.std(z)
			vmaxval = np.mean(z) + 2*np.std(z)
		elif hlos_perts == 'perts':
			vminval, vmaxval = -10, 10 # Settings for S-G Perts	
		# ~ vminval, vmaxval = -20, 20 # Manually set vmin/vmax
		else:
			print("Error: Bad value for hlos_perts")
			sys.exit(1)
		
		# Apply NaN mask
		if ae_era5 == 'era5':
			for o in range(len(hatcharray)):
				for p in range(len(hatcharray[0])):
					if o > 1:
						hatcharray[o,p] = 0
					if o > 0 and p > 17:
						hatcharray[o,p] = 0

		"""=========================================================="""
		"""=======================Plotting==========================="""
		"""=========================================================="""
		print("Current working directory:", os.getcwd())
		os.chdir('Plots') # Change directory to Plots folder
		# Save data for MATLAB (toggle on and off)
		# ~ savemat("data.mat", {'data':z})
		# ~ np.save("data.npy", [x, y, z], allow_pickle=True)
		# ~ np.save("rayleigh_times.npy", rayleigh_times_new, allow_pickle=True)
		# ~ np.save("date_time.npy", date_time, allow_pickle=True)
		# ~ np.save("data_track.npy", [data_lat_new, data_lon_new], allow_pickle=True)
		
		YYYY = str(filename)[6:10]
		MM = str(filename)[11:13]

		# ==========================Plotting data============================ #
		# Initialise figure
		print('Initialise plot...')
		fig = plt.figure()
		gridspec.GridSpec(5,5) # Initialise gridspec
		
		# Main plot
		ax1 = plt.subplot2grid((5,5), (0,2), colspan=3, rowspan=4)
		
		# Flatten x_meshgrid and convert from datetime format to mpl_compatible nums
		x_lims = [np.ndarray.flatten(x)[0], np.ndarray.flatten(x)[-1]]
		x_lims = dates.date2num(x_lims)
		# Y limits
		y_lim_max = (maxheight + (vert_res/2)) / 1000
		y_lim_min = (0 - (vert_res/2)) / 1000
		y_lims = [y_lim_max, y_lim_min]
		
		# Choose pc or im
		if pc_or_im == 'pc':
			# Plot pcolormesh plot
			print('Plotting using pcolormesh...\n')
			cs = plt.pcolormesh(x, y/1000, z, cmap=rdbucmap, vmin=vminval, vmax=vmaxval)
			fixmydateaxis(ax1, date_time) # Uses my own function 'fixmydateaxis'
			
			# NaN mask
			mask = plt.pcolormesh(x, y/1000, hatcharray, cmap=grayhatchescmap)
		
		elif pc_or_im == 'im':
			# Plot imshow plot
			print('Plotting using imshow...\n')
			# Plot profiles using imshow
			cs = plt.imshow(z, aspect='auto', cmap=rdbucmap, extent=[x_lims[0],
				x_lims[1], y_lims[0], y_lims[1]], vmin=vminval, vmax=vmaxval,
				interpolation=im_interp)
			
			# NaN mask
			mask = plt.imshow(hatcharray, aspect='auto', cmap=grayhatchescmap,
				extent=[x_lims[0], x_lims[1], y_lims[0], y_lims[1]], interpolation=im_interp)
		
		else:
			print("\nError: Bad value for pc_or_im")
			sys.exit(1)
			
		# Finalising imshow plot
		ax1.set_xlim(x_lims[0], x_lims[1])
		ax1.set_ylim(y_lims[0], y_lims[1])
		ax1.xaxis_date() # Initialises date axis
		date_form = dates.DateFormatter('%H:%M') # Sets date format
		ax1.xaxis.set_major_formatter(date_form)
		plt.gca().invert_yaxis() # Invert axis for imshow
		
		# Topography
		topography = ax1.fill_between(x_topo, y_topo/1000, -10, color= '#987f52')
		topography_line = ax1.plot(x_topo, y_topo/1000, color = 'k', linewidth=0.7)
		
		# ================================================================= #

		# Fix axes

		### X axis
		ax1.set_xlabel('Time')
		# Ensure the number of date ticks is sensible
		timerange = date_time[-1] - date_time[0]
		date_intvl = int(np.ceil(timerange.seconds/(5*60)))
		ax1.xaxis.set_major_locator(plt.MaxNLocator(5)) # Maximum number of date ticks
		ax1.xaxis.set_major_locator(dates.MinuteLocator(interval=date_intvl))
		# ~ ax1.grid(color='k', linestyle = 'dashed', linewidth = 0.25, axis='x') # Enable vertical grid lines
		ax1b = ax1.twiny()
		ax1b.xaxis.tick_top()
		ax1b_end = km_counter[-1] + (km_counter[-1] - km_counter[-2])/2
		ax1b.set_xticks(np.arange(0,ax1b_end,500))
		ax1b.set_xlim(0, ax1b_end)
		ax1b.set_xlabel('Along Track Distance / km', labelpad = 7)

		### Y axis
		ax1.set_ylabel('Altitude / km')
		ax1.set_yticks(np.arange(np.ceil(maxh_km)+1))
		ax1.yaxis.set_major_locator(plt.MaxNLocator(np.ceil(maxh_km)/2 + 1))
		ax1.yaxis.set_minor_locator(plt.MaxNLocator(np.ceil(maxh_km)+1))
		ax1.tick_params(axis='y', which='both', left=True, right=True, labelright=True)
		ax1.grid(color='gray', linestyle = 'dotted', linewidth = 0.25, axis='y',
			which='both')
		ax1.yaxis.set_label_position("right")

		# ======================= Satellite Track plot ======================= #
		ax2 = plt.subplot2grid((5,5), (0,0), colspan=2, rowspan=4)
		map = Basemap(projection='cyl',llcrnrlat=bmlowlat,urcrnrlat=bmupperlat,\
					llcrnrlon=bmleftlon,urcrnrlon=bmrightlon,resolution='i', ax=ax2)
		map.fillcontinents(color='#ffdd99', lake_color='#cceeff', zorder = 2, alpha=1)
		map.drawmapboundary(linewidth=0.75, fill_color='#cceeff')
		map.drawcoastlines(linewidth=0.25, color='#666666', zorder=5)
		if region == 'andes':
			map.drawmeridians([-70, -60, -50], linewidth=0.3, zorder=5)
			map.drawparallels([-70, -60, -50], linewidth=0.3, zorder=5)
		elif region == 'azores':
			map.drawmeridians([-60, -40, -20], linewidth=0.3, zorder=5)
			map.drawparallels([20, 30, 40, 50, 60], linewidth=0.3, zorder=5)
		
		# Map topography
		topo_levels = [25.100,200,500,1000,1500,2000,2500]
		trial_levels = np.geomspace(150, 3100, num=15)
		map_topo_iso = map.contour(x_map_topo, y_map_topo, z_map_topo, levels=topo_levels, colors='k', linewidths = 0.15, zorder=4, alpha = 0.4)
		map_topo = map.contourf(x_map_topo, y_map_topo, z_map_topo, levels=topo_levels, cmap='copper_r', zorder=3, alpha = 1)
		map_topo2 = map.contourf(x_map_topo, y_map_topo, z_map_topo, levels=topo_levels, cmap='Oranges', zorder=3, alpha = 0.4)
		map_topo3 = map.contourf(x_map_topo, y_map_topo, z_map_topo, levels=topo_levels, cmap='Greys', zorder=3, alpha = 0.25)
		
		if era5_type in ['cloud', 'mslp+cloud']:
			# Simulated IR Cloud Fraction
			era5_cloud = map.contourf(x_era5_cloud, y_era5_cloud, z_era5_cloud*100, cmap='Greys_r',
				zorder=3, alpha=0.35, levels=[0,20,40,60,80,100])
		
		if era5_type in ['mslp', 'mslp+cloud']:
			# MSLP plot
			small_mslp_levels = np.arange(920, 1060, 2)
			large_mslp_levels = np.arange(920, 1060, 10)
			for i in large_mslp_levels: # Remove multiples of ten for better isobars
				small_mslp_levels = np.delete(small_mslp_levels, np.where(small_mslp_levels == i))
			small_era5_mslp = map.contour(x_era5_mslp, y_era5_mslp, z_era5_mslp/100, levels=small_mslp_levels, linewidths=0.2, colors='b', zorder=5, alpha = 1)
			large_era5_mslp = map.contour(x_era5_mslp, y_era5_mslp, z_era5_mslp/100, levels=large_mslp_levels, linewidths=0.5, colors='b', zorder=5, alpha = 1)
			ax2.clabel(large_era5_mslp, large_era5_mslp.levels, inline=True, fmt= '%d', fontsize=5)
		
		# Track Plot
		map.scatter(data_lon_new - 360, data_lat_new, marker = 'x', color = 'red',
			s=0.3, zorder=6)
		map.plot(data_lon_new - 360, data_lat_new, color = 'red', linewidth = '0.5', zorder=6)
		midpointlonindex = int(np.floor(len(data_lon_new)/2))
		midpointlatindex = int(np.floor(len(data_lat_new)/2))

		# Add lines at time intervals corresponding to date_time ticks
		date_ticks = dates.num2date(ax1.xaxis.get_ticklocs())
		for i in range(len(date_ticks)):
			# Co-locate between ticks and rayleigh_times array
			nearest_date = find_nearest(rayleigh_times_new, date_ticks[i])
			idx = (np.where(dates.date2num(rayleigh_times_new) == nearest_date)[0])[0]
			str_time = date_ticks[i].strftime('%H:%M') # Sets string form of each tick
			# Location of each line
			lon_loc = data_lon_new[idx]-360
			lat_loc = data_lat_new[idx]
			# ~ map.scatter([lon_loc], [lat_loc], s=50, marker='o', zorder=111,
				# ~ color='black') # Large circles for timestamps rather than quivers
			# Plots line on map and annotates with text
			try: 
				data_lat_new[idx+30] == False
				data_lat_new[idx-30] == False
				data_lon_new[idx+30] == False
				data_lon_new[idx-30] == False
			except:
				continue
			map.quiver(lon_loc, lat_loc, data_lat_new[idx-30] - data_lat_new[idx+30],
				data_lon_new[idx+30] - data_lon_new[idx-30], angles='xy', color='k',
				zorder=7, width = 0.5, headaxislength=0, headlength=0, pivot='middle', units='xy')
			ax2.annotate(str_time, (lon_loc+1, lat_loc+1), fontsize=8, zorder=7, color = 'k', weight='demibold')
		
		if era5_type in ['mslp', 'cloud', 'mslp+cloud']:
			# Adding ERA5 timestamp to plot
			ax2.text(-40,-40, tstamp, size=8, color = 'blue', ha="right", va="bottom", zorder = 10)
		
		# ==================================================================== #

		# Fix axes
		if region == 'andes':
			ax2.set_xticks([-80, -70, -60, -50, -40])
			ax2.set_yticks([-80, -70, -60, -50, -40])
		elif region == 'azores':
			ax2.set_xticks([-80, -60, -40, -20, 0])
			ax2.set_yticks([10, 20, 30, 40, 50, 60, 70])
		ax2.set_xlabel('Longitude / deg')
		ax2.set_ylabel('Latitude / deg')
		ax2.set_aspect('auto') # Stretch map to fill subplot
		for axis in ['top','bottom','left','right']: # Set axes thickness
			ax1.spines[axis].set_linewidth(0.75)
			ax1.spines[axis].set_zorder(10)
			ax2.spines[axis].set_linewidth(0.75)
			ax2.spines[axis].set_zorder(10)
		
		# Add colorbar to Aeolus figure
		fig.subplots_adjust(bottom=0.2, right=0.88, left=0.12, wspace=0.55)
		cbar_ax = fig.add_axes([0.4475, 0.2, 0.4325, 0.025])
		if hlos_perts == 'hlos':
			fig.colorbar(cs, cmap=rdbucmap, orientation='horizontal', # RAW
				label='HLOS Rayleigh Wind Speed / ms$^{-1}$', cax=cbar_ax)
		if hlos_perts == 'perts':
			fig.colorbar(cs, cmap=rdbucmap, orientation='horizontal', # S-G Perts
				label='HLOS Rayleigh Wind Perturbation / ms$^{-1}$', cax=cbar_ax)
		
		# Add colorbar to Track figure
		if era5_type in ['cloud', 'mslp+cloud']:
			cbar_ax = fig.add_axes([0.12, 0.2, 0.27, 0.025])
			ax3 = fig.colorbar(era5_cloud, cmap='Greys_r', orientation='horizontal', # S-G Perts
				label='Cloud fraction / %', cax=cbar_ax,
				ticks=[0,20,40,60,80,100], spacing='uniform')
			ax3.set_label('Cloud fraction / %', labelpad = 5)

		# Set figure title
		print("\nFilename check:", str(filename)[6:-3])
		strdate = str(filename)[6:16]
		print("strdate:", strdate)
		strtime = str(filename)[17:19] + ':' + str(filename)[19:21] + ':' + \
			str(filename)[21:23]
		print("strtime:", strtime)
		if ae_era5 == 'ae':
			if hlos_perts == 'hlos':
				str_plt_title = 'Aeolus: HLOS Rayleigh Wind Cross-section'
			elif hlos_perts == 'perts':
				str_plt_title = 'Aeolus: HLOS Rayleigh Wind Perturbations'
		elif ae_era5 == 'era5':
			if hlos_perts == 'hlos':
				str_plt_title = 'ERA5 interpolated onto Aeolus: HLOS Rayleigh Wind Cross-section'
			elif hlos_perts == 'perts':
				str_plt_title = 'ERA5 interpolated onto Aeolus: HLOS Rayleigh Wind Perturbations'
		else:
			print("\nError: Bad value for ae_era5")
			sys.exit(1)
		str_plt_title += '\n' + 'Orbit: ' + strdate + ' ' + strtime
		plt.suptitle(str_plt_title, y=0.09, ha='center', fontweight='semibold')

		# Saving figure
		pngsavename = str(filename)[:-3]
		if complete_boxes != 0:
			pngsavename += '_orb' + str(complete_boxes+1)
		pngsavename += '.png'
		print("\npngsavename:", pngsavename, '\n\nSaving figure...')
		plt.savefig(pngsavename, dpi=600)
		print(pngsavename, "saved in:", os.getcwd())
		
		# Return to programs directory
		os.chdir('..')
		os.chdir('Programs')
		
		# Time taken for the file
		fduration = datetime.now() - fstartTime
		print('\nThat file took', fduration, 'seconds to analyse')
		
		# Reset after box completed
		start_elmnt = 0
		end_elmnt = 0
		complete_boxes += 1
		
		plt.close()

# Time taken for the file
pduration = datetime.now() - pstartTime
print('\nThat program took', pduration, 'seconds to run')
