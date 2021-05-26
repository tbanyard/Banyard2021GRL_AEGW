========================================================================
This GIT repository contains all the files required for producing the
plots in the following paper:
------------------------------------------------------------------------
Banyard, T. P. et al. (2021), Atmospheric Gravity Waves in Aeolus Wind 
Lidar Observations, Geophysical Research Letters
------------------------------------------------------------------------
Aeolus data were provided by the European Space Agency, and can be
accessed via https://aeolus-ds.eo.esa.int/oads/access/. The AIRS data
were provided by NASA; L1 radiance data can be acquired via
https://disc.gsfc.nasa.gov/, and were retrieved to L2 temperatures using
the method described by Hoﬀmann and Alexander (2009). ERA5 data can be
accessed from the Copernicus Climate Data Store,
https://cds.climate.copernicus.eu/. CORAL Lidar data are not routinely
archived publically, but the data used for this study have been archived
at https://halo-db.pa.op.dlr.de/dataset/7620.
========================================================================

========================================================================
------------------------------IMPORTANT---------------------------------
========================================================================
In order for these programs to work, the user will need to merge the
following data files into their GIT repository. Caution: These files
are very large, necessitating the use of the below link.
------------------------------------------------------------------------
tinyurl.com/DataforBanyard2021GRL
========================================================================


========================================================================
The following files and directories can be found in this repository:
------------------------------------------------------------------------
AE_Data| Holds a netCDF file which contains the Aeolus L2B HLOS Rayleigh
|||||||| Wind Data used in this study.
AIRS_Data| Holds netCDF and MATLAB files which contain the AIRS data
|||||||||| used in this study.
CORAL_Data| Holds a netCDF file which contains the CORAL data used in
||||||||||| this study.
ERA5_andessfcvars| Holds a netCDF file which contains ERA5 data on a
|||||||||||||||||| single level, used for cloud fraction and MSLP in the
|||||||||||||||||| 1st and 3rd figures of the GRL paper.
ERA5_Data| Holds 2 netCDF files which contain the ERA5 data for the day
|||||||||| of and day after the Aeolus overpass used in this case study.
ERA5_special| Holds a netCDF file which contains ERA5 data for the
||||||||||||| specific time snapshots used in the profiles plot in the
||||||||||||| 2nd figure of the GRL paper.
Plots| Holds the plots created from python programs in the Programs 
|||||| folder, as well as the figures from the research paper itself.
Programs| Holds all programs used in this study.
anacondaenv.yaml| Anaconda python environment file to run programs in.
elev.0.25-deg.nc| netCDF file containing TBASE topography data.
ERA5_altitudes.txt| List of altitudes corresponding to each ERA5 model
||||||||||||||||||| level, used as an approximation for interpolating
||||||||||||||||||| ERA5 data.
ERA5_pressure_levels.txt| Pressure levels corresponding to each ERA5
||||||||||||||||||||||||| model level, as above.
ERA5_pressure_levels_labels| List of model level labels.
README.md| This file.
========================================================================

========================================================================
Details for each python program can be found below:
________________________________________________________________________
AEGW_functions.py| Contains all global functions required in this study
Aeolus_ncload_v3.0.py| Produces Figs 1 & 3 in the paper
AIRS_ncload_v2.0.py| Produces AIRS plots in Fig 2 in the paper
cmapstore.py| Contains customized colormaps used in the paper
CORAL_ncload_v2.0.py| Produces CORAL plots in Fig 2 in the paper
profilesplot_v3.py| Produces the Profile plots in Fig 2 in the paper
========================================================================

========================================================================
Individual details for each can be found below:
________________________________________________________________________


========================================================================
AEGW_functions.py
========================================================================
------------------------------------------------------------------------
------------------------No version control------------------------------
------------------------------------------------------------------------
========================================================================
Provides functions for all files for the Aeolus project
------------------------------------------------------------------------
This program contains all the functions that are called in each plotting
program within this project, as well as some additional useful tools.
Listed in order, there are a variety of functions for dealing with
datasets in general, functions for dealing with the specific datasets
used in this study - specifically the ERA5, Aeolus and TBASE datasets,
a function for navigating directories in python, functions for plotting,
and functions for loading and converting Aeolus data from its original
file format.
========================================================================


========================================================================
Aeolus_ncload_v3.0.py
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
x--v3.0---Final Published Version--------------------------------------x
------------------------------------------------------------------------
========================================================================
Reads .nc files converted from DBL files from the Aeolus database and 
produces plots of key parameters from the datasets
------------------------------------------------------------------------
This program produces plots corresponding to Figures 1 and 3 in the GRL
paper. Aeolus data has been extracted from the original DBL files and
saved in netCDF format for each individual orbit. This program produces
a timeseries of profiles from Aeolus over the Southern Andes for an
orbit on the 26th July 2019, and calculates the wind perturbations from
the HLOS Rayleigh Wind product. This program can be run in an Aeolus or
ERA5 mode, and a raw HLOS wind or perturbations mode. The user can alter
these settings within the code.
========================================================================


========================================================================
AIRS_ncload_v2.0.py
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Additional Flexibility and uses the Lambert conformal conic---
----------projection.
---v1.2-----------------------------------------------------------------
---v1.3---Fixed for final GRL submission--------------------------------
x--v2.0---Final Published Version--------------------------------------x
------------------------------------------------------------------------
========================================================================
Reads .mat files of AIRS data and plots along scan track over the Andes
------------------------------------------------------------------------
This program produces plots of the AIRS temperature perturbations at an
altitude of 30 km, corresponding to the Aeolus orbit over the Andes of
26th July 2019. Two MATLAB files are loaded for a before and after
picture of the gravity wave activity around the time of the Aeolus orbit
for validation purposes.
========================================================================


========================================================================
cmapstore.py
========================================================================
------------------------------------------------------------------------
----------------------------No version control--------------------------
------------------------------------------------------------------------
========================================================================
Contains customized colormaps used in the Aeolus GRL paper.
------------------------------------------------------------------------
This program contains a function which is called to load a selection of
customized colormaps which are used in the Aeolus GRL paper. The RdBunew
colormap is used in the figures in the paper.
========================================================================


========================================================================
CORAL_ncload_v2.0.py
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v1.1---Timeseries using qbocmap and raw data-------------------------
---v1.2---Improving plot------------------------------------------------
---v1.3---Time domain flexibility and improved appearance for fig2 GRL--
x--v2.0---Final Published Version--------------------------------------x
------------------------------------------------------------------------
========================================================================
Reads .nc files of data from the CORAL lidar in Argentina, and plots it.
------------------------------------------------------------------------
This program produces a set of timeseries from the CORAL lidar which
correspond to the Aeolus orbit over the Andes of 26th July 2019. A plot
of the evolution of the raw temperatures is produced, as well as the
temperature perturbations. Only the latter appears in the Aeolus GRL
paper.
========================================================================


========================================================================
profilesplot_v3.py
========================================================================
------------------------------------------------------------------------
---v1.0---Initial_File--------------------------------------------------
---v2.0---Updated for GRL submission: Fixed AIRS profiles and axes------
x--v3.0---Final Published Version--------------------------------------x
------------------------------------------------------------------------
========================================================================
Reads data files for CORAL, AIRS, AEOLUS and ERA5 profiles figure.
------------------------------------------------------------------------
This program produces a triplet of profiles for before, during and after
the Aeolus orbit over the Andes of 26th July 2019. These profiles are
co-located in AIRS and ERA5 to the Aeolus orbit marked with an arrow in
figures 1 and 3, and a star in figure 2b and 2c, in the GRL paper.
========================================================================
________________________________________________________________________
Any questions regarding this GIT repository should be directed to the
corresponding author, Timothy P. Banyard, at tpb38@bath.ac.uk.
________________________________________________________________________

END
