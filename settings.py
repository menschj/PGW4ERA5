#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
description     Settings namelist for all routines in PGW for ERA5
authors		    Before 2022: original developments by Roman Brogli
                Since 2022:  upgrade to PGW for ERA5 by Christoph Heim 
"""
##############################################################################
##############################################################################

### GENERAL SETTINGS
##############################################################################
# debug output level
i_debug = 2 # [0-2]

# File naming convention for GCM climate deltas
# ({} is placeholder for variable name).
climate_delta_file_name_base = '{}_delta.nc'

# File naming convention for GCM files for the ERA climatology
# ({} is placeholder for variable name).
# This is required exclusively for the surface pressure (ps).
era_climate_file_name_base = '{}_historical.nc'

# File naming convention for ERA5 files to be read in and written out.
era5_file_name_base = 'cas{:%Y%m%d%H}0000.nc'
#era5_file_name_base = 'caf{:%Y%m%d%H}.nc'


### 02 PREPROCESS DELTAS 
##############################################################################

### SMOOTHING
####################################

### REGRIDDING
####################################
# depending on whether the xesmf pacakge is installed, it can be used
# for interpolation. Else, an xarray-based method is used.
# the latter should be identical to XESMF
# except for tiny differences that appear to originate from
# numerical precision
i_use_xesmf_regridding = 0


### SURFACE PRESSURE ADJUSTMENT SETTINGS 
##########################################################################

### SURFACE PRESSURE ADJUSTMENT SETTINGS 
##########################################################################
# re-interpolation turned on/off
i_reinterp = 0
# reference pressure level (None is stronlgy recommended)
#p_ref_inp = 50000
p_ref_inp = None
# surface pressure adjustment factor in the iterative routine
adj_factor = 0.95
# convergence threshold (maximum geopotential error)
# if procedure does not converge, raise this value a little bit.
thresh_phi_ref_max_error = 0.10
# maximum number of iterations before error is raised.
max_n_iter = 10
##########################################################################

# dimension names in ERA5 file
TIME_ERA        = 'time'
VERT_ERA        = 'level'
VERT_HL_ERA     = 'level1'
SOIL_HL_ERA     = 'soil1'
LON_ERA         = 'lon'
LAT_ERA         = 'lat'

# dimension names in GCM
TIME_GCM        = 'time'
PLEV_GCM        = 'plev'
LON_GCM         = 'lon'
LAT_GCM         = 'lat'

# map from naming convention variable names
# to variable names in ERA5 files to process.
var_name_map = {
    # surface geopotential
    'zgs'  :'FIS',
    # geopotential
    'zg'   :'PHI',
    # air temperature
    'ta'   :'T',   
    # near-surface temperature
    'tas'  :'T_SKIN',
    # soil temperature
    'st'   :'T_SO',
    # air relative humidity
    'hur'  :'RELHUM',
    # air specific humidity
    'hus'  :'QV',
    # lon-wind speed
    'ua'   :'U',
    # lat-wind speed 
    'va'   :'V',
    # surface pressure
    'ps'   :'PS',
}
