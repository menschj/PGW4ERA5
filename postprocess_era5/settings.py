# debug output level
i_debug = 2 # [0-2]

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
