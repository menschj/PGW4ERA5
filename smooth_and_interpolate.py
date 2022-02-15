# -*- coding: utf-8 -*- 
###############################################################################
import subprocess, os, argparse
import xarray as xr
#import xesmf as xe
from interpolate import interpannualcycle
#from era5_vertical import vertinterpol
from package.utilities import cd, cdo_remap
from pathlib import Path
###############################################################################

## input arguments
parser = argparse.ArgumentParser(description =
                'Prepare PGW climate deltas for adding to BC files.')
# variable to plot
parser.add_argument('var_names', type=str,
                    help='separate multiple with ","')
# number of parallel processes
parser.add_argument('-p', '--n_par', type=int, default=1)
# number of BC files per day
parser.add_argument('-n', '--n_per_day', type=int, default=8)
args = parser.parse_args()
print(args)

performsmooth   = 1
regridhorinew   = 0

# number of output time steps in total
n_out_time_steps = 366 * args.n_per_day

var_name_map = {
        'hur'   :'RELHUM',
        'hus'   :'QV',
        'hurs'  :'RELHUM_S',
        'ta'    :'T', 
        'tas'   :'T_S', 
        'ps'    :'PS', 
        'ua'    :'U', 
        'va'    :'V',
        'zg'    :'PHI',
        'orog'  :'HSURF',
        }
var_names = args.var_names.split(',')
for var_name in var_names:
    if var_name not in var_name_map:
            raise ValueError('Input var_name {} is invalid.'.format(
                            var_name))


if performsmooth:
#settings for timesmoothing script:

    #list of the files that contain the variables to be smoothed
    #samepath = '/project/pr94/robro/inputfiles_for_surrogate_hadgem/input_github/'
    samepath = '/scratch/snx3000/heimc/pgw/deltas/day/MPI-ESM1-2-HR/'
    annualcycleraw = [
    samepath+'tas_delta.nc',
    #samepath+'Diff_HadGEM2-ES_RCP85_PP.nc',
    #samepath+'Diff_HadGEM2-ES_RCP85_QV.nc',
    #samepath+'Diff_HadGEM2-ES_RCP85_QV_S.nc',
    #samepath+'Diff_HadGEM2-ES_RCP85_T.nc',
    #samepath+'Diff_HadGEM2-ES_RCP85_T_S.nc',
    #samepath+'Diff_HadGEM2-ES_RCP85_T_SO.nc',
    #samepath+'Diff_HadGEM2-ES_RCP85_U.nc',
    #samepath+'Diff_HadGEM2-ES_RCP85_V.nc',
    #samepath+'Diff_HadGEM2-ES_RCP85_RELHUM.nc',
    #samepath+'Diff_HadGEM2-ES_RCP85_RELHUM_S.nc'
    ]
    #list of variablenames
    #variablename_to_smooth = ['PP', 'QV', 'QV_S', 'T', 'T_S', 'T_SO', 'U', 'V','RELHUM','RELHUM_S']
    variablename_to_smooth = ['tas']
    #path to put the output netcdf
    outputpath = '/scratch/snx3000/heimc/pgw/test_smoothing/'


    #enter the command to run the script:
    for num,pathi in enumerate(annualcycleraw):
            commandsmooth = f"python timesmoothing.py {pathi} {variablename_to_smooth[num]} {outputpath} > outputfile_smooth.txt &"
            subprocess.run(commandsmooth, shell=True)
            print(commandsmooth)
    

if regridhorinew:
    ###########################################################################
    ### Namelist
    ###########################################################################
    gcm_data_path='/scratch/snx3000/heimc/pgw/deltas/Emon/MPI-ESM1-2-HR'
    gcm_data_path='/scratch/snx3000/heimc/pgw/deltas/Emon_RHint/MPI-ESM1-2-HR'
    #gcm_data_path='/scratch/snx3000/heimc/pgw/deltas/Amon/MPI-ESM1-2-HR'
    delta_inp_name_base='{}_delta.nc'
    delta_out_name_base='{}_delta.nc'

    #delta_inp_name_base='{}_historical.nc'
    #delta_out_name_base='{}_historical.nc'

    #delta_inp_name_base='{}_ssp585.nc'
    #delta_out_name_base='{}_ssp585.nc'

    #gcm_data_path='/scratch/snx3000/heimc/pgw/deltas/test2/MPI-ESM1-2-HR'
    #delta_name_base='plev_{}_delta.nc'

    #out_dir = '/scratch/snx3000/heimc/pgw/regridded/Emon/extpar_SA_3'
    #out_dir = '/scratch/snx3000/heimc/pgw/regridded/Emon_RHint/MPI-ESM1-2-HR'
    out_dir = '/scratch/snx3000/heimc/pgw/regridded/Emon_xr/MPI-ESM1-2-HR'
    #out_dir = '/scratch/snx3000/heimc/pgw/regridded/Amon/MPI-ESM1-2-HR'
    #out_grid_file = 'target_grid_SA_3'
    #out_grid_file = 'target_grid_extpar_SA_3'
    #out_grid_file = 'target_grid_era5'
    #out_grid_file = 'target_grid_era5_test'
    #out_grid_file = 'target_grid_era5_test2'
    #out_grid_file = 'target_grid_era5_test3'
    target_file_path = '/scratch/snx3000/heimc/lmp/wd/06080100_SA_3_ctrl/int2lm_in/cas20060801000000.nc'
    ###########################################################################

    Path(out_dir).mkdir(parents=True, exist_ok=True)

    #for var_name in var_names:
    #    print('regrid horizontal {}'.format(var_name))
    #    inp_file = os.path.join(gcm_data_path, delta_inp_name_base.format(var_name))
    #    out_file = os.path.join(out_dir, delta_out_name_base.format(var_name))
    #    #print(inp_file)
    #    #print(out_file)
    #    cdo_remap(out_grid_file, inp_file, out_file, method='bil')

    #    #target_file = os.path.join(target_file_path, 'cas20060805000000.nc')
    #    #ds_gcm = xr.open_dataset(inp_file)
    #    #ds_target = xr.open_dataset(target_file)
    #    #regridder = xe.Regridder(ds_gcm, ds_target, "bilinear")
    #    ##dr_out = regridder(dr)
    #    #quit()

    for var_name in var_names:
        print('regrid horizontal {}'.format(var_name))
        inp_file = os.path.join(gcm_data_path, delta_inp_name_base.format(var_name))
        out_file = os.path.join(out_dir, delta_out_name_base.format(var_name))
        #print(inp_file)
        #print(out_file)
        targ_lon = xr.open_dataset(target_file_path).lon
        targ_lat = xr.open_dataset(target_file_path).lat
        ds_in = xr.open_dataset(inp_file)
        ds_out = ds_in.interp(lon=targ_lon, lat=targ_lat)
        ds_out.to_netcdf(out_file)

        #test = xr.open_dataset('/scratch/snx3000/heimc/pgw/regridded/Emon/MPI-ESM1-2-HR/tas_delta.nc')
        #test = test.assign_coords({'lat':ds_out.lat.values})
        #test = test - ds_out
        #test.to_netcdf('test.nc')
        #quit()


