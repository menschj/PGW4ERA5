#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
description     PGW for ERA5 preprocessing of climate deltas
authors		Before 2022: original developments by Roman Brogli
                Since 2022:  upgrade to PGW for ERA5 by Christoph Heim 
"""
##############################################################################
import os, argparse
import xarray as xr
import numpy as np
from pathlib import Path
from functions import filter_data, regrid_lat_lon
from settings import (
    i_debug,
    era_climate_file_name_base,
    climate_delta_file_name_base,
    LON_ERA, LAT_ERA,
    TIME_GCM, PLEV_GCM, LON_GCM, LAT_GCM,
    i_use_xesmf_regridding,
)
##############################################################################

## input arguments
parser = argparse.ArgumentParser(description =
            'PGW for ERA5: Preprocess climate deltas before adding ' +
            'them to ERA5. The main routine (step_03_apply_to_era.py) ' +
            'requires ' +
            'climate deltas for [ta,hur,ua,va,zg,hurs,tas], as well ' +
            'as the ERA climatological value of ps. The script should thus ' +
            'be run for variables ps,ta,hur,ua,va,zg,hurs,tas and it will ' +
            'automatically look for the file ' +
            'ps_$era_climate_file_name_base.nc ' +
            'and e.g. ta_$climate_delta_file_name_base.nc, where ' +
            '$era_climate_file_name_base and $climate_delta_file_name_base ' +
            'can be set in settings.py, among other things.')

# processing step to perform during script execution
parser.add_argument('processing_step', type=str, 
            choices=['smoothing','regridding'],
            help='Possible processing steps are: ' +
            'smoothing: [For daily climate deltas, a smoothing of ' +
            'the annual cycle should be applied. For monthly ' +
            'climate deltas this is not necessary.] ' +
            'regridding: [If the climate deltas are not on the same ' +
            'horizontal grid as ERA5, they can be regridded here. '+
            'WARNING: The script assumes regular (thus non-rotated) ' +
            'lat/lon grid for ' +
            'input (climate delta) and output (ERA5)!]')

# variable(s) to process
parser.add_argument('var_names', type=str,
            help='Variable names (e.g. tas) to process. Separate ' +
            'multiple variable names with "," (e.g. tas,ta).')

# input directory
parser.add_argument('-i', '--input_dir', type=str,
            help='Directory with input climate delta ' +
            'files for selected ' +
            'processing step.')

# output directory
parser.add_argument('-o', '--output_dir', type=str,
            help='Directory where output climate delta files of selected ' +
            'processing step should be stored.')

# target ERA5 example file to take grid information
parser.add_argument('-e', '--era5_file_path', type=str, default=None,
            help='Path to example ERA5 file ' +
            'from which to take grid information for regridding.')

args = parser.parse_args()
print(args)
##############################################################################

# make sure required input arguments are set.
if args.input_dir is None:
    raise ValueError('Input directory (-i) is required.')
if args.output_dir is None:
    raise ValueError('Output directory (-o) is required.')
if (args.processing_step == 'regridding') and (args.era5_file_path is None):
    raise ValueError('era5_file_path is required for regridding step.')

# create output directory
Path(args.output_dir).mkdir(exist_ok=True, parents=True)

# set up list of variable names
var_names = args.var_names.split(',')
print('Run {} for variable names {}.'.format(
        args.processing_step, var_names))


##############################################################################
for var_name in var_names:
    print(var_name)

    # for ps, take ERA climate mean value (i.e. e.g. ps_historical.nc)
    if var_name == 'ps':
        var_file_name = era_climate_file_name_base.format(var_name)
        #var_file_name = climate_delta_file_name_base.format(var_name)
    # for all other variables, take climate delta (i.e. e.g. tas_delta.nc)
    else:
        var_file_name = climate_delta_file_name_base.format(var_name)

    inp_file = os.path.join(args.input_dir, var_file_name)
    out_file = os.path.join(args.output_dir, var_file_name)


    ## smoothing
    if args.processing_step == 'smoothing':

        filter_data(inp_file, var_name, out_file)

    ## regridding
    elif args.processing_step == 'regridding':
        
        # TODO debug
        #ds_with = xr.open_dataset(
        #    '/net/o3/hymet_nobackup/heimc/data/pgw/local/regridded/with_xesmf/tas_delta.nc')
        #ds_without = xr.open_dataset(
        #    '/net/o3/hymet_nobackup/heimc/data/pgw/local/regridded/without_xesmf/tas_delta.nc')
        #print(ds_without - ds_with)
        #(ds_without - ds_with).to_netcdf('test.nc')
        #quit()
        # TODO debug

        ds_gcm = xr.open_dataset(inp_file)
        ds_era5 = xr.open_dataset(args.era5_file_path)

        ds_gcm = regrid_lat_lon(ds_gcm, ds_era5, var_name,
                                method='bilinear',
                                i_use_xesmf=i_use_xesmf_regridding)

        ds_gcm.to_netcdf(out_file)

