#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
description     PGW for ERA5
authors		    Before 2022: original developments by Roman Brogli
                After 2022:  updates by Christoph Heim 
"""
##############################################################################
import os, argparse
import xarray as xr
from pathlib import Path
from functions import filter_data
from settings import *
##############################################################################

## input arguments
parser = argparse.ArgumentParser(description =
            'PGW for ERA5: Preprocess climate deltas before adding ' +
            'them to ERA5. The main routine (pgw_for_era5.py) requires ' +
            'climate deltas for [ta,hur,ua,va,zg,hurs,tas], as well '+
            'as the ERA climatological value of ps.')

# processing step to perform during script execution
parser.add_argument('processing_step', type=str, 
            choices=['smoothing','regridding'],
            help='Possible processing steps are: ' +
            'smoothing: [For daily climate deltas, a smoothing of ' +
            'the annual cycle should be applied. For monthly ' +
            'climate deltas this is not necessary.] ' +
            'regridding: [If the climate deltas are not on the same ' +
            'horizontal grid as ERA5, they can be regridded here. '+
            'WARNING: The script assumes regular lat/lon grid for '+
            'input (climate delta) and output (ERA5)!]')

# variable(s) to process
parser.add_argument('var_names', type=str,
            help='Variable names (e.g. tas) to process. Separate ' +
            'multiple variable names with "," (e.g. tas,ta).')

# input directory
parser.add_argument('-i', '--input_directory', type=str,
            help='Directory with input files for selected ' +
            'processing step.')

# output directory
parser.add_argument('-o', '--output_directory', type=str,
            help='Directory where output files of selected ' +
            'processing step should be stored.')

# target ERA5 example file to take grid information
parser.add_argument('-e', '--era5_file_path', type=str, default=None,
            help='Path to example ERA5 file ' +
            'from which to take grid information for regridding.')

args = parser.parse_args()
print(args)

# make sure required input arguments are set.
if args.input_directory is None:
    raise ValueError('Input directory (-i) is required.')
if args.output_directory is None:
    raise ValueError('Output directory (-o) is required.')
if (args.processing_step == 'regridding') and (args.era5_file_path is None):
    raise ValueError('era5_file_path is required for regridding step.')

# create output directory
Path(args.output_directory).mkdir(exist_ok=True, parents=True)

# set up list of variable names
var_names = args.var_names.split(',')
print('Run {} for variable names {}.'.format(args.processing_step, var_names))


for var_name in var_names:
    print(var_name)

    # for ps, take ERA climate mean value (i.e. e.g. ps_historical.nc)
    if var_name == 'ps':
        var_file_name = era_climate_file_names.format(var_name)
    # for all other variables, take climate delta (i.e. e.g. tas_delta.nc)
    else:
        var_file_name = climate_delta_file_names.format(var_name)

    inp_file = os.path.join(args.input_directory, var_file_name)
    out_file = os.path.join(args.output_directory, var_file_name)


    ## smoothing
    if args.processing_step == 'smoothing':

        filter_data(inp_file, var_name, out_file)

    ## regridding
    elif args.processing_step == 'regridding':

        targ_lon = xr.open_dataset(args.era5_file_path).lon
        targ_lat = xr.open_dataset(args.era5_file_path).lat
        ds_in = xr.open_dataset(inp_file)
        ds_out = ds_in.interp(lon=targ_lon, lat=targ_lat)
        ds_out.to_netcdf(out_file)

