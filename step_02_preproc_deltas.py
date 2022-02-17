#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
description     PGW for ERA5 preprocessing of climate deltas
authors		    Before 2022: original developments by Roman Brogli
                Since 2022:  upgrade to PGW for ERA5 by Christoph Heim 
"""
##############################################################################
from settings import i_use_xesmf_regridding
import os, argparse
import xarray as xr
if i_use_xesmf_regridding:
    import xesmf as xe
import numpy as np
from pathlib import Path
from functions import filter_data
from settings import (
    i_debug,
    era_climate_file_name_base,
    climate_delta_file_name_base,
    LON_ERA, LAT_ERA,
    TIME_GCM, PLEV_GCM, LON_GCM, LAT_GCM,
)
##############################################################################

## input arguments
parser = argparse.ArgumentParser(description =
            'PGW for ERA5: Preprocess climate deltas before adding ' +
            'them to ERA5. The main routine (pgw_for_era5.py) requires ' +
            'climate deltas for [ta,hur,ua,va,zg,hurs,tas], as well ' +
            'as the ERA climatological value of ps. Note that that some ' +
            'additional settings can be made in settings.py.')

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
parser.add_argument('-i', '--input_dir', type=str,
            help='Directory with input files for selected ' +
            'processing step.')

# output directory
parser.add_argument('-o', '--output_dir', type=str,
            help='Directory where output files of selected ' +
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

        targ_lon = xr.open_dataset(args.era5_file_path)[LON_ERA]
        targ_lat = xr.open_dataset(args.era5_file_path)[LAT_ERA]
        ds_gcm = xr.open_dataset(inp_file)

        ## determine if GCM data set is periodic
        dlon_gcm = np.median(np.diff(ds_gcm[LON_GCM].values))
        dlat_gcm = np.median(np.diff(ds_gcm[LAT_GCM].values))
        if (dlon_gcm + np.max(ds_gcm[LON_GCM].values) - 
                      np.min(ds_gcm[LON_GCM].values)) >= 359.9:
            periodic_lon = True
            if i_debug >= 1:
                print('Regridding: Use periodic boundary conditions as GCM ' +
                       'data appears to be periodic in longitudinal ' +
                       'direction.')
        else:
            periodic_lon = False


        ## Interpolate with XESMF
        ## XESMF alters the exact values of the latitude coordinate a little
        ## bit which was found to be problematic. Therefore, there is an
        ## xarray-only implmenetation below.
        if i_use_xesmf_regridding:
            ds_targ = xr.open_dataset(args.era5_file_path)
            ds_in = ds_gcm
            regridder = xe.Regridder(ds_in, ds_targ, "bilinear", 
                                     periodic=periodic_lon)
            print(regridder)
            ds_out = regridder(ds_in[var_name])
            ds_gcm = ds_out.to_dataset(name=var_name)
            # keep attributes of variables and coordinates
            for field in [var_name, TIME_GCM, PLEV_GCM, LAT_GCM, 
                          LON_GCM, 'height']:
                if field in ds_in:
                    ds_gcm[field].attrs = ds_in[field].attrs
            # keep global attributes
            ds_gcm.attrs = ds_in.attrs

        ## Interpolate without XESMF. The results should be identical to XESMF
        ## except for tiny differences that appear to originate from
        ## numerical precision.
        else:
            #### LATITUDE INTERPOLATION
            ######################################################################
            ## make sure latitude is increasing with index
            if (ds_gcm[LAT_GCM].isel({LAT_GCM:0}).values > 
                ds_gcm[LAT_GCM].isel({LAT_GCM:-1}).values):
                if i_debug >= 1:
                    print('Regridding: GCM data has opposite ' +
                          'order of latitude. Apply reindexing.')
                # flip latitude dimension
                ds_gcm = ds_gcm.reindex(
                        {LAT_GCM:list(reversed(ds_gcm[LAT_GCM]))})

            ## in case ERA5 northmost/southmost grid points are just slightly
            ## closer to the pole than in the GCM dataset, add one grid point
            ## North and South in the GCM dataset with boundary values.
            north = ds_gcm.isel({LAT_GCM:-1})
            north[LAT_GCM].values += dlat_gcm
            #TODO debug
            #print(north[var_name].mean(dim=[LON_GCM]))
            #print(north[var_name])
            north[LAT_GCM].values = 90
            north[var_name] = north[var_name].mean(dim=[LON_GCM])

            south = ds_gcm.isel({LAT_GCM:0})
            south[LAT_GCM].values -= dlat_gcm

            south[LAT_GCM].values = -90
            south[var_name] = south[var_name].mean(dim=[LON_GCM])

            ds_gcm = xr.concat([south,ds_gcm,north], dim=LAT_GCM)

            ## make sure there is no extrapolation to the North and South
            if ( (np.max(targ_lat.values) > np.max(ds_gcm[LAT_GCM].values)) |
                 (np.min(targ_lat.values) < np.min(ds_gcm[LAT_GCM].values))):
                print('GCM lat: min {} max {}'.format(
                                np.min(ds_gcm[LAT_GCM].values),
                                np.max(ds_gcm[LAT_GCM].values))) 
                print('ERA5 lat: min {} max {}'.format(
                                np.min(targ_lat.values),
                                np.max(targ_lat.values))) 
                raise ValueError('ERA5 dataset extends further North or South ' +
                                  'than GCM dataset!. Perhaps consider using ' +
                                  'ERA5 on a subdomain only if global coverage ' +
                                  'is not required?') 

            ## run interpolation
            ds_gcm = ds_gcm.interp({LAT_GCM:targ_lat})

            #### LONGITUDE INTERPOLATION
            ######################################################################
            ## in case ERA5 westmost/eastmost grid points are just slightly
            ## furhter west/east than in the GCM dataset, add one grid point
            ## West and East in the GCM dataset with i) boundary values if data
            ## set ist not periodic, or ii) periodic values if data set is
            ## periodic.
            if periodic_lon:
                west = ds_gcm.isel({LON_GCM:0})
                east = ds_gcm.isel({LON_GCM:-1})
            else:
                west = ds_gcm.isel({LON_GCM:-1})
                east = ds_gcm.isel({LON_GCM:0})
            west[LON_GCM].values = (ds_gcm[LON_GCM].isel({LON_GCM:-1}).values + 
                                    dlon_gcm)
            east[LON_GCM].values = (ds_gcm[LON_GCM].isel({LON_GCM:0}).values -
                                    dlon_gcm)
            ds_gcm = xr.concat([east,ds_gcm,west], dim=LON_GCM)

            ## make sure there is no extrapolation to the East and West
            if ( (np.max(targ_lon.values) > np.max(ds_gcm[LON_GCM].values)) |
                 (np.min(targ_lon.values) < np.min(ds_gcm[LON_GCM].values))):
                print('GCM lon: min {} max {}'.format(
                                np.min(ds_gcm[LON_GCM].values),
                                np.max(ds_gcm[LON_GCM].values))) 
                print('ERA5 lon: min {} max {}'.format(
                                np.min(targ_lon.values),
                                np.max(targ_lon.values))) 
                raise ValueError('ERA5 dataset extends further East or West ' +
                                  'than GCM dataset!. Perhaps consider using ' +
                                  'ERA5 on a subdomain only if global coverage ' +
                                  'is not required?') 

            ## run interpolation
            ds_gcm = ds_gcm.interp({LON_GCM:targ_lon})

        ## test for NaN
        if np.sum(np.isnan(ds_gcm[var_name])).values > 0:
            raise ValueError('NaN in GCM dataset after interpolation.')

        ds_gcm.to_netcdf(out_file)

