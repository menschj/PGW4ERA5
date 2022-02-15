#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
description     Preprocess model data while loading. Do model specific
                stuff.
author			Christoph Heim
date created    09.07.2019
date changed    15.12.2021
usage           no args
"""
###############################################################################
import os#, copy, warnings, random
import numpy as np
import xarray as xr
from scipy.interpolate import interp1d
#from datetime import datetime, timedelta
#from numba import jit, njit
#from package.nl_models import nlm, models_cmip6
#from package.nl_variables import (nlv, add_var_attributes,
#                                  dimx,dimy,dimz,dimt)
from package.utilities import (dt64_to_dt, subsel_domain, 
                                select_common_timesteps, Timer)
#from base.nl_global import inp_glob_base_dir, model_specifics_path
#from package.constants import CON_G, CON_RD
###############################################################################

#MODEL_PP        = 'model_pp'
#MODEL_PP_DONE   = 'done'
#
#inp_base_dir = inp_glob_base_dir
#
#debug_level_2 = 2
#debug_level_4 = 4



def interp_plev_to_alt_new(x_out, x_in, data_in):
    # TODO: ACCESS models have nan for data within orography --> whole
    #       column becomes nan over orography.
    # TODO: implement adjustment to orography.
    #out = np.interp(x_out, x_in, data_in)
    #print(out)
    ##return(out)
    #quit()
    #f = interp1d(x_in, data_in, kind='linear', 
    #            fill_value=np.nan, bounds_error=False)
    #f = interp1d(x_in, data_in, kind='linear', 
    #            fill_value='extrapolate', bounds_error=False)
    x_in2 = x_in.copy()
    #if np.sum(np.isnan(data_in)) > 0:
    #    print(~np.isnan(data_in))
    #    quit()
    #    #x_in2 = x_in.copy()
    #    #x_in[np.isnan(x_in)] = 0
    #    #x_in = np.where(~np.isnan(x_in), x_in, 0)
    #    x_in2[0] = 1
    #    x_in2[1] = 2
    #    x_in2[2] = 3
    #    x_in2[3] = 4
    #    #print(x_in2)
    #    #quit()
    #print(x_in2)
    #print(x_out)
    data_in2 = data_in.copy()
    data_in2[np.isnan(data_in2)] = 0
    #print(data_in)
    f = interp1d(x_in2, data_in2, kind='cubic', 
                fill_value='extrapolate', bounds_error=False)
    out = f(x_out)
    #print(out)
    #if np.sum(np.isnan(data_in)) > 0:
    #    quit()
    #print(out)
    #quit()
    return(out)




def compute_P_plev(inp_file_path, out_file_path):

    ds = xr.open_dataset(inp_file_path)

    vert_key = 'plev'
    lat_key = 'lat'
    lon_key = 'lon'
    fact = 1
    # take pressure from vertical coordinate and expand to 3D field
    ds = ds.rename({'ta':'pa'})
    ds['pa'].values = (fact * ds[vert_key].expand_dims(dim={
                            'time':ds.time,
                            lat_key:ds[lat_key],
                            lon_key:ds[lon_key]}).transpose(
                                'time', vert_key, lat_key, lon_key)
                    ).values
    ds.to_netcdf(out_file_path)



def org_interp(var_name, inp_var_file, inp_alt_file, 
                out_var_file):
    """
    """

    ds_plev = xr.open_dataset(inp_var_file)
    alt = xr.open_dataset(inp_alt_file).zg

    #targ_alt = np.sort(np.loadtxt(os.path.join(model_specifics_path,
    #            'MPI-ESM1-2-HR_plev_mean_alt_historical_1985_2015.dat')))
    targ_alt = np.sort(np.loadtxt('target_alt.dat'))
    #targ_alt[0] = 500
    #print(targ_alt)
    #quit()
    #alt, var = select_common_timesteps(alt, var)

    ds_alt = ds_plev.copy()


    ## interpolate plev to alt
    ######################################################################
    ###### NEW WAY 
    #print(var)
    #quit()
    var_alt = xr.apply_ufunc(
        interp_plev_to_alt_new,
        targ_alt,
        alt,
        ds_plev[var_name],
        input_core_dims=[["plev"], ["plev"], ["plev"]],  # list with one entry per arg
        output_core_dims=[["alt"]],  # returned data has one dimension
        exclude_dims=set(("plev",)),  # dimensions allowed to change size. Must be a set!
        vectorize=True,  # loop over non-core dims
    ).transpose('time', 'alt', 'lat', 'lon')
    #print(var_alt)
    #print()

    ds_alt[var_name].values = var_alt.values
    ds_alt = ds_alt.rename({'plev':'alt'})
    ds_alt = ds_alt.assign_coords({'alt':targ_alt})
    #ds_alt.to_netcdf('test.nc')
    ds_alt.to_netcdf(out_var_file)
    #quit()

    #print(var.lat.values)
    #print(dep_vars['HSURF'].lat.values)
    #print(dep_vars['HSURF'].lat.values - var.lat.values)
    #quit()
    #try:
    #    var = var.where(var.alt > dep_vars['HSURF'], np.nan)
    #except ValueError:
    #    print('Warning: Coords of 3D field and HSURF differ. model_pp {}.'.format(mkey))
    #    dep_vars['HSURF'] = dep_vars['HSURF'].assign_coords({'lon':var.lon.values})
    #    dep_vars['HSURF'] = dep_vars['HSURF'].assign_coords({'lat':var.lat.values})
    #    var = var.where(var.alt > dep_vars['HSURF'], np.nan)





###############################################################################
###############################################################################
###############################################################################


if __name__ == '__main__':
    
    model_path = os.path.join('/net/o3/hymet_nobackup/heimc/data/pgw',
                            'MPI-ESM1-2-HR')
    experiments = ['historical', 'ssp585']

    var_names = ['ta', 'ua', 'va', 'hur', 'pa']
    for var_name in var_names:
        #var_name = 'ta'
        print(var_name)

        for experiment in experiments:
            #experiment = 'historical'
            print(experiment)

            inp_var_file = os.path.join(model_path,'plev_{}_{}.nc'.format(var_name, experiment)) 
            inp_alt_file = os.path.join(model_path,'plev_zg_{}.nc'.format(experiment)) 
            out_var_file = os.path.join(model_path,'alt_{}_{}.nc'.format(var_name, experiment)) 
            org_interp(var_name, inp_var_file, inp_alt_file, out_var_file)




    #for experiment in experiments:
    #    #experiment = 'historical'
    #    inp_file_path = os.path.join(model_path,'plev_ta_{}.nc'.format(experiment)) 
    #    out_file_path = os.path.join(model_path,'plev_pa_{}.nc'.format(experiment)) 
    #    compute_P_plev(inp_file_path, out_file_path)
