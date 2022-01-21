# -*- coding: utf-8 -*-

import xarray as xr
import numpy as np
import sys
import time
from pathlib import Path
import os
from scipy.interpolate import RegularGridInterpolator
from functions import get_alt_half_level, get_alt_full_level
from constants import CON_G

def vertinterpol(constant_path, datapath, var_name, 
                    outputpath, inputtimesteps, start_time=0):

    for stepnum in range(start_time,inputtimesteps):
        print('step {}'.format(stepnum))
        ds_in = xr.open_dataset(f"{datapath}/{var_name}{stepnum:05d}.nc")

        ds_const = xr.open_dataset(constant_path)
        print(ds_const)
        quit()

        cosmo_alt_hl = get_alt_half_level(vcoord, hsurf)
        cosmo_alt = get_alt_full_level(cosmo_alt_hl)

        #data = data.interp(alt=cosmo_alt)
        #print(data)
        #data.to_netcdf('test.nc')
        #quit()
        #data['plev'] = hlevels_flat
        #data = data.rename({'plev':'level'})

        #newdata = np.zeros((len(hlevels_flat), terrain.shape[0], terrain.shape[1]))
        ##compute the actual height of all model grid points
        #newheights = hlevels_flat[:,None,None] + terrain.values[None,:,:] * smoothing[:,None,None]

        ##interpolater needs to be ascending, therefore multiply by -1
        #neg_hlev_flat = hlevels_flat * -1
        #neg_newheigths = newheights * -1

        #grid dimensions for interpolation function
        xx = np.arange(len(hsurf.rlon))
        yy = np.arange(len(hsurf.rlat))

        #get index for efficient computation (avoid loops)
        yid, xid = np.ix_(yy,xx)
        #print(yid.shape)
        #print(xid.shape)
        #print(ds_in[var_name].values.shape)
        #print(cosmo_alt.values.shape)
        #quit()

        ## create output dataset
        ds_out = ds_in.copy()
        var_out = xr.DataArray(
            data = np.zeros((len(cosmo_alt.level),len(yy),len(xx))),
            dims = ['level','rlat','rlon'],
        )
        ds_out[var_name] = var_out

        #ds_out = ds_out.assign_coords({'level':cosmo_alt.level})
        ds_out = ds_out.assign_coords(
                {'level':vcoord.values[:-1] + np.diff(vcoord.values)/2})
                #dtype='int32')})
        ds_out = ds_out.drop('alt')
        #print(ds_in)
        #print()
        #print(ds_out)
        #quit()

        #get the 3D interpolation fucntion
        fn = RegularGridInterpolator((ds_in.alt.values, yy, xx),
                                    ds_in[var_name].isel(time=0).values)
                                    #bounds_error=False)

        #interpolate the data to the actual height in the model
        #data.values = fn((cosmo_alt.values, yid, xid))
        ds_out[var_name].values = fn((cosmo_alt.values, yid, xid))
        ds_out[var_name] = ds_out[var_name].expand_dims({'time':ds_in.time}, axis=0)
        #print(ds_out)

        #data.assign_coords(level=outlevels)#small bug fixed

        #data = data.to_dataset(name=outvar)
        #try to fix weird saving problem; shold not be necessary if not more one job per node is requested.
        ds_out.to_netcdf(f"{outputpath}/{var_name}{stepnum:05d}.nc")






if __name__ == "__main__":
	constant_path=str(sys.argv[1])
	datapath=str(sys.argv[2])
	var_name=str(sys.argv[3])
	#outvar=str(sys.argv[4])
	outputpath=str(sys.argv[4])
	#vcflat=int(sys.argv[5])
	inputtimesteps=int(sys.argv[5])
	start_time=int(sys.argv[6])


	vertinterpol(constant_path, datapath, var_name, 
                    outputpath, inputtimesteps, start_time=start_time)
