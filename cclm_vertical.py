# -*- coding: utf-8 -*-

import xarray as xr
import numpy as np
import sys
import time
from pathlib import Path
import os
from scipy.interpolate import RegularGridInterpolator
from functions import get_alt_half_level, get_alt_full_level

def vertinterpol(constant_path, datapath, var_name, 
                    outputpath, inputtimesteps, starttime=0):

    for stepnum in range(starttime,inputtimesteps):
        print('step {}'.format(stepnum))
        ds_in = xr.open_dataset(f"{datapath}/{var_name}{stepnum:05d}.nc")

        ds_const = xr.open_dataset(constant_path)
        vcoord = ds_const.vcoord
        hsurf = ds_const.HSURF

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










def vertinterpol_ORIG(terrainpath, datapath, variablename, outvar, 
    outputpath, vcflat, inputtimesteps, starttime=0):
    """
    Note: This function is specific to a certain verical grid in the cclm regional climate model and can only be used for that specific grid!
    If heights are inputed on file, make sure to copy it to the folder where the script will execute (in SCRATCH)
    """

    if os.path.exists('heights.txt'):
        #YUSPECIF gives height of half-levels (levels1; one more level than full)
        hlevels_flat_half=np.genfromtxt('heights.txt',skip_header=1)[:,1]
        plevels_cclm_half=np.genfromtxt('heights.txt',skip_header=1)[:,2]

        hlevels_flat=np.zeros(hlevels_flat_half.shape[0]-1)
        plevels_cclm=np.zeros(hlevels_flat_half.shape[0]-1)

        hlevels_flat = hlevels_flat_half[1:] + \
        (0.5 * (hlevels_flat_half[:-1] - hlevels_flat_half[1:]) )
        plevels_cclm = plevels_cclm_half[1:] + \
        (0.5 * (plevels_cclm_half[:-1] - plevels_cclm_half[1:]) )

        plevels_cclm = plevels_cclm * 100. #convert to pascal
    else:#use default levels
        print('Using default levels')
        plevels_cclm =np.asanyarray([13.2831,16.9293,21.1897,26.0707,31.5651,37.6459,44.2959,51.4832,59.0873,67.2557,75.7827,84.8931,94.3575,104.4696,114.9722,126.2605,138.0695,150.8587,164.3444,179.0299,194.4467,211.1740,228.6557,247.5546,267.2140,288.3890,310.3115,333.8365,358.0724,383.9844,411.6542,439.9525, 469.7307, 499.4267,529.9663,558.8535,587.7976,616.6264,645.1941,673.3691,701.0281,727.9977,754.1439,779.3544,803.5253,826.5168,848.2252,868.5767,887.5140,904.9537,920.8323,935.1420,947.8856,959.0633,968.6701,976.7777,983.8329,989.6878,994.2697,997.6310,1000.0000
])
        plevels_cclm = plevels_cclm * 100. #convert to pascal

        hlevels_flat = np.asanyarray([29798.5000  ,28255.9297 ,26824.0391 ,25497.3906 ,24269.3496 ,23134.0293 ,22081.8809 ,21105.6406 ,20207.7891 ,19360.5996  ,18576.5195 ,17827.8594  ,17127.9805 ,16451.0996 , 15811.5195  ,15183.5303  ,14581.3701  ,13982.0400 ,13400.0303 , 12815.3701 , 12248.1904  ,11678.5400 ,11126.4805  ,10572.0898  ,10035.4297  ,9496.5801 ,8975.5898  ,8452.5498  ,7947.5200  ,7440.5498  ,6931.9600 ,6442.3799 ,5956.4102  ,5498.0400  ,5050.8398  ,4647.9600  ,4261.9102  ,3893.2600   ,3542.1499  ,3208.5200 , 2892.2300  ,2593.7100 , 2312.9500   ,2049.7500  ,1803.8900  ,1575.5699   ,1364.6801  ,1170.9000   ,993.8400  ,833.4400  ,689.5300  ,561.5200 , 448.8200  , 350.9500 , 267.5500 ,197.6700   ,137.2300  ,87.3300  ,48.4400 ,20.0000 ,0.0000 ])


    outlevels = np.arange(len(hlevels_flat),dtype='int32')

    terrain = xr.open_dataset(terrainpath)['HSURF'].squeeze()
    terrain.values[terrain.values < 0] = 0

    smoothing = (vcflat - hlevels_flat) / vcflat
    smoothing = np.where(smoothing > 0, smoothing, 0)

    Path(outputpath).mkdir(parents=True, exist_ok=True)

    for stepnum in range(starttime,inputtimesteps):
        print('step {}'.format(stepnum))
        data = xr.open_dataset(f"{datapath}/{variablename}{stepnum:05d}.nc")[variablename]

        data = data.interp(plev=plevels_cclm)
        data['plev'] = hlevels_flat
        data = data.rename({'plev':'level'})

        newdata = np.zeros((len(hlevels_flat), terrain.shape[0], terrain.shape[1]))
        #compute the actual height of all model grid points
        newheights = hlevels_flat[:,None,None] + terrain.values[None,:,:] * smoothing[:,None,None]

        #interpolater needs to be ascending, therefore multiply by -1
        neg_hlev_flat = hlevels_flat * -1
        neg_newheigths = newheights * -1

        #grid dimensions for interpolation function
        xx = np.arange(terrain.shape[0])
        yy = np.arange(terrain.shape[1])

        #get index for efficient computation (avoid loops)
        xid, yid = np.ix_(xx,yy)

        #get the 3D interpolation fucntion
        fn = RegularGridInterpolator((neg_hlev_flat, xx, yy), data.values)

        #interpolate the data to the actual height in the model
        data.values = fn((neg_newheigths, xid, yid))

#		data.level.data.assign(outlevels)
        data.assign_coords(level=outlevels)#small bug fixed

        data = data.to_dataset(name=outvar)
        #try to fix weird saving problem; shold not be necessary if not more one job per node is requested.
        data.to_netcdf(f"{outputpath}/{outvar}{stepnum:05d}.nc")
        #dummy=0
        #while dummy < 1000:
        #    try:
        #            data.to_netcdf(f"{outputpath}/{outvar}{stepnum:05d}.nc")
        #            data.close()
        #            dummy = 2000
        #    except:
        #            time.sleep(0.5)
        #            dummy = dummy+1


if __name__ == "__main__":
	constant_path=str(sys.argv[1])
	datapath=str(sys.argv[2])
	var_name=str(sys.argv[3])
	#outvar=str(sys.argv[4])
	outputpath=str(sys.argv[4])
	#vcflat=int(sys.argv[5])
	inputtimesteps=int(sys.argv[5])
	starttime=int(sys.argv[6])


	vertinterpol(constant_path, datapath, var_name, outputpath, inputtimesteps, starttime=starttime)
	#vertinterpol_ORIG(terrainpath, datapath, variablename, outvar, outputpath, vcflat, inputtimesteps, starttime=starttime)
