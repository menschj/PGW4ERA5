#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
description     Auxiliary functions for Postprocess_CLM
author		Christoph Heim based on developments by Roman Brogli
date created    12.01.2022
"""
##############################################################################
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from numba import njit, jit
from datetime import datetime,timedelta
from dateutil.relativedelta import relativedelta
from package.utilities import dt64_to_dt
from constants import CON_RD, CON_G
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator
##############################################################################

def hour_of_year(dt): 
    beginning_of_year = datetime(dt.year, 1, 1, tzinfo=dt.tzinfo)
    return(int((dt - beginning_of_year).total_seconds() // 3600))

def specific_to_relative_humidity(QV, P, T):
    """
    Compute relative humidity from specific humidity.
    """
    RH = 0.263 * P * QV *(np.exp(17.67*(T - 273.15)/(T-29.65)))**(-1)
    return(RH)

def relative_to_specific_humidity(RH, P, T):
    """
    Compute specific humidity from relative humidity.
    """
    QV = (RH  * np.exp(17.67 * (T - 273.15)/(T - 29.65))) / (0.263 * P)
    return(QV)


def fix_grid_coord_diffs(fix_ds, ref_ds):
    ### There may be small grid inconsistencies that have to be fixed... 
    ### the reason is that the cdo remapping produes them.
    if np.max(np.abs(fix_ds.rlon.values - ref_ds.rlon.values)) > 1E-5:
        raise ValueError('Grids differ!')
    if np.mean(np.abs(fix_ds.rlon.values - ref_ds.rlon.values)) > 1E-8:
        #print('fix small differences in inputs grids.')
        fix_ds['rlon'] = ref_ds.rlon
        fix_ds['rlat'] = ref_ds.rlat


def load_delta(delta_inp_path, var_name, date_time, 
                laf_time, diff_time_step,
                name_base='{}_delta.nc'):

    #delta_year = xr.open_dataset(os.path.join(delta_inp_path,
    #                        name_base.format(var_name)))
    ## replace delta year values with year of current date_time
    #for i in range(len(delta_year.time)):
    #    delta_year.time.values[i] = dt64_to_dt(
    #                delta_year.time[i]).replace(year=date_time.year)
    ## interpolate in time and select variable
    #delta = delta_year[var_name].interp(time=date_time, 
    #                            method='linear', 
    #                        ).expand_dims(dim='time', axis=0)
    ## make sure time is in the same format as in laf file
    #delta['time'] = laf_time


    name_base = name_base.split('.nc')[0] + '_{:05d}.nc'
    delta = xr.open_dataset(os.path.join(delta_inp_path,
            name_base.format(var_name, diff_time_step)))[var_name]
    # make sure time is in the same format as in laf file
    delta['time'] = laf_time
    
    return(delta)


def get_delta_era5(var, var_name, laffile, delta_inp_path,
                    diff_time_step, date_time):
    delta = load_delta(delta_inp_path, var_name, date_time, laffile.time,
                        diff_time_step)
    delta = delta.assign_coords({'lat':var.lat.values})
    return(delta)


def get_delta_interp_era5(var, target_P, var_name, laffile, delta_inp_path,
                        diff_time_step, date_time):
    delta = load_delta(delta_inp_path, var_name, date_time, laffile.time,
                        diff_time_step)
    delta = delta.assign_coords({'lat':var.lat.values})

    if var_name == 'ta':
        sfc_var_name = 'tas'
        delta_sfc = load_delta(delta_inp_path, sfc_var_name, 
                            date_time, laffile.time,
                            diff_time_step)
        delta_sfc = delta_sfc.assign_coords({'lat':var.lat.values})
        ps_hist = load_delta(delta_inp_path, 'ps', 
                            date_time, laffile.time,
                            diff_time_step,
                            name_base='{}_historical.nc')
        ps_hist = ps_hist.assign_coords({'lat':var.lat.values})
    else:
        delta_sfc = None
        ps_hist = None

    # interpolate delta onto ERA5 vertical grid
    delta = vert_interp_delta(delta, target_P, delta_sfc, ps_hist)
    return(delta)


def replace_delta_sfc(source_P, ps_hist, delta, delta_sfc):
    #ps_hist = 101300
    #print(source_P)
    #print(ps_hist)
    #print(delta)
    #print(delta_sfc)
    out_source_P = source_P.copy()
    out_delta = delta.copy()
    if ps_hist > np.max(source_P):
        sfc_ind = len(source_P) - 1
        out_source_P[sfc_ind] = ps_hist
        out_delta[sfc_ind] = delta_sfc
    elif ps_hist < np.min(source_P):
        raise ValueError()
    else:
        sfc_ind = np.max(np.argwhere(ps_hist > source_P))
        #print(sfc_ind)
        out_delta[sfc_ind:] = delta_sfc
        out_source_P[sfc_ind] = ps_hist

    #print(out_source_P)
    #print(out_delta)
    #quit()
    return(out_source_P, out_delta)


def vert_interp_delta(delta, target_P, delta_sfc, ps_hist):

    #print(delta)
    if delta.name in ['hus', 'QV']:
        print('WARNING: DEBUG MODE FOR variable hus model top!!!')
        top = xr.zeros_like(delta.sel(plev=100000))
        top['plev'] = 500
        delta = xr.concat([delta, top], dim='plev').transpose(
                                    'time', 'plev', 'lat', 'lon')
        #print(delta.mean(dim=['lat','lon','time']))
        #print(delta.plev)
        #quit()

    # sort delta dataset from top to bottom (pressure ascending)
    delta = delta.reindex(plev=list(reversed(delta.plev)))

    # create 4D source pressure with GCM pressure levels
    source_P = delta.plev.expand_dims(
                    dim={'lon':delta.lon,
                         'lat':delta.lat,
                         'time':delta.time}).transpose(
                                    'time', 'plev', 'lat', 'lon')

    #lon_ind = 22
    #lat_ind = 50
    #print(source_P.isel(lon=lon_ind, lat=lat_ind).values)
    #print(target_P.sel(level=np.max(target_P.level)).isel(lon=lon_ind, lat=lat_ind).values)
    #print(delta.isel(lon=lon_ind, lat=lat_ind).values)
    ##print(ps_hist.isel(lon=lon_ind, lat=lat_ind).values)
    #print(delta_sfc.isel(lon=lon_ind, lat=lat_ind).values)


    ## if surface values are given, replace them at the
    ## level of the surface pressure
    if delta_sfc is not None:
        source_P, delta = xr.apply_ufunc(
                replace_delta_sfc, source_P, 
                ps_hist, 
                #target_P.sel(level=np.max(target_P.level)), 
                delta, delta_sfc,
                input_core_dims=[['plev'],[],['plev'],[]],
                output_core_dims=[['plev'],['plev']],
                vectorize=True)
        source_P = source_P.transpose('time', 'plev', 'lat', 'lon')
        delta = delta.transpose('time', 'plev', 'lat', 'lon')


    #print(source_P.isel(lon=lon_ind, lat=lat_ind).values)
    #print(delta.isel(lon=lon_ind, lat=lat_ind).values)
    #quit()

    if source_P.dims != ('time', 'plev', 'lat', 'lon'):
        raise ValueError()
    if delta.dims != ('time', 'plev', 'lat', 'lon'):
        raise ValueError()
    if target_P.dims != ('time', 'level', 'lat', 'lon'):
        raise ValueError()

    delta_interp = interp_nonfixed(delta, source_P, target_P,
                        'plev', 'level',
                        extrapolate='constant')
    return(delta_interp)




def interp_nonfixed(var, source_P, targ_P,
                    inp_vdim_name, out_vdim_name,
                    extrapolate='off'):
    """
    extrapolate:
        - off: no extrapolation
        - linear: linear extrapolation
        - constant: constant extrapolation
    """
    if extrapolate not in ['off', 'linear', 'constant']:
        raise ValueError()

    targ = xr.zeros_like(targ_P)
    tmp = np.zeros_like(targ.values.squeeze())
    #interp_vprof(var.values.squeeze(), source_P.values.squeeze(),
    #                    targ_P.squeeze().values, 
    interp_vprof(var.values.squeeze(), np.log(source_P.values.squeeze()),
                        np.log(targ_P.squeeze()).values, 
                        tmp,
                        len(var.lat), len(var.lon),
                        extrapolate)
    tmp = np.expand_dims(tmp, axis=0)
    targ.values = tmp
    return(targ)



@njit()
def interp_extrap_1d(src_x, src_y, targ_x, extrapolate):
    targ_y = np.zeros(len(targ_x))
    for ti in range(len(targ_x)):
        i1 = -1
        i2 = -1
        require_extrap = False
        for si in range(len(src_x)):
            #print(src_x[si])
            ty = np.nan
            # extrapolate lower end
            if (si == 0) and src_x[si] > targ_x[ti]:
                if extrapolate == 'linear':
                    i1 = si
                    i2 = si + 1
                elif extrapolate == 'constant':
                    i1 = si
                    i2 = si
                require_extrap = True
                break
            # exact match
            elif src_x[si] == targ_x[ti]:
                i1 = si
                i2 = si
                break
            # upper src_x found (interpolation)
            elif src_x[si] > targ_x[ti]:
                i1 = si - 1
                i2 = si
                break
            # we are still smaller than targ_x[ti] 
            else:
                pass

        # extrapolate upper end
        if i1 == -1:
            if extrapolate == 'linear':
                i1 = len(src_x) - 2
                i2 = len(src_x) - 1
            elif extrapolate == 'constant':
                i1 = len(src_x) - 1 
                i2 = len(src_x) - 1
            require_extrap = True

        if require_extrap and extrapolate == 'off':
            raise ValueError('Extrapolation turned off but data out of bounds.')

        if i1 == i2:
            targ_y[ti] = src_y[i1]
        else:
            targ_y[ti] = (
                src_y[i1] + (targ_x[ti] - src_x[i1]) * 
                (src_y[i2] - src_y[i1]) / (src_x[i2] - src_x[i1])
            )

    return(targ_y)


@njit()
def interp_vprof(orig_array, src_p,
                targ_p, interp_array,
                nlat, nlon,
                extrapolate):
    """
    Helper function for compute_VARNORMI. Speedup of ~100 time
    compared to pure python code!
    """
    for lat_ind in range(nlat):
        for lon_ind in range(nlon):
            #print('lat {} lon {}'.format(lat_ind, lon_ind))
            src_val_col = orig_array[:, lat_ind, lon_ind]
            src_p_col = src_p[:, lat_ind, lon_ind]
            targ_p_col = targ_p[:, lat_ind, lon_ind]
            #print(targ_p_col.shape)
            ## np.interp does not work for extrapolation 
            #interp_col = np.interp(targ_p_col, src_p_col, src_val_col)
            ## scipty with extrapolation but slow
            #f = interp1d(src_p_col, src_val_col, fill_value='extrapolate')
            #interp_col = f(targ_p_col)
            ## faster implementation with numba

            #print(src_p_col)
            #print(src_val_col)
            #print(targ_p_col)
            interp_col = interp_extrap_1d(src_p_col, src_val_col, 
                                        targ_p_col, extrapolate)
            #print(interp_col)
            #quit()
            if np.any(np.isnan(interp_col)):
                #raise ValueError('Interpolated data contains NaN either due to '+
                #                'extrapolation turned off but data out of bounds, ' +
                #                'or because NaN are inherent to data!')
                raise ValueError('Interpolated data contains NaN!')
            interp_array[:, lat_ind, lon_ind] = interp_col
    #return(interp_array)


def interp(var, P, targ_p, inp_vdim_name):
    targ = xr.zeros_like(var).isel({inp_vdim_name:range(len(targ_p))})
    tmp = np.zeros_like(targ.values.squeeze())
    interp_vprof(var.values.squeeze(), np.log(P.values.squeeze()),
                        np.log(np.expand_dims(targ_p, axis=(1,2))), 
                        tmp,
                        len(var.lat), len(var.lon))
    tmp = np.expand_dims(tmp, axis=0)
    targ.values = tmp
    targ = targ.rename({inp_vdim_name:'plev'})
    targ = targ.assign_coords(plev=targ_p)
    return(targ)

def debug_interp(var, P=None):
    plevs = [102000, 100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000,
    20000, 15000, 10000]

    #lon = -15
    #lat = -15
    #var = var.isel({'time':0})
    #var = var.interp({'lon':lon,'lat':lat})
    #if P is not None:
    #    P = P.isel({'time':0})
    #    P = P.interp({'lon':lon,'lat':lat})
    #    var = var.assign_coords(level1=np.log(P.values))
    #    var = var.rename({'level1':'plev'})
    #    var = var.interp(plev=np.log(plevs), kwargs={'fill_value':'extrapolate'})
    #    var = var.assign_coords(plev=np.exp(var.plev))



    if P is not None:
        if 'level1' in var.dims:
            inp_vdim_name = 'level1'
        elif 'level' in var.dims:
            inp_vdim_name = 'level'
        else:
            raise NotImplementedError()
        var = interp(var, P, np.asarray(plevs), inp_vdim_name)
        var = var.mean(dim=['time','lon','lat'])
    else:
        var = var.mean(dim=['time','lon','lat'])

    return(var)


def determine_p_ref(p_era, p_pgw, p_ref_last, p_ref_opts):
    for p in p_ref_opts:
        if (p_era > p) & (p_pgw > p):
        #if x / p > 1.01:
            if p_ref_last is None:
                return(p)
            else:
                return(min(p, p_ref_last))

def integ_geopot_era5(P_hl, FIS, T, QV, level1, p_ref):
    """
    """
    # take log half-level pressure difference (located at full levels)
    dlnP = np.log(P_hl).diff(
                dim='level1', label='lower').rename({'level1':'level'})

    # create geopotential array and fill with surface geopotential
    PHI_hl = FIS.expand_dims(dim={'level1':level1}).copy()
    #PHI = laffile['FIS'].expand_dims(dim={'level':laffile.level}).copy()

    # compute virtual temperature
    TV = T * (1 + 0.61 * QV)

    ## integrate over model half levels
    for l in sorted(TV.level.values, reverse=True):
        #print(l)
        # geopotential at full level
        PHI_hl.loc[dict(level1=l)] = (
                PHI_hl.sel(level1=l+1) +
                (CON_RD * TV.sel(level=l) * dlnP.sel(level=l))
        )
        #print('{}   {}'.format(l, PHI_hl.loc[dict(level1=l)].mean().values/CON_G))

        ## geopotential at half level
        #alpha = 1. - (
        #        (laffile['P_hl'].sel(level1=l) / 
        #        (laffile['P_hl'].sel(level1=l+1) - laffile['P_hl'].sel(level1=l))) * 
        #        dlnP.sel(level=l)
        #)
        #PHI.loc[dict(level=l)] = (
        #        PHI_hl.sel(level1=l+1) +
        #        (CON_RD * TV.sel(level=l) * alpha)
        #)

            
    #PHI = PHI.transpose('time', 'level', 'lat', 'lon')
    PHI_hl = PHI_hl.transpose('time', 'level1', 'lat', 'lon')

    ## integrate from last half-level below reference pressure
    ## up to reference pressure
    p_diff = P_hl - p_ref
    p_diff = p_diff.where(p_diff >= 0, np.nan)
    ind_ref_star = p_diff.argmin(dim='level1')
    hl_ref_star = p_diff.level1.isel(level1=ind_ref_star)
    #hl_ref_star.to_netcdf('test.nc')
    #quit()
    p_ref_star = P_hl.sel(level1=hl_ref_star)
    phi_ref_star = PHI_hl.sel(level1=hl_ref_star)

    phi_ref = (
            phi_ref_star -
            (CON_RD * TV.sel(level=hl_ref_star-1)) * 
            (np.log(p_ref) - np.log(p_ref_star))
    )
    if 'level1' in phi_ref.coords:
        del phi_ref['level1']
    if 'level' in phi_ref.coords:
        del phi_ref['level']
    if 'plev' in phi_ref.coords:
        del phi_ref['level']
    #print(phi_ref)
    #quit()
    #del phi_ref['level1']
    #del phi_ref['level']

    return(PHI_hl, phi_ref)


def integ_pressure_upward_era5(laffile, phi_ref, hl_ind=None):
    """
    Function to compute pressure change from temperature and humidity
    changes to maintain hydrostatic balance.

    Hydrostatic equation to be solved:
    dln(p) / dz = -g / (R * Tv)
    """
    TV = laffile['T'] * (1 + 0.61 * laffile['QV'])

    lnP_hl = np.log(laffile['PS'].expand_dims(dim={'level1':laffile.level1}).copy())
    
    # compute altitude change accross full level (dz)
    dz = laffile['PHI_hl'].diff('level1', label='lower').rename(
            {'level1':'level'}) / CON_G

    # vertically integrate pressure
    for l in sorted(TV.level.values, reverse=True):
        #print(l)

        # change of log pressure with level
        dlnP = (
                - CON_G / (CON_RD * TV.sel(level=l)) * 
                dz.sel(level=l)
        )
        lnP_hl.loc[dict(level1=l)] = lnP_hl.sel(level1=l+1) - dlnP


    lnP_hl = lnP_hl.transpose('time', 'level1', 'lat', 'lon')
    

    ## integrate from last half-level below reference pressure
    ## up to reference pressure
    phi_diff = phi_ref - laffile['PHI_hl']
    phi_diff = phi_diff.where(phi_diff > 0, np.nan)
    #plt.plot(test.mean(dim=['lon','lat','time']),
    #        laffile['P_hl'].mean(dim=['lon','lat','time']))
    #plt.show()
    ind_ref_star = phi_diff.argmin(dim='level1')
    #ind_ref_star.to_netcdf('test.nc')
    #print(ind_ref_star)
    hl_ref_star = phi_diff.level1.isel(level1=ind_ref_star)
    #hl_ref_star.to_netcdf('test.nc')
    #print(hl_ref_star)
    phi_ref_star = laffile['PHI_hl'].sel(level1=hl_ref_star)
    #print((phi_ref-phi_ref_star).mean(dim=['lon','lat','time']).values)
    #p_ref_star.to_netcdf('test.nc')
    #print(p_ref_star)
    lnP_ref_star = lnP_hl.sel(level1=hl_ref_star)
    #phi_ref_star.to_netcdf('test.nc')
    #print(phi_ref_star)

    dlnP = (
            - CON_G / (CON_RD * TV.sel(level=hl_ref_star-1)) * 
            (phi_ref - phi_ref_star)/CON_G
    )
    lnP_ref = lnP_ref_star + dlnP
    del lnP_ref['level1']
    del lnP_ref['level']

    P_ref = np.exp(lnP_ref)
    P_ref_star = np.exp(lnP_ref_star)

    #print(P_ref.mean(dim=['lon','lat','time']).values)
    #print(P_ref_star.mean(dim=['lon','lat','time']).values)
    #quit()
    #print(P_ref)
    #quit()


    if hl_ind is not None:
        P_hl = np.exp(lnP_hl)
        P_ref = P_hl.sel(level1=hl_ind)


    return(P_ref)






def integ_geopot_downward_era5(laffile, p_ref, p_ref_star, phi_ref, hl_ref_star):
    """
    """
    # compute virtual temperature
    TV = laffile['T'] * (1 + 0.61 * laffile['QV'])

    ## integrate from reference pressure down to
    ## last half-level below reference pressure
    phi_ref_star = (
            phi_ref +
            (CON_RD * TV.sel(level=hl_ref_star-1)) * 
            (np.log(p_ref) - np.log(p_ref_star))
    )
    print(phi_ref_star)
    phi_ref_star.to_netcdf('test.nc')
    raise NotImplementedError()
    #quit()

    # take log half-level pressure difference (located at full levels)
    dlnP = np.log(laffile['P_hl']).diff(
            dim='level1', label='lower').rename({'level1':'level'})

    # create geopotential array and fill with surface geopotential
    PHI_hl = laffile['PHI_hl'].copy()


    for l in range(start_level,np.max(TV.level.values)+1):
        print(l)
        # geopotential at full level
        PHI_hl.loc[dict(level1=l+1)] = (
                PHI_hl.sel(level1=l) -
                (CON_RD * TV.sel(level=l) * dlnP.sel(level=l))
        )
        print(np.sum(PHI_hl.sel(level1=l+1) - laffile['FIS'] < 0).values)
        print((PHI_hl.sel(level1=l+1) - laffile['FIS']).mean(
            dim=['lon','lat','time']).values/CON_G)
    quit()

            
    PHI = PHI.transpose('time', 'level', 'lat', 'lon')
    PHI_hl = PHI_hl.transpose('time', 'level1', 'lat', 'lon')

    return(PHI, PHI_hl)




def integ_pressure_era5(laffile):
    """
    Function to compute pressure change from temperature and humidity
    changes to maintain hydrostatic balance.

    Hydrostatic equation to be solved:
    dln(p) / dz = -g / (R * Tv)
    """
    start_level = 30 
    #print(laffile['P_hl'].sel(level1=start_level))
    #quit()

    TV = laffile['T'] * (1 + 0.61 * laffile['QV'])

    lnP_hl = np.log(laffile['P_hl'])
    
    # compute altitude change accross full level (dz)
    dz = laffile['PHI_hl'].diff('level1', label='lower').rename(
            {'level1':'level'}) / CON_G

    # vertically integrate pressure
    for l in range(start_level,np.max(dz.level.values)+1):
        #print(l)

        # change of log pressure with level
        dlnP = (
                - CON_G / (CON_RD * TV.sel(level=l)) * 
                dz.sel(level=l)
        )
        lnP_hl.loc[dict(level1=l+1)] = lnP_hl.sel(level1=l) + dlnP

    
    P_hl = np.exp(lnP_hl)

    #print(P_hl.mean(dim=['time','lon','lat']))
    #print(laffile['P_hl'].mean(dim=['time','lon','lat']))
    #print((P_hl-laffile['P_hl']).mean(dim=['time','lon','lat']))
    
    #plt.plot((P_hl-laffile['P_hl']).mean(dim=['time','lon','lat']))
    #plt.show()

    # take surface pressure
    PS = P_hl.sel(level1=np.max(P_hl.level1))
    return(PS)
















def add_delta(var, laffile, delta_inp_path, diff_time_step):
    delta = xr.open_dataset(
                f'{delta_inp_path}/{var}{diff_time_step:05d}.nc')[var]
    ## HCH2021 start
    ## There are small grid inconsistencies that have to be fixed... 
    fix_grid_coord_diffs(delta, laffile)
    fix_grid_coord_diffs(delta, laffile)
    ## HCH 2021 stop

    if 'time' in laffile[var].dims:
        laffile[var].data[0,::] = (laffile[var].data[0,::] + 
                            delta.data.astype('float32'))
    else:
        laffile[var].data = (laffile[var].data + 
                            delta.data.astype('float32'))
    #quit()



def get_alt_half_level(vcoord, hsurf):
    """
    Compute altitude of half levels (Gal-Chen coordinate)
    """
    height_flat_half_level = vcoord #these are half levels

    smoothing = (vcoord.vcflat - height_flat_half_level) / vcoord.vcflat
    smoothing = smoothing.where(smoothing > 0, 0)

    hsurf = hsurf.squeeze()

    #alt_half_level = (
    #    height_flat_half_level.values[:,None,None] + 
    #    hsurf.values[None,:,:] * smoothing.values[:,None,None]
    #)
    #print(alt_half_level.mean((1,2)))
    alt_half_level = (
            height_flat_half_level + 
            (hsurf * smoothing).sortby('level1', ascending=True)
    )
    return(alt_half_level)



def get_alt_full_level(alt_half_level):
    """
    Approximate full level based on half levels
    """
    # approximate full level by 1/2 half level
    alt_full_level = (
            alt_half_level.isel(level1=range(0,len(alt_half_level.level1)-1)) +
            alt_half_level.diff('level1') / 2
    )
    alt_full_level = alt_full_level.rename({'level1':'level'})
    return(alt_full_level)



def get_pref(vcoord, hsurf, hl_or_fl):
    """
    Compute reference pressure for COSMO.
    hl_or_fl:   if ='hl', reference pressure at half levels is computed.
                if ='fl', reference pressure at full levels is computed.
    """
    # altitude of model half levels
    if hl_or_fl == 'hl':
        alt = get_alt_half_level(vcoord, hsurf)
    elif hl_or_fl == 'fl':
        alt = get_alt_half_level(vcoord, hsurf)
        alt = get_alt_full_level(alt)
    else:
        raise ValueError()

    ##### Code taken from cosmo source code: vgrid_refatm_utils.f90
    ##### subroutine reference_atmosphere_2 and is
    ##### only valid for irefatm=2 and ivctype=2

    # Constants
    p0sl = vcoord.p0sl # sea-level pressure [Pa]
    t0sl = vcoord.t0sl   # sea-level temperature [K]
    dt0lp = vcoord.dt0lp   # Beta value (dT0/dln(p)) [K/Pa]
    # irefatm = 2
    delta_t = vcoord.delta_t
    h_scal = vcoord.h_scal
    # Source: COSMO description Part VII, page 12
    t00 = t0sl - delta_t

    pref = p0sl * np.exp (- CON_G / CON_RD * h_scal / t00 * \
               np.log((np.exp(alt / h_scal) * t00 + delta_t) / \
                      (t00 + delta_t)) )

    return(pref)


def get_pref_sfc(vcoord, hsurf):
    """
    Compute reference pressure for COSMO at surface.
    hl_or_fl:   if ='hl', reference pressure at half levels is computed.
                if ='fl', reference pressure at full levels is computed.
    """
    hsurf = hsurf.isel(time=0)

    ##### Code taken from cosmo source code: vgrid_refatm_utils.f90
    ##### subroutine reference_atmosphere_2 and is
    ##### only valid for irefatm=2 and ivctype=2

    # Constants
    p0sl = vcoord.p0sl # sea-level pressure [Pa]
    t0sl = vcoord.t0sl   # sea-level temperature [K]
    dt0lp = vcoord.dt0lp   # Beta value (dT0/dln(p)) [K/Pa]
    # irefatm = 2
    delta_t = vcoord.delta_t
    h_scal = vcoord.h_scal
    # Source: COSMO description Part VII, page 12
    t00 = t0sl - delta_t

    pref_sfc = p0sl * np.exp (- CON_G / CON_RD * h_scal / t00 * \
               np.log((np.exp(hsurf / h_scal) * t00 + delta_t) / \
                      (t00 + delta_t)) )

    return(pref_sfc)












def interp_cubic(x_out, x_in, data_in):
    f = interp1d(x_in, data_in, kind='cubic')
                #fill_value='extrapolate', bounds_error=False)
    out = f(x_out)
    return(out)





























def adjust_pressure_to_new_climate2(delta_inp_dir, delta_time_step, 
                    BC_file, alt_hl, alt, pref_hl, pref):
    """
    Function to compute pressure change from temperature and humidity
    changes to maintain hydrostatic balance.

    Hydrostatic equation to be solved:
    dln(p) / dz = -g / (R * Tv)
    """

    ## Load variables in base climate state
    nlevels = len(BC_file['PP'].level)
    P0  = BC_file['PP'].isel(time=0).assign_coords(level=range(nlevels)) + pref
    T0  = BC_file['T' ].isel(time=0).assign_coords(level=range(nlevels))
    QV0 = BC_file['QV'].isel(time=0).assign_coords(level=range(nlevels))
    #QC0 = BC_file['QC'].isel(time=0).assign_coords(level=range(nlevels))
    #QI0 = BC_file['QI'].isel(time=0).assign_coords(level=range(nlevels))
    RH0 = specific_to_relative_humidity(QV0, P0, T0)

    #plt.plot(np.log(P0.mean(dim=['rlon','rlat'])))
    #plt.show()

    ## Load change of variables between climate states
    delta_T = xr.open_dataset(os.path.join(delta_inp_dir,
                    'T{:05d}.nc'.format(delta_time_step)))['T']
    delta_T = delta_T.assign_coords(level=range(nlevels))
    if 'time' in delta_T.dims:
        delta_T = delta_T.isel(time=0)
    delta_RH = xr.open_dataset(os.path.join(delta_inp_dir,
                    'RELHUM{:05d}.nc'.format(delta_time_step)))['RELHUM']
    delta_RH = delta_RH.assign_coords(level=range(nlevels))
    if 'time' in delta_RH.dims:
        delta_RH = delta_RH.isel(time=0)
    delta_P = xr.open_dataset(os.path.join(delta_inp_dir,
                    'PP{:05d}.nc'.format(delta_time_step)))['PP']
    delta_P = delta_P.assign_coords(level=range(nlevels))
    if 'time' in delta_P.dims:
        delta_P = delta_P.isel(time=0)
    #print(delta_P.sel(level=0).mean(dim=['rlon','rlat']).values)
    #quit()


    ## There may be small grid inconsistencies that have to be fixed... 
    fix_grid_coord_diffs(delta_T, BC_file)
    fix_grid_coord_diffs(delta_RH, BC_file)
    fix_grid_coord_diffs(delta_P, BC_file)

    ## Compute variables in scenario climate state
    T1 = T0 + delta_T
    QV1 = relative_to_specific_humidity(RH0+delta_RH, P0, T1)

    ## Interpolate full level fields to half levels
    P0_hl = (P0.sel(level=range(1,len(P0.level))) -
            P0.diff('level') / 2).rename({'level':'level1'})
    T0_hl = (T0.sel(level=range(1,len(T0.level))) -
            T0.diff('level') / 2).rename({'level':'level1'})
    QV0_hl = (QV0.sel(level=range(1,len(QV0.level))) -
            QV0.diff('level') / 2).rename({'level':'level1'})
    #QC0_hl = (QC0.sel(level=range(1,len(QC0.level))) -
    #        QC0.diff('level') / 2).rename({'level':'level1'})
    #QI0_hl = (QI0.sel(level=range(1,len(QI0.level))) -
    #        QI0.diff('level') / 2).rename({'level':'level1'})
    T1_hl = (T1.sel(level=range(1,len(T1.level))) -
            T1.diff('level') / 2).rename({'level':'level1'})
    QV1_hl = (QV1.sel(level=range(1,len(QV1.level))) -
            QV1.diff('level') / 2).rename({'level':'level1'})

    TV0_hl = T0_hl * (1 + 0.61 * QV0_hl)
    TV1_hl = T1_hl * (1 + 0.61 * QV1_hl)
    # using QC and QI does not reduce the error (also not increase)
    #TV0_hl = T0_hl * (1 + 0.61 * QV0_hl - QC0_hl - QI0_hl)

    #P0_comp_dir = xr.zeros_like(P0)
    lnP0_comp_ln = xr.zeros_like(P0)
    #P1_comp_dir = xr.zeros_like(P0)
    #lnP1_comp_ln = xr.zeros_like(P0)
    lnP1_comp_ln_2 = xr.zeros_like(P0)

    # take upper-most value from BC
    top_pressure_0 = P0.sel(level=0)
    top_pressure_1 = P0.sel(level=0)
    #print(top_pressure_1.mean(dim=['rlon','rlat']).values)
    top_pressure_1_2 = P0.sel(level=0) + 0.86*delta_P.sel(level=0)
    #top_pressure_1_2 = P0.sel(level=0) + delta_P.sel(level=0)

    #P0_comp_dir.loc[dict(level=0)] = top_pressure_0
    lnP0_comp_ln.loc[dict(level=0)] = np.log(top_pressure_0)
    #P1_comp_dir.loc[dict(level=0)] = top_pressure_1 
    #lnP1_comp_ln.loc[dict(level=0)] = np.log(top_pressure_1)
    lnP1_comp_ln_2.loc[dict(level=0)] = np.log(top_pressure_1_2)
    
    # compute altitude change accross full level (dz)
    dz = alt_hl.diff('level1').rename({'level1':'level'})
    dz_hl = alt.diff('level').rename({'level':'level1'})
    dz_hl['level1'] = dz_hl['level1'] + 1

    #print(dz.mean(dim=['rlon','rlat']))
    #print(dz.level)
    #print(dz_hl.mean(dim=['rlon','rlat']))
    #print(dz_hl.level1)

    ## Method 1
    for lvi in range(0,np.max(dz.level.values)):
        #print(lvi)
        ## vertically integrate pressure

        ## base state
        dlnP = (
                - CON_G / (CON_RD * TV0_hl.sel(level1=lvi+1)) * 
                dz_hl.sel(level1=lvi+1)
        )
        lnP0_comp_ln.loc[dict(level=lvi+1)] = lnP0_comp_ln.sel(level=lvi) + dlnP

        ### future state
        #dlnP = (
        #        - CON_G / (CON_RD * TV1_hl.sel(level1=lvi+1)) * 
        #        dz_hl.sel(level1=lvi+1)
        #)
        #lnP1_comp_ln.loc[dict(level=lvi+1)] = lnP1_comp_ln.sel(level=lvi) + dlnP

        ## future state (test)
        dlnP = (
                - CON_G / (CON_RD * TV1_hl.sel(level1=lvi+1)) * 
                dz_hl.sel(level1=lvi+1)
        )
        lnP1_comp_ln_2.loc[dict(level=lvi+1)] = lnP1_comp_ln_2.sel(level=lvi) + dlnP

    ## convert back to real pressure
    P0_comp_ln = np.exp(lnP0_comp_ln)
    #P1_comp_ln = np.exp(lnP1_comp_ln)
    P1_comp_ln_2 = np.exp(lnP1_comp_ln_2)


    ## Method 2
    #P1 = P0 + P1_comp_ln-P0_comp_ln
    #P1_hl = (P1.sel(level=range(1,len(P1.level))) -
    #        P1.diff('level') / 2).rename({'level':'level1'})
    #for lvi in range(0,np.max(dz.level.values)):
    #    print(lvi)
    #    ## vertically integrate pressure

    #    ## base state
    #    rho0_hl = P0_hl.sel(level1=lvi+1) / (CON_RD * TV0_hl.sel(level1=lvi+1))
    #    dP = - rho0_hl * CON_G * dz_hl.sel(level1=lvi+1)
    #    P0_comp_dir.loc[dict(level=lvi+1)] = P0_comp_dir.sel(level=lvi) + dP

    #    ## future state
    #    #rho1_hl = P0_hl.sel(level1=lvi+1) / (CON_RD * TV1_hl.sel(level1=lvi+1))
    #    rho1_hl = P1_hl.sel(level1=lvi+1) / (CON_RD * TV1_hl.sel(level1=lvi+1))
    #    dP = - rho1_hl * CON_G * dz_hl.sel(level1=lvi+1)
    #    P1_comp_dir.loc[dict(level=lvi+1)] = P1_comp_dir.sel(level=lvi) + dP


    ## compute pressure change
    #delta_P_dir = P1_comp_dir-P0_comp_dir
    #delta_P_ln = P1_comp_ln-P0_comp_ln
    delta_P_ln_2 = P1_comp_ln_2-P0_comp_ln


    ## plot results (error)
    #handles = []
    ##handle, = plt.plot((P0_comp_dir-P0).mean(dim=['rlon', 'rlat']),
    ##        alt.mean(dim=['rlon', 'rlat']), label='direct')
    ##handles.append(handle)
    #handle, = plt.plot((P0_comp_ln-P0).mean(dim=['rlon', 'rlat']),
    #        alt.mean(dim=['rlon', 'rlat']), label='ln')
    #plt.ylabel('altitude [m]')
    #plt.xlabel('error [Pa]')
    #handles.append(handle)
    #plt.legend(handles=handles)
    #plt.show()
    #quit()


    ## plot results (change)
    #handles = []
    ##handle, = plt.plot(delta_P_dir.mean(dim=['rlon', 'rlat']),
    ##        alt.mean(dim=['rlon', 'rlat']), label='direct')
    ##handles.append(handle)
    #handle, = plt.plot(delta_P_ln_2.mean(dim=['rlon', 'rlat']),
    #        alt.mean(dim=['rlon', 'rlat']), label='ln p_top(GCM)')
    #handles.append(handle)
    #handle, = plt.plot(delta_P_ln.mean(dim=['rlon', 'rlat']),
    #        alt.mean(dim=['rlon', 'rlat']), label='ln')
    #handles.append(handle)
    #handle, = plt.plot(delta_P.mean(dim=['rlon', 'rlat']),
    #        alt.mean(dim=['rlon', 'rlat']), label='GCM')
    #handles.append(handle)
    #plt.ylabel('altitude [m]')
    #plt.xlabel('change [Pa]')
    #plt.legend(handles=handles)
    #plt.show()
    #quit()

    # replace PP in BC_file
    BC_file['PP'].data = BC_file['PP'] + delta_P_ln_2
    #print(BC_file['PP'])

    return(BC_file)








def adjust_pressure_to_new_climate(delta_inp_dir, delta_time_step, 
                    BC_file, alt_hl, alt, pref_hl, pref):
    """
    Function to compute pressure change from temperature and humidity
    changes to maintain hydrostatic balance.

    Hydrostatic equation to be solved:
    dln(p) / dz = -g / (R * Tv)
    """

    ## Load variables in base climate state
    P0 = BC_file['PP'].isel(time=0) + pref
    T0 = BC_file['T'].isel(time=0)
    QV0 = BC_file['QV'].isel(time=0)
    #QC0 = BC_file['QC'].isel(time=0)
    #QI0 = BC_file['QI'].isel(time=0)
    RH0 = specific_to_relative_humidity(QV0, P0, T0)
    TV0 = T0 * (1 + 0.61 * QV0)
    #TV0 = T0 * (1 + 0.61 * QV0 - QC0 -QI0)

    #plt.plot(np.log(P0.mean(dim=['rlon','rlat'])))
    #plt.show()

    ## Load change of variables between climate states
    delta_T = xr.open_dataset(os.path.join(delta_inp_dir,
                    'T{:05d}.nc'.format(delta_time_step)))['T']
    delta_T['level'] = range(len(delta_T.level))
    if 'time' in delta_T.dims:
        delta_T = delta_T.isel(time=0)
    delta_RH = xr.open_dataset(os.path.join(delta_inp_dir,
                    'RELHUM{:05d}.nc'.format(delta_time_step)))['RELHUM']
    delta_RH['level'] = range(len(delta_RH.level))
    if 'time' in delta_RH.dims:
        delta_RH = delta_RH.isel(time=0)

    ## There may be small grid inconsistencies that have to be fixed... 
    fix_grid_coord_diffs(delta_T, BC_file)
    fix_grid_coord_diffs(delta_RH, BC_file)

    ## Compute variables in scenario climate state
    T1 = T0 + delta_T
    QV1 = relative_to_specific_humidity(RH0+delta_RH, P0, T1)
    TV1 = T1 * (1 + 0.61 * QV1)

    ## 
    #lnP_hl = xr.zeros_like(alt_hl)
    P0_hl = xr.ones_like(alt_hl)
    P1_hl = xr.ones_like(alt_hl)
    # extrapolate P to obtain P_hl at model top
    P0_hl.loc[0,:] = (
            P0.sel(level=0) +
            (P0.sel(level=0) - P0.sel(level=1)) /
            (alt.sel(level=0) - alt.sel(level=1)) *
            (alt_hl.sel(level1=0) - alt.sel(level=0))
            )
    # boundary condition: set model top pressure in future
    # to model top pressure of base climate
    P1_hl.loc[0,:] = (
            P0.sel(level=0) +
            (P0.sel(level=0) - P0.sel(level=1)) /
            (alt.sel(level=0) - alt.sel(level=1)) *
            (alt_hl.sel(level1=0) - alt.sel(level=0))
            )
    # compute ln(p)
    lnP0 = np.log(P0)
    lnP0_hl = np.log(P0_hl)
    lnP1_hl = np.log(P1_hl)

    # compute altitude change accross full level (dz)
    dz = alt_hl.diff('level1').rename({'level1':'level'})

    for lvi in range(len(dz.level)):
        ## vertically integrate pressure

        # base state
        dlnP = (
                - CON_G / (CON_RD * TV0.sel(level=lvi)) * dz.sel(level=lvi)
        )
        lnP0_hl.loc[lvi+1,:] = lnP0_hl.loc[lvi,:] + dlnP

        rho0 = P0.sel(level=lvi) / (CON_RD * TV0.sel(level=lvi))
        dP = - rho0 * CON_G * dz.sel(level=lvi)
        P0_hl.loc[lvi+1,:] = P0_hl.loc[lvi,:] + dP

        # future state
        dlnP = (
                - CON_G / (CON_RD * TV1.sel(level=lvi)) * dz.sel(level=lvi)
        )
        lnP1_hl.loc[lvi+1,:] = lnP1_hl.loc[lvi,:] + dlnP

        rho1 = P0.sel(level=lvi) / (CON_RD * TV1.sel(level=lvi))
        dP = - rho1 * CON_G * dz.sel(level=lvi)
        P1_hl.loc[lvi+1,:] = P1_hl.loc[lvi,:] + dP


    # compute ln(P) at full levels
    lnP0_comp = (
            lnP0_hl.sel(level1=range(0,len(lnP0_hl.level1)-1)) +
            lnP0_hl.diff('level1') / 2
    ).rename({'level1':'level'})
    lnP1_comp = (
            lnP1_hl.sel(level1=range(0,len(lnP1_hl.level1)-1)) +
            lnP1_hl.diff('level1') / 2
    ).rename({'level1':'level'})

    P0_comp_dir = (
            P0_hl.sel(level1=range(0,len(P0_hl.level1)-1)) +
            P0_hl.diff('level1') / 2
    ).rename({'level1':'level'})
    P1_comp_dir = (
            P1_hl.sel(level1=range(0,len(P1_hl.level1)-1)) +
            P1_hl.diff('level1') / 2
    ).rename({'level1':'level'})

    #print(lnP0.mean(dim=['rlon', 'rlat']))
    #print(lnP0_comp.mean(dim=['rlon', 'rlat']))
    #plt.plot(lnP0.mean(dim=['rlon', 'rlat']))
    #plt.plot(lnP0_comp.mean(dim=['rlon', 'rlat']))
    #plt.show()
    
    # convert back to real pressure
    P0 = np.exp(lnP0)
    P0_comp_ln = np.exp(lnP0_comp)
    P1_comp_ln = np.exp(lnP1_comp)

    #print(P0.mean(dim=['rlon', 'rlat']))
    #print(P0_comp.mean(dim=['rlon', 'rlat']))
    #plt.plot(P0.mean(dim=['rlon', 'rlat']), alt.mean(dim=['rlon', 'rlat']))
    #plt.plot(P0_comp_ln.mean(dim=['rlon', 'rlat']), alt.mean(dim=['rlon', 'rlat']))
    #plt.plot(P0_comp_dir.mean(dim=['rlon', 'rlat']), alt.mean(dim=['rlon', 'rlat']))
    #plt.show()

    print((P0_comp_dir-P0).mean(dim=['rlon', 'rlat']))
    print((P0_comp_ln-P0).mean(dim=['rlon', 'rlat']))
    plt.plot((P0_comp_dir-P0).mean(dim=['rlon', 'rlat']), alt.mean(dim=['rlon', 'rlat']))
    plt.plot((P0_comp_ln-P0).mean(dim=['rlon', 'rlat']), alt.mean(dim=['rlon', 'rlat']))
    plt.show()
    quit()

    delta_P_dir = P1_comp_dir-P0_comp_dir
    delta_P_ln = P1_comp_ln-P0_comp_ln

    #plt.plot(delta_P_dir.mean(dim=['rlon', 'rlat']), alt.mean(dim=['rlon', 'rlat']))
    #plt.plot(delta_P_ln.mean(dim=['rlon', 'rlat']), alt.mean(dim=['rlon', 'rlat']))
    #plt.show()
    #quit()

    # replace PP in BC_file
    BC_file['PP'].data = BC_file['PP'] + delta_P_dir
    #print(BC_file['PP'])

    return(BC_file)



#def adjust_pressure_to_new_climate_OLD(delta_inp_dir, delta_time_step, 
#                    BC_file, pref, alt_half_level, alt_full_level):
#    """
#    Function to compute pressure change from temperature and humidity
#    changes to maintain hydrostatic balance.
#
#    Equation to be solved:
#    (hydrostatic equation and equation of state)
#    dΔp = -Δρ * g * dz
#    where Δ is the change between two climate states
#    and p0 is the pressure in the base climate state.
#    """
#
#    ## Load variables in base climate state
#    P0 = BC_file['PP'].isel(time=0) + pref
#    T0 = BC_file['T'].isel(time=0)
#    QV0 = BC_file['QV'].isel(time=0)
#    RH0 = specific_to_relative_humidity(QV0, P0, T0)
#    TV0 = T0 * (1 + 0.61 * QV0)
#
#    ## Load change of variables between climate states
#    dT = xr.open_dataset(os.path.join(delta_inp_dir,
#                    'T{:05d}.nc'.format(delta_time_step)))['T']
#    dRH = xr.open_dataset(os.path.join(delta_inp_dir,
#                    'RELHUM{:05d}.nc'.format(delta_time_step)))['RELHUM']
#
#    ## There may be small grid inconsistencies that have to be fixed... 
#    fix_grid_coord_diffs(dT, BC_file)
#    fix_grid_coord_diffs(dRH, BC_file)
#
#    ## Compute variables in scenario climate state
#    T1 = T0 + dT
#    QV1 = relative_to_specific_humidity(RH0+dRH, P0, T1)
#    TV1 = T1 * (1 + 0.61 * QV1)
#
#    # pressure change integration accross full model levels
#    drho = xr.zeros_like(alt_full_level)
#    dP_hl = xr.zeros_like(alt_half_level)
#
#    # compute altitude change accross full level (dz)
#    dz_full_level = alt_half_level.diff('level1').rename({'level1':'level'})
#
#    for lvi in range(len(dz_full_level.level)):
#        #print(lvi)
#        # compute density change between climate states
#        drho.loc[lvi,:] = (
#                P0.isel(level=lvi) / (CON_RD * TV1.isel(level=lvi)) -
#                P0.isel(level=lvi) / (CON_RD * TV0.isel(level=lvi))
#        )
#        ## vertically integrate pressure change between climate states
#        dP_hl.loc[lvi+1,:] = (
#            dP_hl.loc[lvi,:] -
#            drho.isel(level=lvi) * CON_G * dz_full_level.isel(level=lvi)
#        )
#    #plt.plot(drho.mean(dim=['rlon','rlat']), alt_full_level.mean(dim=['rlon','rlat']))
#    #plt.show()
#    #plt.plot(dP_hl.mean(dim=['rlon','rlat']), alt_half_level.mean(dim=['rlon','rlat']))
#    #plt.show()
#    #quit()
#
#    ### interpolate pressure from half levels to full levels
#    ## Simple linear interpolations
#    dP = dP_hl.isel(level1=range(len(dP_hl)-1)) + dP_hl.diff('level1')/2
#    dP = dP.rename({'level1':'level'})
#    ## This would be the more precise way of interpolating. But takes very long!!
#    #dP_good = xr.apply_ufunc(
#    #    interp_cubic,
#    #    alt_full_level,
#    #    alt_half_level,
#    #    dP_hl,
#    #    input_core_dims=[["level"], ["level1"], ["level1"]],  # list with one entry per arg
#    #    output_core_dims=[["level"]],  # returned data has one dimension
#    #    exclude_dims=set(("level1",)),  # dimensions allowed to change size. Must be a set!
#    #    vectorize=True,  # loop over non-core dims
#    #)
#
#    #plt.plot(dP_hl.mean(dim=['rlon','rlat']), alt_half_level.mean(dim=['rlon','rlat']))
#    #plt.plot(dP.mean(dim=['rlon','rlat']), alt_full_level.mean(dim=['rlon','rlat']))
#    ##plt.plot(dP_good.mean(dim=['rlon','rlat']), alt_full_level.mean(dim=['rlon','rlat']))
#    #plt.show()
#    ##plt.plot((dP-dP_good).mean(dim=['rlon','rlat']), alt_full_level.mean(dim=['rlon','rlat']))
#    ##plt.show()
#    #quit()
#
#    #plt.plot((BC_file['PP']).mean(dim=['rlon','rlat','time']), alt_full_level.mean(dim=['rlon','rlat']))
#    #plt.plot((BC_file['PP'].isel(time=0)+dP).mean(dim=['rlon','rlat']), alt_full_level.mean(dim=['rlon','rlat']))
#    #plt.show()
#    #quit()
#
#    # replace PP in BC_file
#    BC_file['PP'].data = BC_file['PP'] + dP
#
#    return(BC_file)








####################################################################################
## Original functions from Roman Brogli
####################################################################################
####################################################################################



def pressure_recompute(laf_file, pref, height_array, height_flat):
    """
    Original version by Roman Brogli.
    """
    #function to compute pressure field in a differen climate using the barometric
    #formula (maintaining hydrostatic balance)
    #temperature changes
    dT_sfc = xr.open_dataset(f'{delta_inp_path}/T_S{delta_time_step:05d}.nc')['T_S']
    dT_atmos = xr.open_dataset(f'{delta_inp_path}/T{delta_time_step:05d}.nc')['T']
    ## HCH2021 start
    ## There are small grid inconsistencies that have to be fixed... 
    fix_grid_coord_diffs(dT_sfc, laf_file)
    fix_grid_coord_diffs(dT_atmos, laf_file)
    ## HCH 2021 stop

    #get pressure field
    pressure_original = laffile['PP'] + pref
    pressure_new = pressure_original.copy()

    #get height difference between model levels
    dz = height_array[:-1] - height_array[1:]
    #print(dz.shape)
    #quit()

    temperature = laffile['T']
    sfc_temperature = laffile['T_S']

    #define barometric height formula
    def barometric(reference_pressure, reference_temperature, dz, lapse_rate):
        R = 8.3144598 #universal gas constant
        M = 0.0289644 # molar mass of air #standard lapse rate
        g = 9.80665
        #lapse_rate = - 0.0065
        exo = - g * M / (R * lapse_rate) #exponent in barometric formula

        #print(reference_pressure.shape)
        #print(reference_temperature.shape)
        #print(dz)
        #print('stop')
        #quit()
        pressure = reference_pressure * ( (reference_temperature + \
        (lapse_rate * dz)) / reference_temperature )**exo

        return pressure

    #print('start')
    #print(temperature.shape)
    #print(height_flat.shape)

    #print(temperature[:,-1,:,:].shape)
    #print(height_flat[-1].shape)

    #compute surface pressure
    surface_press = barometric(pressure_original[:,-1,:,:], temperature[:,-1,:,:], -height_flat[-1], 0.0065)

    #print('start')
    #print(sfc_temperature.shape)
    #print(dT_sfc.shape)
    #print((sfc_temperature+dT_sfc).shape)
    #quit()

    #get the lowest model level in warmer climate
    pressure_new[:,-1,:,:] = barometric(surface_press, sfc_temperature+dT_sfc, height_flat[-1], -0.0065)
    #get the rest (loop from ground up)
    for level in range(len(dz)-1, -1, -1):
            pressure_new[:,level,:,:] = barometric(pressure_new[:,level+1,:,:], \
            temperature[:,level+1,:,:]+dT_atmos[level+1,:,:], dz[level,:,:], -0.0065)

    new_pp = pressure_new.data - pref
    #convert to PP
    laffile['PP'].data = new_pp.astype('float32')

    return laffile




def get_pref_old(vcflat, terrainpath, height_flat):
    """
    Reference pressure for COSMO as developped by Roman Brogli.
    """
    smoothing = (vcflat - height_flat) / vcflat
    smoothing = np.where(smoothing > 0, smoothing, 0)

    const = xr.open_dataset(terrainpath)
    hsurf = const['HSURF'].squeeze()

    # the height at which the reference pressure needs to be computed needs to be 
    # derived form the terrain following coordinates:
    newheights = np.zeros((len(height_flat), hsurf.shape[0], hsurf.shape[1]))

    newheights = (
            height_flat.values[:,None,None] + 
            hsurf.values[None,:,:] * smoothing[:,None,None]
    )

    #New formulation as researched by Christian Steger (untested)
    # Constants
    p0sl = height_flat.p0sl # sea-level pressure [Pa]
    t0sl = height_flat.t0sl   # sea-level temperature [K]
    # Source: COSMO description Part I, page 29
    g = 9.80665     # gravitational acceleration [m s-2]
    R_d = 287.05    # gas constant for dry air [J K-1 kg-1]
    # Source: COSMO source code, data_constants.f90

    # irefatm = 2
    delta_t = height_flat.delta_t
    h_scal = height_flat.h_scal
    # Source: COSMO description Part VII, page 66
    t00 = t0sl - delta_t

    pref = p0sl * np.exp (-g / R_d * h_scal / t00 * \
               np.log((np.exp(newheights / h_scal) * t00 + delta_t) / \
                      (t00 + delta_t)) )
    pref_sfc = p0sl * np.exp (-g / R_d * h_scal / t00 * \
               np.log((np.exp(hsurf.data / h_scal) * t00 + delta_t) / \
                      (t00 + delta_t)) )

    return pref, pref_sfc, newheights
