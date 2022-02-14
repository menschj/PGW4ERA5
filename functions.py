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


def load_delta(delta_inp_path, var_name, delta_date_time, 
                laf_time, diff_time_step,
                name_base='{}_delta.nc'):

    delta_year = xr.open_dataset(os.path.join(delta_inp_path,
                            name_base.format(var_name)))
    if delta_date_time is not None:
        # replace delta year values with year of current delta_date_time
        for i in range(len(delta_year.time)):
            delta_year.time.values[i] = dt64_to_dt(
                        delta_year.time[i]).replace(year=delta_date_time.year)
        # interpolate in time and select variable
        delta = delta_year[var_name].interp(time=delta_date_time, 
                                    method='linear', 
                                ).expand_dims(dim='time', axis=0)
        # make sure time is in the same format as in laf file
        delta['time'] = laf_time
    else:
        delta = delta_year[var_name]

    return(delta)


def get_delta_era5(var, var_name, laffile, delta_inp_path,
                    diff_time_step, date_time=None):
    delta = load_delta(delta_inp_path, var_name, date_time, laffile.time,
                        diff_time_step)
    #delta = delta.assign_coords({'lat':var.lat.values})
    return(delta)


def get_delta_interp_era5(var, target_P, var_name, laffile, delta_inp_path,
                        diff_time_step, date_time):
    delta = load_delta(delta_inp_path, var_name, date_time, laffile.time,
                        diff_time_step)
    #delta = delta.assign_coords({'lat':var.lat.values})

    if var_name in ['ta','hur']:
        sfc_var_name = var_name + 's'
        delta_sfc = load_delta(delta_inp_path, sfc_var_name, 
                            date_time, laffile.time,
                            diff_time_step)
        #delta_sfc = delta_sfc.assign_coords({'lat':var.lat.values})
        ps_hist = load_delta(delta_inp_path, 'ps', 
                            date_time, laffile.time,
                            diff_time_step,
                            name_base='{}_historical.nc')
        #ps_hist = ps_hist.assign_coords({'lat':var.lat.values})
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

    ##print(delta)
    #if delta.name in ['hus', 'QV']:
    #    print('WARNING: DEBUG MODE FOR variable hus model top!!!')
    #    top = xr.zeros_like(delta.sel(plev=100000))
    #    top['plev'] = 500
    #    delta = xr.concat([delta, top], dim='plev').transpose(
    #                                'time', 'plev', 'lat', 'lon')
    #    #print(delta.mean(dim=['lat','lon','time']))
    #    #print(delta.plev)
    #    #quit()

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

    if np.min(target_P) < np.min(source_P):
        raise ValueError('ERA5 top pressure is lower than '+
                         'climate delta top pressure!')

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
            raise ValueError('Extrapolation deactivated but data out of bounds.')

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


def determine_p_ref(p_era, p_pgw, p_ref_last, p_ref_opts):
    for p in p_ref_opts:
        if (p_era > p) & (p_pgw > p):
        #if x / p > 1.01:
            if p_ref_last is None:
                return(p)
            else:
                return(min(p, p_ref_last))


def integ_geopot(P_hl, FIS, T, QV, level1, p_ref):
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

        ## geopotential at full level
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

