#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
description     Auxiliary functions for PGW for ERA5
author		Christoph Heim based on original developments by Roman Brogli
date created    12.01.2022
"""
##############################################################################
import os
import xarray as xr
import numpy as np
from numba import njit
from datetime import datetime,timedelta
from constants import CON_RD, CON_G
from settings import *
##############################################################################


##############################################################################
##### ARBITRARY FUNCTIONS
##############################################################################
def dt64_to_dt(date):
    """
    Converts a numpy datetime64 object to a python datetime object 
    Input:
      date - a np.datetime64 object
    Output:
      DATE - a python datetime object
    source: https://gist.github.com/blaylockbk/1677b446bc741ee2db3e943ab7e4cabd
    """
    timestamp = ((date - np.datetime64('1970-01-01T00:00:00'))
                 / np.timedelta64(1, 's'))
    return datetime.utcfromtimestamp(timestamp)


##############################################################################
##### PHYSICAL COMPUTATIONS
##############################################################################
def specific_to_relative_humidity(hus, pa, ta):
    """
    Compute relative humidity from specific humidity.
    """
    hur = 0.263 * pa * hus *(np.exp(17.67*(ta - 273.15)/(ta-29.65)))**(-1)
    return(hur)


def relative_to_specific_humidity(hur, pa, ta):
    """
    Compute specific humidity from relative humidity.
    """
    hus = (hur  * np.exp(17.67 * (ta - 273.15)/(ta - 29.65))) / (0.263 * pa)
    return(hus)


def integ_geopot(pa_hl, zgs, ta, hus, level1, p_ref):
    """
    Integrate ERA5 geopotential from surfce to a reference pressure
    level p_ref.
    """
    # take log half-level pressure difference (located at full levels)
    dlnpa = np.log(pa_hl).diff(
                dim=VERT_HL_ERA, 
                label='lower').rename({VERT_HL_ERA:VERT_ERA})

    # create geopotential array and fill with surface geopotential
    phi_hl = zgs.expand_dims(dim={VERT_HL_ERA:level1}).copy()

    # compute virtual temperature
    tav = ta * (1 + 0.61 * hus)

    ## integrate over model half levels
    for l in sorted(tav[VERT_ERA].values, reverse=True):
        # geopotential at full level
        phi_hl.loc[{VERT_HL_ERA:l}] = (
                phi_hl.sel({VERT_HL_ERA:l+1}) +
                (CON_RD * tav.sel({VERT_ERA:l}) * dlnpa.sel({VERT_ERA:l}))
        )

            
    phi_hl = phi_hl.transpose(TIME_ERA, VERT_HL_ERA, LAT_ERA, LON_ERA)

    ## integrate from last half level below reference pressure
    ## up to reference pressure
    # determine level below reference pressure
    p_diff = pa_hl - p_ref
    p_diff = p_diff.where(p_diff >= 0, np.nan)
    ind_ref_star = p_diff.argmin(dim=VERT_HL_ERA)
    hl_ref_star = p_diff[VERT_HL_ERA].isel({VERT_HL_ERA:ind_ref_star})
    # get pressure and geopotential of that level
    p_ref_star = pa_hl.sel({VERT_HL_ERA:hl_ref_star})
    phi_ref_star = phi_hl.sel({VERT_HL_ERA:hl_ref_star})

    # finally interpolate geopotential to reference
    # pressure level
    phi_ref = (
            phi_ref_star -
            (CON_RD * tav.sel({VERT_ERA:hl_ref_star-1})) * 
            (np.log(p_ref) - np.log(p_ref_star))
    )

    # remove multi-dimensional coordinates
    if VERT_HL_ERA in phi_ref.coords:
        del phi_ref[VERT_HL_ERA]
    if VERT_ERA in phi_ref.coords:
        del phi_ref[VERT_ERA]
    if PLEV_GCM in phi_ref.coords:
        del phi_ref[VERT_ERA]

    return(phi_ref)


##############################################################################
##### CLIMATE DELTA COMPUTATION AND INTERPOLATION
##############################################################################
def load_delta(delta_inp_path, var_name, era_date_time, 
               delta_date_time=None, name_base='{}_delta.nc'):
    """
    Load a climate delta and if delta_date_time is given,
    interpolate it to that date and time of the year.
    """
    ## full climate delta (either daily or monthly)
    full_delta = xr.open_dataset(os.path.join(delta_inp_path,
                            name_base.format(var_name)))

    ## if climate delta should be interpolated to a specific time
    if delta_date_time is not None:
        # replace delta year values with year of current delta_date_time
        for i in range(len(full_delta.time)):
            full_delta.time.values[i] = dt64_to_dt(
                        full_delta.time[i]).replace(
                                year=delta_date_time.year)
        # interpolate in time and select variable
        delta = full_delta[var_name].interp(time=delta_date_time, 
                                    method='linear', 
                                ).expand_dims(dim='time', axis=0)

        # make sure time is in the same format as in ERA5 file
        # ERA5 has "seconds since xyz" while delta has np.datetime64
        delta['time'] = era_date_time

    ## if full climate delta should be returned without 
    ## time interpolation
    else:
        delta = full_delta[var_name]

    return(delta)


def load_delta_interp(delta_inp_path, var_name, target_P,
                        era_date_time, delta_date_time):
    """
    Does the following:
        - load a climate delta
        - for specific variables (ta and hur) also load surface value
          as well as historical surface pressure. This is to extend
          the 3D climate deltas with surface values which makes
          the interpolation to the ERA5 model levels more precise.
        - vertically interpolate climate deltas to ERA5 model levels
    """
    delta = load_delta(delta_inp_path, var_name, era_date_time, delta_date_time)

    ## for specific variables also load climate delta for surface
    ## values and the historical surface pressure.
    if var_name in ['ta','hur']:
        sfc_var_name = var_name + 's'
        delta_sfc = load_delta(delta_inp_path, sfc_var_name, 
                            era_date_time, delta_date_time)
        ps_hist = load_delta(delta_inp_path, 'ps', 
                            era_date_time, delta_date_time,
                            name_base='{}_historical.nc')
    else:
        delta_sfc = None
        ps_hist = None

    # interpolate climate delta onto ERA5 model levels
    delta = vert_interp_delta(delta, target_P, delta_sfc, ps_hist)
    return(delta)


def replace_delta_sfc(source_P, ps_hist, delta, delta_sfc):
    """
    In the 3D climate deltas, replace the value just below
    the surface by the surface climate delta value and insert
    it a historical surface pressure. This improves the precision
    of the climate deltas during interpolation to the ERA5 model levels.
    All 3D climate delta values below the historical surface pressure
    are set to the surface value (constant extrapolation). This is
    because within the orography the GCM climate delta is assumed
    to be incorrect.
    """
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
        out_delta[sfc_ind:] = delta_sfc
        out_source_P[sfc_ind] = ps_hist
    return(out_source_P, out_delta)


def vert_interp_delta(delta, target_P, delta_sfc=None, ps_hist=None):
    """
    Vertically interpolate climate delta onto ERA5 model levels.
    If delta_sfc and ps_hist are given, surface values will
    be inserted into the 3D climate delta at the height of
    the surface pressure. This gives a more precise interpolation.
    Climate delta values below the surface are set to the surface
    climate delta because below the surface, the GCM climate delta
    is considered unreliable and thus constant extrapolation
    seems more reasonable.
    """

    # sort delta dataset from top to bottom (pressure ascending)
    delta = delta.reindex(
                {PLEV_GCM:list(reversed(delta[PLEV_GCM]))})

    # create 4D source pressure with GCM pressure levels
    source_P = delta[PLEV_GCM].expand_dims(
                    dim={LON_GCM:delta[LON_GCM],
                         LAT_GCM:delta[LAT_GCM],
                         TIME_GCM:delta[TIME_GCM]}).transpose(
                                TIME_GCM, PLEV_GCM, LAT_GCM, LON_GCM)

    ## if surface values are given, replace them at the
    ## level of the surface pressure
    if delta_sfc is not None:
        source_P, delta = xr.apply_ufunc(
                replace_delta_sfc, source_P, 
                ps_hist, 
                delta, delta_sfc,
                input_core_dims=[[PLEV_GCM],[],[PLEV_GCM],[]],
                output_core_dims=[[PLEV_GCM],[PLEV_GCM]],
                vectorize=True)
        source_P = source_P.transpose(TIME_GCM, PLEV_GCM, LAT_GCM, LON_GCM)
        delta = delta.transpose(TIME_GCM, PLEV_GCM, LAT_GCM, LON_GCM)

    # make sure all arrays contain the required dimensions
    if source_P.dims != (TIME_GCM, PLEV_GCM, LAT_GCM, LON_GCM):
        raise ValueError()
    if delta.dims != (TIME_GCM, PLEV_GCM, LAT_GCM, LON_GCM):
        raise ValueError()
    if target_P.dims != (TIME_ERA, VERT_ERA, LAT_ERA, LON_ERA):
        raise ValueError()

    # make sure there is no extrapolation at the model top
    if np.min(target_P) < np.min(source_P):
        raise ValueError('ERA5 top pressure is lower than '+
                         'climate delta top pressure!')

    # run interpolation
    delta_interp = interp_logp_3d(delta, source_P, target_P,
                        extrapolate='constant')
    return(delta_interp)


def interp_logp_3d(var, source_P, targ_P, extrapolate='off'):
    """
    Interpolate 3D array in vertical (pressure) dimension using the
    logarithm of pressure.
    extrapolate:
        - off: no extrapolation
        - linear: linear extrapolation
        - constant: constant extrapolation
    """
    if extrapolate not in ['off', 'linear', 'constant']:
        raise ValueError()

    targ = xr.zeros_like(targ_P)
    tmp = np.zeros_like(targ.values.squeeze())
    interp_1d_for_latlon(var.values.squeeze(),
                np.log(source_P.values.squeeze()),
                np.log(targ_P.squeeze()).values, 
                tmp,
                len(targ_P[LAT_ERA]), len(targ_P[LON_ERA]),
                extrapolate)
    tmp = np.expand_dims(tmp, axis=0)
    targ.values = tmp
    return(targ)



@njit()
def interp_1d_for_latlon(orig_array, src_p, targ_p, interp_array,
                        nlat, nlon, extrapolate):
    """
    Vertical interpolation helper function with numba njit for 
    fast performance.
    Loop over lat and lon dimensions and interpolate each column
    individually
    extrapolate:
        - off: no extrapolation
        - linear: linear extrapolation
        - constant: constant extrapolation
    """
    for lat_ind in range(nlat):
        for lon_ind in range(nlon):
            src_val_col = orig_array[:, lat_ind, lon_ind]
            src_p_col = src_p[:, lat_ind, lon_ind]
            targ_p_col = targ_p[:, lat_ind, lon_ind]

            # call 1D interpolation function for current column
            interp_col = interp_extrap_1d(src_p_col, src_val_col, 
                                        targ_p_col, extrapolate)
            interp_array[:, lat_ind, lon_ind] = interp_col


@njit()
def interp_extrap_1d(src_x, src_y, targ_x, extrapolate):
    """
    Numba helper function for interpolation of 1d vertical column.
    Does constant extrapolation which is used for the climate deltas.
    extrapolate:
        - off: no extrapolation
        - linear: linear extrapolation
        - constant: constant extrapolation
    """
    targ_y = np.zeros(len(targ_x))
    for ti in range(len(targ_x)):
        i1 = -1
        i2 = -1
        require_extrap = False
        for si in range(len(src_x)):
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

        # raise value if extrapolation is required but not enabled.
        if require_extrap and extrapolate == 'off':
            raise ValueError('Extrapolation deactivated but data out of bounds.')

        # interpolate/extrapolate values
        if i1 == i2:
            targ_y[ti] = src_y[i1]
        else:
            targ_y[ti] = (
                src_y[i1] + (targ_x[ti] - src_x[i1]) * 
                (src_y[i2] - src_y[i1]) / (src_x[i2] - src_x[i1])
            )

    return(targ_y)


def determine_p_ref(ps_era, ps_pgw, p_ref_opts, p_ref_last=None):
    """
    Find lowest GCM pressure level among p_ref_opts that lies above 
    surface (surface pressure) in both ERA and PGW climate.
    Also ensure that during the iterations, no reference pressure level 
    at lower altitude than during last iterations is used. This is to
    prevent the iteration algorithm to oscillate between two reference
    pressure levels and not converge.
    """
    for p in p_ref_opts:
        if (ps_era > p) & (ps_pgw > p):
            if p_ref_last is None:
                return(p)
            else:
                return(min(p, p_ref_last))