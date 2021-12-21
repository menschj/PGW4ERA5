#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
description     Auxiliary functions for Postprocess_CLM
author		Christoph Heim based on developments by Roman Brogli
date created    16.12.2021
"""
##############################################################################
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
##############################################################################

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




def get_pref(vcflat, terrainpath, height_flat):
    """
    Reference pressure for COSMO as developped by Roman Brogli.
    """
    smoothing = (vcflat - height_flat) / vcflat
    smoothing = np.where(smoothing > 0, smoothing, 0)

    const = xr.open_dataset(terrainpath)
    hsurf = const['HSURF'].squeeze()

    #the height at which the reference pressure needs to be computed needs to be derived form the terrain   following coordinates:
    newheights = np.zeros((len(height_flat), hsurf.shape[0], hsurf.shape[1]))

    #avoid forloop
    newheights = height_flat.values[:,None,None] + hsurf.values[None,:,:] * smoothing[:,None,None]

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


def fix_grid_coord_diffs(fix_ds, ref_ds):
    ### There may be small grid inconsistencies that have to be fixed... 
    ### the reason is that the cdo remapping produes them.
    if np.max(np.abs(fix_ds.rlon.values - ref_ds.rlon.values)) > 1E-5:
        raise ValueError('Grids differ!')
    if np.mean(np.abs(fix_ds.rlon.values - ref_ds.rlon.values)) > 1E-8:
        print('fix small differences in inputs grids.')
        fix_ds['rlon'] = ref_ds.rlon
        fix_ds['rlat'] = ref_ds.rlat


def get_alt_half_level(vcoord, terrainpath):
    """
    Compute altitude of half levels (Gal-Chen coordinate)
    """
    height_flat_half_level = vcoord #these are half levels

    smoothing = (vcoord.vcflat - height_flat_half_level) / vcoord.vcflat
    smoothing = smoothing.where(smoothing > 0, 0)

    const = xr.open_dataset(terrainpath)
    hsurf = const['HSURF'].squeeze()

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



def interp_cubic(x_out, x_in, data_in):
    f = interp1d(x_in, data_in, kind='cubic')
                #fill_value='extrapolate', bounds_error=False)
    out = f(x_out)
    return(out)



def adjust_pressure_to_new_climate(delta_inp_dir, delta_time_step, 
                    BC_file, pref, alt_half_level, alt_full_level):
    """
    Function to compute pressure change from temperature and humidity
    changes to maintain hydrostatic balance.

    Equation to be solved:
    (hydrostatic equation and equation of state)
    dΔp = -Δρ * g * dz
    dΔp = - p0/(Rd * ΔTv) * g * dz
    where Δ is the change between two climate states
    and p0 is the pressure in the base climate state.
    """

    # gas constant for dry air
    CON_RD = 287.06 # [J kg-1 K-1]
    # gravity constant (assumed constant over the entire profile)
    CON_G = 9.81 # [m s-2]

    ## Load variables in base climate state
    P0 = BC_file['PP'].isel(time=0) + pref
    T0 = BC_file['T'].isel(time=0)
    QV0 = BC_file['QV'].isel(time=0)
    RH0 = specific_to_relative_humidity(QV0, P0, T0)
    TV0 = T0 * (1 + 0.61 * QV0)

    ## Load change of variables between climate states
    dT = xr.open_dataset(os.path.join(delta_inp_dir,
                    'T{:05d}.nc'.format(delta_time_step)))['T']
    dRH = xr.open_dataset(os.path.join(delta_inp_dir,
                    'RELHUM{:05d}.nc'.format(delta_time_step)))['RELHUM']

    ## There may be small grid inconsistencies that have to be fixed... 
    fix_grid_coord_diffs(dT, BC_file)
    fix_grid_coord_diffs(dRH, BC_file)

    ## Compute variables in scenario climate state
    T1 = T0 + dT
    QV1 = relative_to_specific_humidity(RH0+dRH, P0, T1)
    TV1 = T1 * (1 + 0.61 * QV1)

    # pressure change integration accross full model levels
    drho = xr.zeros_like(alt_full_level)
    dP_hl = xr.zeros_like(alt_half_level)

    # compute altitude change accross full level (dz)
    dz_full_level = alt_half_level.diff('level1').rename({'level1':'level'})

    for lvi in range(len(dz_full_level.level)):
        #print(lvi)
        # compute density change between climate states
        drho.loc[lvi:] = (
                P0.isel(level=lvi) / (CON_RD * TV1.isel(level=lvi)) -
                P0.isel(level=lvi) / (CON_RD * TV0.isel(level=lvi))
        )
        ## vertically integrate pressure change between climate states
        dP_hl.loc[lvi+1,:] = (
            dP_hl.loc[lvi,:] -
            drho.isel(level=lvi) * CON_G * dz_full_level.isel(level=lvi)
        )
    #plt.plot(drho.mean(dim=['rlon','rlat']), alt_full_level.mean(dim=['rlon','rlat']))
    #plt.show()
    #plt.plot(dP_hl.mean(dim=['rlon','rlat']), alt_half_level.mean(dim=['rlon','rlat']))
    #plt.show()
    #quit()

    ### interpolate pressure from half levels to full levels
    ## Simple linear interpolations
    dP = dP_hl.isel(level1=range(len(dP_hl)-1)) + dP_hl.diff('level1')/2
    dP = dP.rename({'level1':'level'})
    ## This would be the more precise way of interpolating. But takes very long!!
    #dP_good = xr.apply_ufunc(
    #    interp_cubic,
    #    alt_full_level,
    #    alt_half_level,
    #    dP_hl,
    #    input_core_dims=[["level"], ["level1"], ["level1"]],  # list with one entry per arg
    #    output_core_dims=[["level"]],  # returned data has one dimension
    #    exclude_dims=set(("level1",)),  # dimensions allowed to change size. Must be a set!
    #    vectorize=True,  # loop over non-core dims
    #)

    #plt.plot(dP_hl.mean(dim=['rlon','rlat']), alt_half_level.mean(dim=['rlon','rlat']))
    #plt.plot(dP.mean(dim=['rlon','rlat']), alt_full_level.mean(dim=['rlon','rlat']))
    ##plt.plot(dP_good.mean(dim=['rlon','rlat']), alt_full_level.mean(dim=['rlon','rlat']))
    #plt.show()
    ##plt.plot((dP-dP_good).mean(dim=['rlon','rlat']), alt_full_level.mean(dim=['rlon','rlat']))
    ##plt.show()
    #quit()

    #plt.plot((BC_file['PP']).mean(dim=['rlon','rlat','time']), alt_full_level.mean(dim=['rlon','rlat']))
    #plt.plot((BC_file['PP'].isel(time=0)+dP).mean(dim=['rlon','rlat']), alt_full_level.mean(dim=['rlon','rlat']))
    #plt.show()
    #quit()

    # replace PP in BC_file
    BC_file['PP'].data = BC_file['PP'] + dP

    return(BC_file)
