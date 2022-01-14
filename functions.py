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
from constants import CON_RD, CON_G
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


def fix_grid_coord_diffs(fix_ds, ref_ds):
    ### There may be small grid inconsistencies that have to be fixed... 
    ### the reason is that the cdo remapping produes them.
    if np.max(np.abs(fix_ds.rlon.values - ref_ds.rlon.values)) > 1E-5:
        raise ValueError('Grids differ!')
    if np.mean(np.abs(fix_ds.rlon.values - ref_ds.rlon.values)) > 1E-8:
        print('fix small differences in inputs grids.')
        fix_ds['rlon'] = ref_ds.rlon
        fix_ds['rlat'] = ref_ds.rlat


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



def get_pref(vcoord, terrainpath, hl_or_fl):
    """
    Compute reference pressure for COSMO.
    hl_or_fl:   if ='hl', reference pressure at half levels is computed.
                if ='fl', reference pressure at full levels is computed.
    """
    # altitude of model half levels
    if hl_or_fl == 'hl':
        alt = get_alt_half_level(vcoord, terrainpath)
    elif hl_or_fl == 'fl':
        alt = get_alt_half_level(vcoord, terrainpath)
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




def interp_cubic(x_out, x_in, data_in):
    f = interp1d(x_in, data_in, kind='cubic')
                #fill_value='extrapolate', bounds_error=False)
    out = f(x_out)
    return(out)




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

    ## 
    #lnP_hl = xr.zeros_like(alt_hl)
    P0_hl = xr.ones_like(alt_hl)
    P1_hl = xr.ones_like(alt_hl)
    # extrapolate P to obtain P_hl at model top
    P0_hl.loc[0,:] = (
            P0.isel(level=0) +
            (P0.isel(level=0) - P0.isel(level=1)) /
            (alt.isel(level=0) - alt.isel(level=1)) *
            (alt_hl.isel(level1=0) - alt.isel(level=0))
            )
    # boundary condition: set model top pressure in future
    # to model top pressure of base climate
    P1_hl.loc[0,:] = (
            P0.isel(level=0) +
            (P0.isel(level=0) - P0.isel(level=1)) /
            (alt.isel(level=0) - alt.isel(level=1)) *
            (alt_hl.isel(level1=0) - alt.isel(level=0))
            )
    # compute ln(p)
    lnP0 = np.log(P0)
    lnP0_hl = np.log(P0_hl)
    lnP1_hl = np.log(P1_hl)

    # compute altitude change accross full level (dz)
    dz = alt_hl.diff('level1').rename({'level1':'level'})

    for lvi in range(len(dz.level)):
        #print(lvi)

        ## vertically integrate pressure

        # base state
        dlnP = (
                - CON_G / (CON_RD * TV0.isel(level=lvi)) * dz.isel(level=lvi)
        )
        lnP0_hl.loc[lvi+1,:] = lnP0_hl.loc[lvi,:] + dlnP

        rho0 = P0.isel(level=lvi) / (CON_RD * TV0.isel(level=lvi))
        dP = - rho0 * CON_G * dz.isel(level=lvi)
        P0_hl.loc[lvi+1,:] = P0_hl.loc[lvi,:] + dP

        # future state
        dlnP = (
                - CON_G / (CON_RD * TV1.isel(level=lvi)) * dz.isel(level=lvi)
        )
        lnP1_hl.loc[lvi+1,:] = lnP1_hl.loc[lvi,:] + dlnP

        rho1 = P0.isel(level=lvi) / (CON_RD * TV1.isel(level=lvi))
        dP = - rho1 * CON_G * dz.isel(level=lvi)
        P1_hl.loc[lvi+1,:] = P1_hl.loc[lvi,:] + dP


    # compute ln(P) at full levels
    lnP0_comp = (
            lnP0_hl.isel(level1=range(0,len(lnP0_hl.level1)-1)) +
            lnP0_hl.diff('level1') / 2
    ).rename({'level1':'level'})

    P0_comp2 = (
            P0_hl.isel(level1=range(0,len(P0_hl.level1)-1)) +
            P0_hl.diff('level1') / 2
    ).rename({'level1':'level'})

    lnP1_comp = (
            lnP1_hl.isel(level1=range(0,len(lnP1_hl.level1)-1)) +
            lnP1_hl.diff('level1') / 2
    ).rename({'level1':'level'})

    P1_comp2 = (
            P1_hl.isel(level1=range(0,len(P1_hl.level1)-1)) +
            P1_hl.diff('level1') / 2
    ).rename({'level1':'level'})

    #print(lnP0.mean(dim=['rlon', 'rlat']))
    #print(lnP0_comp.mean(dim=['rlon', 'rlat']))
    #plt.plot(lnP0.mean(dim=['rlon', 'rlat']))
    #plt.plot(lnP0_comp.mean(dim=['rlon', 'rlat']))
    #plt.show()
    
    # convert back to real pressure
    P0 = np.exp(lnP0)
    P0_comp = np.exp(lnP0_comp)
    P1_comp = np.exp(lnP1_comp)

    #print(P0.mean(dim=['rlon', 'rlat']))
    #print(P0_comp.mean(dim=['rlon', 'rlat']))
    #plt.plot(P0.mean(dim=['rlon', 'rlat']), alt.mean(dim=['rlon', 'rlat']))
    #plt.plot(P0_comp2.mean(dim=['rlon', 'rlat']), alt.mean(dim=['rlon', 'rlat']))
    #plt.plot(P1_comp2.mean(dim=['rlon', 'rlat']), alt.mean(dim=['rlon', 'rlat']))
    #plt.show()

    #print((P0_comp-P0).mean(dim=['rlon', 'rlat']))
    #print((P0_comp2-P0).mean(dim=['rlon', 'rlat']))
    #plt.plot((P0_comp2-P0).mean(dim=['rlon', 'rlat']), alt.mean(dim=['rlon', 'rlat']))
    #plt.plot((P0_comp-P0).mean(dim=['rlon', 'rlat']), alt.mean(dim=['rlon', 'rlat']))
    #plt.show()

    #plt.plot((P1_comp2-P0_comp2).mean(dim=['rlon', 'rlat']), alt.mean(dim=['rlon', 'rlat']))
    #plt.plot((P1_comp-P0_comp).mean(dim=['rlon', 'rlat']), alt.mean(dim=['rlon', 'rlat']))
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









def adjust_pressure_to_new_climate_OLD(delta_inp_dir, delta_time_step, 
                    BC_file, pref, alt_half_level, alt_full_level):
    """
    Function to compute pressure change from temperature and humidity
    changes to maintain hydrostatic balance.

    Equation to be solved:
    (hydrostatic equation and equation of state)
    dΔp = -Δρ * g * dz
    where Δ is the change between two climate states
    and p0 is the pressure in the base climate state.
    """

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
        drho.loc[lvi,:] = (
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













def pressure_recompute(laf_file, pref, height_array, height_flat):
    """
    Original version by Roman Brogli.
    """
    #function to compute pressure field in a differen climate using the barometric
    #formula (maintaining hydrostatic balance)
    #temperature changes
    dT_sfc = xr.open_dataset(f'{Diffspath}/T_S{laftimestep:05d}.nc')['T_S']
    dT_atmos = xr.open_dataset(f'{Diffspath}/T{laftimestep:05d}.nc')['T']
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
