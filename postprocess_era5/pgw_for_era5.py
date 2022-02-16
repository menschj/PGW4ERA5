#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
description     PGW for ERA5
authors		    Before 2022: original developments by Roman Brogli
                After 2022:  updates by Christoph Heim 
"""
##############################################################################
import argparse, os
import xarray as xr
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
from base.functions import (
        specific_to_relative_humidity,
        relative_to_specific_humidity,
        load_delta,
        load_delta_interp,
        integ_geopot,
        interp_logp_3d,
        determine_p_ref,
        )
from constants import CON_G, CON_RD
from parallel import IterMP
from settings import *
##############################################################################

def pgw_for_era5(inp_era_file_path, out_era_file_path,
                delta_input_dir, era_step_dt):
    """
    ##########################################################################

    Main function to update ERA5 files with the PGW signal.
    The terminology used is the "ERA climate state" referring to
    the climate state in the ERA5 files, as well as the "PGW 
    climate state" referring to the future (or past) climate state.
    The script adds (and requires) climate deltas for:
        - ua
        - va
        - ta (using tas)
        - hus (computed using a hur and hurs climate delta)
        - surface and soil temperature
    and consequently iteratively updates ps to maintain hydrostatic
    balance. During this, the climate delta for zg is additionally required.
    Finally, the GCM ps from the historical simulation is also needed.

    ##########################################################################

    If the variable names in the ERA5 files to be processed deviate from
    the usual naming convention, the dict 'var_name_map' in the file 
    settings.py allows to map the usual names to the names in the ERA5
    used here. Also the coordinate names in the ERA5 or the GCM climate
    delta files can be changed there, if required.

    ##########################################################################

    Some more information about the iterative surface pressure
    adjustment:

    - As a default option, the climate deltas are interpolated to
    the ERA5 model levels of the ERA climate state (i_reinterp = 0).
    There is an option implemented (i_reinterp = 1) in which the
    deltas are re-interpolated on the updated ERA5 model levels
    with each iteration of surface pressure adjustment. This was
    found to lead more balanced PGW climate states if the climate
    deltas have coarse vertical resolution. However, it also
    implies that the ERA5 fields are extrapolated at the surface
    (if the surface pressure increases) the effect of which was not
    tested in detail. The extrapolation is done assuming that the
    boundary values are constant, which is not ideal for height-dependent
    variables like e.g. temperature.

    - The procedure requires a reference pressure level (e.g. 500 hPa) for
    which the geopotential is computed. Based on the deviation between the
    computed and the GCM reference pressure geopotential, the surface pressure
    is adjusted. Since the climate deltas may not always be available at 
    native vertical GCM resolution, but the climate delta for the geopotential
    on one specific pressure level itself is computed by the GCM using data
    from all GCM model levels below, this introduces an error in the surface
    pressure adjustment used here.
    The higher the reference pressure is chosen, the large this error may
    get. To alleviate this problem, the default option is that the reference
    pressure is determined locally as the lowest possible pressure above
    the surface for which a climate delta for the geopotential is available.
    It is strongly recommended to use this option (p_ref = None).

    - If the iteration does not converge, 'thresh_phi_ref_max_error' in
    the file settings.py may have to be raised a little bit. Setting
    i_debug = 2 may help to diagnose if this helps.

    ##########################################################################

    """
    if i_debug >= 0:
        print('Start working on file {}'.format(inp_era_file_path))

    # open data set
    laffile = xr.open_dataset(inp_era_file_path, decode_cf=False)

    # compute pressure on ERA5 levels and half-levels
    pa_era = (laffile.akm + 
              laffile[var_name_map['ps']] * laffile.bkm).transpose(
                TIME_ERA, VERT_ERA, LAT_ERA, LON_ERA)
    pa_hl_era = (laffile.ak + 
                laffile[var_name_map['ps']] * laffile.bk).transpose(
                TIME_ERA, VERT_HL_ERA, LAT_ERA, LON_ERA)

    # compute relative humidity in ERA climate state
    laffile[var_name_map['hur']] = specific_to_relative_humidity(
                        laffile[var_name_map['hus']], 
                        pa_era, laffile[var_name_map['ta']]).transpose(
                        TIME_ERA, VERT_ERA, LAT_ERA, LON_ERA)

    # update surface temperature using near-surface temperature delta
    var_name = 'tas'
    if i_debug >= 2:
        print('update {}'.format(var_name))
    delta_tas = load_delta(delta_input_dir, var_name,
                            laffile[TIME_ERA], era_step_dt)
    laffile[var_name_map[var_name]].values += delta_tas.values

    # update temperature of soil layers
    var_name = 'st'
    if i_debug >= 2:
        print('update {}'.format(var_name))
    # set climatological lower soil temperature delta to annual mean
    # climate delta of near-surface temperature.
    delta_st_clim = load_delta(delta_input_dir, 'tas',
                            laffile[TIME_ERA], 
                            delta_date_time=None).mean(dim=[TIME_GCM])
    delta_st = (
            delta_st_clim + np.exp(-laffile.soil1/2.8) * 
                    (delta_tas - delta_st_clim)
    )
    delta_st = delta_st.transpose(TIME_ERA, SOIL_HL_ERA, LAT_ERA, LON_ERA)
    laffile[var_name_map[var_name]].values += delta_st

    # containers for variable computation
    vars_era = {}
    vars_pgw = {}
    deltas = {}

    # If no re-interpolation is done, the final PGW climate state
    # variables can be computed already now, before updating the 
    # surface pressure. This means that the climate deltas or
    # interpolated on the ERA5 model levels of the ERA climate state.
    if not i_reinterp:

        ### interpolate climate deltas onto ERA5 grid
        for var_name in ['ta','hur','ua','va']:
            if i_debug >= 2:
                print('update {}'.format(var_name))

            ## interpolate climate deltas to ERA5 model levels
            ## use ERA climate state
            delta_var = load_delta_interp(delta_input_dir,
                    var_name, pa_era, laffile[TIME_ERA], era_step_dt)
            deltas[var_name] = delta_var

            ## compute PGW climate state variables
            vars_pgw[var_name] = (
                    laffile[var_name_map[var_name]] + 
                    deltas[var_name]
            )

            # convert relative humidity to specific humidity
            # take PGW climate state temperature and relative humidity
            # but assume ERA climate state pressure
            if var_name == 'hur':
                vars_pgw['hus'] = relative_to_specific_humidity(
                                vars_pgw['hur'], pa_era, vars_pgw['ta'])


    #########################################################################
    ### UPDATE SURFACE PRESSURE USING ITERATIVE PROCEDURE
    #########################################################################
    # change in surface pressure between ERA and PGW climate states
    delta_ps = xr.zeros_like(laffile[var_name_map['ps']])
    # increment to adjust delta_ps with each iteration
    adj_ps = xr.zeros_like(laffile[var_name_map['ps']])
    # maximum error in geopotential (used in iteration)
    phi_ref_max_error = np.inf

    it = 1
    while phi_ref_max_error > thresh_phi_ref_max_error:

        # update surface pressure
        delta_ps += adj_ps
        ps_pgw = laffile[var_name_map['ps']] + delta_ps

        # recompute pressure on full and half levels
        pa_pgw = (laffile.akm + ps_pgw * laffile.bkm).transpose(
                    TIME_ERA, VERT_ERA, LAT_ERA, LON_ERA)
        pa_hl_pgw = (laffile.ak + ps_pgw * laffile.bk).transpose(
                    TIME_ERA, VERT_HL_ERA, LAT_ERA, LON_ERA)


        if i_reinterp:
            # interpolate ERA climate state variables as well as
            # climate deltas onto updated model levels, and
            # compute PGW climate state variables
            for var_name in ['ta', 'hur']:
                vars_era[var_name] = interp_logp_3d(
                                laffile[var_name_map[var_name]], 
                                pa_era, pa_pgw, extrapolate='constant')
                deltas[var_name] = load_delta_interp(delta_input_dir,
                                                var_name, pa_pgw,
                                                laffile[TIME_ERA], era_step_dt)
                vars_pgw[var_name] = vars_era[var_name] + deltas[var_name]

            # convert relative humidity to speicifc humidity in pgw
            vars_pgw['hus'] = relative_to_specific_humidity(
                            vars_pgw['hur'], pa_pgw, vars_pgw['ta'])


        # Determine current reference pressure (p_ref)
        if p_ref_inp is None:
            # get GCM pressure levels as candidates for reference pressure
            p_ref_opts = load_delta(delta_input_dir, 'zg',
                                laffile[TIME_ERA], era_step_dt).plev
            # surface pressure in ERA and PGW climate states
            p_sfc_era = pa_hl_era.sel(
                        {VERT_HL_ERA:np.max(pa_hl_era[VERT_HL_ERA])})
            p_sfc_pgw = pa_hl_pgw.sel(
                        {VERT_HL_ERA:np.max(pa_hl_pgw[VERT_HL_ERA])})
            # reference pressure from a former iteration already set?
            try:
                p_ref_last = p_ref
            except UnboundLocalError:
                p_ref_last = None
            # determine local reference pressure
            p_ref = xr.apply_ufunc(determine_p_ref, p_sfc_era, p_sfc_pgw, 
                    p_ref_opts, p_ref_last,
                    input_core_dims=[[],[],['plev'],[]],
                    vectorize=True)
            if VERT_HL_ERA in p_ref.coords:
                del p_ref[VERT_HL_ERA]
        else:
            p_ref = p_ref_inp

        # compute updated geopotential at reference pressure
        phi_ref_pgw = integ_geopot(pa_hl_pgw, laffile.FIS, vars_pgw['ta'], 
                                    vars_pgw['hus'], laffile[VERT_HL_ERA], p_ref)

        # recompute original geopotential
        phi_ref_era = integ_geopot(pa_hl_era, laffile.FIS,
                                    laffile[var_name_map['ta']], 
                                    laffile[var_name_map['hus']], 
                                    laffile[VERT_HL_ERA], p_ref)

        delta_phi_ref = phi_ref_pgw - phi_ref_era

        ## load climate delta for reference pressure level
        climate_delta_phi_ref = load_delta(delta_input_dir, 'zg',
                            laffile[TIME_ERA], era_step_dt) * CON_G
        climate_delta_phi_ref = climate_delta_phi_ref.sel(plev=p_ref)
        del climate_delta_phi_ref['plev']

        # error in future geopotential
        phi_ref_error = delta_phi_ref - climate_delta_phi_ref

        adj_ps = - adj_factor * ps_pgw / (
                CON_RD * 
                vars_pgw['ta'].sel(level=np.max(laffile.level))
            ) * phi_ref_error
        del adj_ps['level']

        phi_ref_max_error = np.abs(phi_ref_error).max().values
        if i_debug >= 2:
            print('### iteration {:03d}, phi max error: {}'.
                            format(it, phi_ref_max_error))

        it += 1

        if it > max_n_iter:
            raise ValueError('ERROR! Pressure adjustment did not converge '+
                  'for file {}.'.format(inp_era_file_path))


    ### TODO DEBUG TESTING
    #dps = ps_pgw-laffile.PS 
    #dps.to_netcdf(os.path.join(Path(out_era_file_path).parents[0], 
    #        'delta_ps_{}_{:%Y%m%d%H}.nc'.format(p_ref_inp, era_step_dt)))
    ### TODO DEBUG TESTING

    ## If re-interpolation is enabled, interpolate climate deltas for
    ## ua and va onto final PGW climate state ERA5 model levels.
    if i_reinterp:
        for var_name in ['ua', 'va']:
            if i_debug >= 2:
                print('add {}'.format(var_name))
            var_era = interp_logp_3d(laffile[var_name_map[var_name]], 
                            pa_era, pa_pgw, extrapolate='constant')
            dvar = load_delta_interp(delta_input_dir,
                    var_name, pa_pgw,
                    laffile[TIME_ERA], era_step_dt)
            vars_pgw[var_name] = var_era + dvar


    ## fields in ERA file
    laffile[var_name_map['ps']] = ps_pgw
    for var_name in ['ta', 'hus', 'ua', 'va']:
        laffile[var_name_map[var_name]] = vars_pgw[var_name]
    del laffile[var_name_map['hur']]


    ## save updated ERA5 file
    laffile.to_netcdf(out_era_file_path, mode='w')
    laffile.close()
    if i_debug >= 1:
        print('Done. Saved to file {}.'.format(out_era_file_path))








if __name__ == "__main__":
    ## input arguments
    parser = argparse.ArgumentParser(description =
                    'Perturb ERA5 with PGW climate deltas.')

    # first bc step to compute 
    parser.add_argument('-f', '--first_era_step', type=str, default='2006080200')

    # last bc step to compute 
    parser.add_argument('-l', '--last_era_step', type=str, default='2006080300')

    # delta hour increments
    parser.add_argument('-H', '--hour_inc_step', type=int, default=3)

    # input era5 directory
    parser.add_argument('-i', '--input_dir', type=str)

    # output era5 directory
    parser.add_argument('-o', '--output_dir', type=str)

    # climate delta directory (already remapped to ERA5 grid)
    parser.add_argument('-d', '--delta_input_dir', type=str)

    # ERA5 file name base
    parser.add_argument('-b', '--file_name_base', type=str, 
                        default='cas{:%Y%m%d%H}0000.nc')

    # number of parallel jobs
    parser.add_argument('-p', '--n_par', type=int, default=1)

    args = parser.parse_args()

    # first date and last date to datetime object
    first_era_step = datetime.strptime(args.first_era_step, '%Y%m%d%H')
    last_era_step = datetime.strptime(args.last_era_step, '%Y%m%d%H')

    # time steps to process
    era_step_dts = np.arange(first_era_step,
                        last_era_step+timedelta(hours=args.hour_inc_step),
                        timedelta(hours=args.hour_inc_step)).tolist()

    # if output directory doesn't exist create it
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    IMP = IterMP(njobs=args.n_par, run_async=True)
    fargs = {'delta_input_dir':args.delta_input_dir}
    step_args = []

    # iterate over time step and prepare function arguments
    for era_step_dt in era_step_dts:
        print(era_step_dt)

        # set output and input ERA5 file
        inp_era_file_path = os.path.join(args.input_dir, 
                args.file_name_base.format(era_step_dt))
        out_era_file_path = os.path.join(args.output_dir, 
                args.file_name_base.format(era_step_dt))

        step_args.append({
            'inp_era_file_path':inp_era_file_path,
            'out_era_file_path':out_era_file_path,
            'era_step_dt':era_step_dt}
        )

    # run in parallel if args.n_par > 1
    IMP.run(pgw_for_era5, fargs, step_args)

