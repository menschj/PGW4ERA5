#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
description     PGW for ERA5
author		Christoph Heim based on original developments by Roman Brogli
date created    12.01.2022
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
        load_delta_old,
        load_delta_interp,
        integ_geopot,
        interp_logp_3d,
        determine_p_ref,
        )
from constants import CON_G, CON_RD
from parallel import IterMP
from settings import *
##############################################################################

def test_delta(inp_era_file_path, out_era_file_path,
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
    inp_era_file_path =  '/net/o3/hymet_nobackup/heimc/data/pgw/ERA/ERA5/cas20060801000000.nc'
    new_out_path = os.path.join('/net/o3/hymet_nobackup/heimc/data/pgw/test/new/',
                                Path(out_era_file_path).name)
    old_out_path = os.path.join('/net/o3/hymet_nobackup/heimc/data/pgw/test/old/',
                                Path(out_era_file_path).name)
    diff_out_path = os.path.join('/net/o3/hymet_nobackup/heimc/data/pgw/test/diff/',
                                Path(out_era_file_path).name)
    old_delta_input_dir = '/net/o3/hymet_nobackup/heimc/data/pgw/regridded_old/Emon_RHint/MPI-ESM1-2-HR'

    laffile = xr.open_dataset(inp_era_file_path, decode_cf=False)

    var_name = 'tas'
    delta_tas_new = load_delta(delta_input_dir, var_name,
                            laffile[TIME_ERA], era_step_dt)
    delta_tas_old = load_delta_old(old_delta_input_dir, var_name,
                            laffile[TIME_ERA], era_step_dt)

    delta_tas_new.to_netcdf(new_out_path)
    delta_tas_old.to_netcdf(old_out_path)
    (delta_tas_new-delta_tas_old).to_netcdf(diff_out_path)
    #quit()


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
    IMP.run(test_delta, fargs, step_args)

