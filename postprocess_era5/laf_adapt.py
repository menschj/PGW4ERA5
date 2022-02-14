import xarray as xr
import sys, argparse
from pathlib import Path
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.interpolate import RegularGridInterpolator
from base.functions import (
        hour_of_year,
        specific_to_relative_humidity,
        relative_to_specific_humidity,
        get_delta_era5,
        get_delta_interp_era5,
        integ_geopot,
        load_delta,
        interp_nonfixed,
        interp_vprof,
        determine_p_ref,
        )
from constants import CON_G, CON_RD
from package.utilities import Timer
from package.mp import IterMP


var_name_map = {
    'PHI':      'zg',
    'T':        'ta',
    'T_SKIN':   'tas',
    'RELHUM':   'hur',
    'QV':       'hus',
    'U':        'ua',
    'V':        'va',
    'PS':       'ps',
}



def lafadapt(inp_laf_path, out_laf_path, delta_inp_path,
            delta_time_step, new_time_string, laf_dt):

    timer = Timer('seconds')
    timer.start('all')

    laffile = xr.open_dataset(inp_laf_path, decode_cf=False)
    laffile.time.attrs['units'] = new_time_string

    # compute pressure on era5 levels
    P_era = (laffile.akm + laffile.PS * laffile.bkm).transpose(
                'time','level','lat','lon')
    P_hl_era = (laffile.ak + laffile.PS * laffile.bk).transpose(
                'time','level1','lat','lon')


    # compute relative humidity in ERA climate state
    laffile['RELHUM'] = specific_to_relative_humidity(
                        laffile.QV, P_era, laffile.T).transpose(
                        'time', 'level', 'lat', 'lon')

    # also update T_SKIN 
    var_name = 'T_SKIN'
    print('update {}'.format(var_name))
    delta_T_SKIN = get_delta_era5(laffile.T_SKIN,
                    var_name_map[var_name], laffile, 
                    delta_inp_path, delta_time_step, laf_dt)
    laffile.T_SKIN.values += delta_T_SKIN.values

    # update T_SO
    var_name = 'T_SO'
    print('update {}'.format(var_name))
    delta_T_SO_clim = get_delta_era5(laffile.T_SKIN,
                    var_name_map['T_SKIN'], 
                    laffile, delta_inp_path,
                    delta_time_step, date_time=None).mean(dim=['time'])
    delta_T_SO = (
            delta_T_SO_clim + np.exp(-laffile.soil1/2.8) * 
                    (delta_T_SKIN - delta_T_SO_clim)
    )
    delta_T_SO = delta_T_SO.transpose('time', 'soil1', 'lat', 'lon')
    laffile.T_SO.values += delta_T_SO


    i_reinterp = 0
    #p_ref_inp = 50000
    p_ref_inp = None
    adj_factor = 0.95
    thresh_phi_ref_max_error = 0.10
    max_n_iter = 10
    extrapolate = 'constant'
    #extrapolate = 'linear'
    #extrapolate = 'off'


    ## determine future climate state surface pressure using iterative
    ## procedure
    delta_PS = xr.zeros_like(laffile.PS)
    adj_PS = xr.zeros_like(laffile.PS)
    delta_PS_gcm = get_delta_era5(laffile.PS,
                    var_name_map['PS'], laffile, 
                    delta_inp_path, delta_time_step, laf_dt)
    #adj_PS = delta_PS_gcm

    phi_ref_max_error = np.inf

    if not i_reinterp:
        ### interpolate climate deltas onto ERA5 grid
        deltas = {}
        for var_name in ['T','RELHUM','U','V']:
            print('update {}'.format(var_name))
            delta_var = get_delta_interp_era5(
                    laffile[var_name],
                    P_era, var_name_map[var_name], laffile, 
                    delta_inp_path, delta_time_step, laf_dt)
            deltas[var_name] = delta_var

        var_name = 'T'
        T_pgw = laffile[var_name] + deltas[var_name]
        var_name = 'U'
        U_pgw = laffile[var_name] + deltas[var_name]
        var_name = 'V'
        V_pgw = laffile[var_name] + deltas[var_name]

        var_name = 'RELHUM'
        RELHUM_pgw = laffile[var_name] + deltas[var_name]
        # convert relative humidity to speicifc humidity in pgw
        QV_pgw = relative_to_specific_humidity(
                    RELHUM_pgw, P_era, T_pgw)


    it = 1
    while phi_ref_max_error > thresh_phi_ref_max_error:

        # update surface pressure
        delta_PS += adj_PS
        #print('{} delta_PS'.format(delta_PS.mean().values))
        PS_pgw = laffile.PS + delta_PS

        # recompute pressure on half levels
        P_pgw = (laffile.akm + PS_pgw * laffile.bkm).transpose(
                    'time','level','lat','lon')
        P_hl_pgw = (laffile.ak + PS_pgw * laffile.bk).transpose(
                    'time','level1','lat','lon')

        if i_reinterp:
            # interpolate variables onto new model levels
            var_name = 'T'
            T_era = interp_nonfixed(laffile.T, 
                            P_era, P_pgw, 'level', 'level',
                            extrapolate=extrapolate)
            dT = get_delta_interp_era5(
                    laffile[var_name],
                    P_pgw, var_name_map[var_name], laffile, 
                    delta_inp_path, delta_time_step, laf_dt)
            T_pgw = T_era + dT

            var_name = 'RELHUM'
            RELHUM_era = interp_nonfixed(laffile.RELHUM, 
                            P_era, P_pgw, 'level', 'level',
                            extrapolate=extrapolate)
            dRELHUM = get_delta_interp_era5(
                    laffile[var_name],
                    P_pgw, var_name_map[var_name], laffile, 
                    delta_inp_path, delta_time_step, laf_dt)
            RELHUM_pgw = RELHUM_era + dRELHUM

            # convert relative humidity to speicifc humidity in pgw
            QV_pgw = relative_to_specific_humidity(
                        RELHUM_pgw, P_pgw, T_pgw)


        # get p_ref
        if p_ref_inp is None:
            p_ref_opts = load_delta(delta_inp_path, var_name_map['PHI'],
                                laf_dt, laffile.time, delta_time_step).plev
            p_sfc_era = P_hl_era.sel(level1=np.max(P_hl_era.level1))
            p_sfc_pgw = P_hl_pgw.sel(level1=np.max(P_hl_pgw.level1))
            try:
                p_ref_last = p_ref
            except UnboundLocalError:
                p_ref_last = None
            p_ref = xr.apply_ufunc(determine_p_ref, p_sfc_era, p_sfc_pgw, 
                    p_ref_last, p_ref_opts,
                    input_core_dims=[[],[],[],['plev']],
                    vectorize=True)
            if 'level1' in p_ref.coords:
                del p_ref['level1']
        else:
            p_ref = p_ref_inp

        # compute updated geopotential at reference pressure
        PHI_hl_pgw, phi_ref_pgw = integ_geopot(P_hl_pgw,
                laffile.FIS, T_pgw, QV_pgw, laffile.level1, p_ref)

        # recompute original geopotential
        PHI_hl_era, phi_ref_era = integ_geopot(P_hl_era,
                laffile.FIS, laffile.T, laffile.QV, laffile.level1, p_ref)

        delta_phi_ref = phi_ref_pgw - phi_ref_era

        ## load climate delta for reference pressure level
        climate_delta_phi_ref = load_delta(delta_inp_path, var_name_map['PHI'],
                            laf_dt, laffile.time, delta_time_step) * CON_G
        climate_delta_phi_ref = climate_delta_phi_ref.assign_coords(
                                    lat=laffile.lat.values)
        climate_delta_phi_ref = climate_delta_phi_ref.sel(plev=p_ref)
        del climate_delta_phi_ref['plev']

        # error in future geopotential
        phi_ref_error = delta_phi_ref - climate_delta_phi_ref

        adj_PS = - adj_factor * PS_pgw / (CON_RD * 
                T_pgw.sel(level=np.max(laffile.level))) * phi_ref_error
        del adj_PS['level']

        phi_ref_max_error = np.abs(phi_ref_error).max().values
        print('### iteration {:03d}, phi max error: {}'.
                            format(it, phi_ref_max_error))

        it += 1

        if it > max_n_iter:
            raise ValueError('ERROR! Pressure adjustment did not converge '+
                  'for file {}.'.format(inp_laf_path))


    ## TODO DEBUG start
    dPS = PS_pgw-laffile.PS 
    dPS.to_netcdf(os.path.join(Path(out_laf_path).parents[0], 
            'delta_ps_{}_{:%Y%m%d%H}.nc'.format(p_ref_inp, laf_dt)))

    #dPS.to_netcdf('delta_ps_{}_{:%Y%m%d%H}.nc'.format(p_ref_inp, laf_dt))
    #dPS.values -= xr.open_dataset('climate_delta_ps.nc').ps.values
    #print('mean error PS: {}'.format(dPS.mean().values))
    #dPS.to_netcdf('error_delta_ps_{}_{:%Y%m%d%H}.nc'.format(p_ref_inp, laf_dt))
    ## TODO DEBUG stop

    if i_reinterp:
        ## interpolate U and V onto final pgw grid and replace in ERA file
        for var_name in ['U', 'V']:
            print('add {}'.format(var_name))
            var_era = interp_nonfixed(laffile[var_name], 
                            P_era, P_pgw, 'level', 'level',
                            extrapolate=extrapolate)
            dvar = get_delta_interp_era5(
                    laffile[var_name],
                    P_pgw, var_name_map[var_name], laffile, 
                    delta_inp_path, delta_time_step, laf_dt)
            if var_name == 'U':
                U_pgw = var_era + dvar
            elif var_name == 'V':
                V_pgw = var_era + dvar


    ## fields in ERA file
    laffile['PS'] = PS_pgw
    laffile['T'] = T_pgw
    laffile['QV'] = QV_pgw
    laffile['U'] = U_pgw
    laffile['V'] = V_pgw
    del laffile['RELHUM']


    laffile.to_netcdf(out_laf_path, mode='w')
    laffile.close()
    print('saved to file {}.'.format(out_laf_path))
    timer.stop('all')
    timer.print_report()








if __name__ == "__main__":
    ## input arguments
    parser = argparse.ArgumentParser(description =
                    'Perturb ERA5 with PGW climate deltas.')
    # delta hour increments
    parser.add_argument('-d', '--delta_hour_inc', type=int, default=3)
    # first bc step to compute 
    parser.add_argument('-f', '--first_bc_step', type=str, default='2006080200')
    # last bc step to compute 
    parser.add_argument('-l', '--last_bc_step', type=str, default='2006080300')
    # sim start date 
    parser.add_argument('-s', '--sim_start_date', type=str, default='2006080100')
    # sim name excluding ctrl/pgw
    parser.add_argument('-n', '--sim_name_base', type=str, default='SA_3')
    # number of parallel jobs
    parser.add_argument('-p', '--n_par', type=int, default=1)
    args = parser.parse_args()
    print(args)

    sim_start_date = datetime.strptime(args.sim_start_date, '%Y%m%d%H')
    first_bc_step = datetime.strptime(args.first_bc_step, '%Y%m%d%H')
    last_bc_step = datetime.strptime(args.last_bc_step, '%Y%m%d%H')
    sim_name_base = args.sim_name_base
    wd_path = '/scratch/snx3000/heimc/lmp/wd'
    changeyears = 0
    #delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded_old/Emon_RHint/MPI-ESM1-2-HR'
    #delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded_old/Emon/MPI-ESM1-2-HR'
    #delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded_old/Amon/MPI-ESM1-2-HR'

    delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded/Emon_xr/MPI-ESM1-2-HR'
    #delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded/Amon/MPI-ESM1-2-HR'

    pgw_sim_name_ending = 'pgw'
    pgw_sim_name_ending = 'pgwtest'

    pgw_sim_start_date = sim_start_date + relativedelta(years=changeyears)
    laf_dts = np.arange(first_bc_step,
                        last_bc_step+timedelta(hours=args.delta_hour_inc),
                        timedelta(hours=args.delta_hour_inc)).tolist()
    out_laf_dir = os.path.join(wd_path, 
                    '{:%y%m%d}00_{}_{}/int2lm_in/'.format(
                        sim_start_date, sim_name_base, pgw_sim_name_ending))

    #if output directory doesn't exist create it
    Path(out_laf_dir).mkdir(parents=True, exist_ok=True)

    IMP = IterMP(njobs=args.n_par, run_async=True)
    fargs = {'delta_inp_path':delta_inp_path}
    step_args = []

    for laf_dt in laf_dts:
        print(laf_dt)

        new_laf_dt = laf_dt + relativedelta(years=changeyears)
        new_time_string = 'seconds since {:%Y-%m-%d %H:%M:%S}'.format(new_laf_dt)

        inp_laf_path = os.path.join(wd_path, 
                '{:%y%m%d}00_{}_ctrl/int2lm_in/cas{:%Y%m%d%H}0000.nc'.format(
                            sim_start_date, sim_name_base, laf_dt))

        delta_time_step = int(hour_of_year(laf_dt)/args.delta_hour_inc)
        print('use time step {}'.format(delta_time_step))

        out_laf_path = os.path.join(out_laf_dir, 
                'cas{:%Y%m%d%H}0000.nc'.format(new_laf_dt))
        print(out_laf_path)

        step_args.append({
            'inp_laf_path':inp_laf_path,
            'out_laf_path':out_laf_path,
            'delta_time_step':delta_time_step,
            'new_time_string':new_time_string,
            'laf_dt':laf_dt}
        )

    IMP.run(lafadapt, fargs, step_args)

