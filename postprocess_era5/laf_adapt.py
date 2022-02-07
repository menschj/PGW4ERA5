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
        integ_geopot_era5,
        integ_geopot_downward_era5,
        integ_pressure_era5,
        integ_pressure_upward_era5,
        load_delta,
        debug_interp,
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
    'T_S':      'tas',
    'RELHUM':   'hur',
    'QV':       'hus',
    'U':        'ua',
    'V':        'va',
    'PS':       'ps',
}

#var_name_map = {
#    'PHI':      'PHI',
#    'T':        'T',
#    'T_S':      'T_S',
#    'RELHUM':   'RELHUM',
#    'QV':       'QV',
#    'U':        'U',
#    'V':        'V',
#    'PS':       'PS',
#}



def lafadapt(inp_laf_path, out_laf_path, delta_inp_path,
            delta_time_step, new_time_string, laf_dt,
            RELHUM_or_QV_delta):

    timer = Timer('seconds')
    timer.start('all')

    laffile = xr.open_dataset(inp_laf_path, decode_cf=False)
    print('load done')

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
    print('rh computed')

    p_ref_inp = 50000
    p_ref_inp = None
    adj_factor = 0.95
    thresh_phi_ref_max_error = 0.10
    #max_n_iter = 6
    max_n_iter = 8
    extrapolate = 'linear'
    extrapolate = 'constant'
    #extrapolate = 'off'


    ## determine future climate state surface pressure using iterative
    ## procedure
    delta_PS = xr.zeros_like(laffile.PS)
    # TODO: make usable without PS climate delta
    #adj_PS = xr.zeros_like(laffile.PS)
     
    delta_PS_gcm = get_delta_era5(laffile.PS,
                    var_name_map['PS'], laffile, 
                    delta_inp_path, delta_time_step, laf_dt)
    adj_PS = delta_PS_gcm

    phi_ref_max_error = np.inf

    ### interpolate climate deltas onto ERA5 grid
    deltas = {}
    for var_name in ['T',RELHUM_or_QV_delta,'U','V']:
        delta_var = get_delta_interp_era5(
                laffile[var_name],
                P_era, var_name_map[var_name], laffile, 
                delta_inp_path, delta_time_step, laf_dt)
        deltas[var_name] = delta_var

    ##### TODO DEBUG start
    #var_name = 'T'
    #dT = get_delta_interp_era5(
    #        laffile[var_name],
    #        P_era, var_name_map[var_name], laffile, 
    #        delta_inp_path, delta_time_step, laf_dt)
    #T_pgw = laffile.T + dT

    #var_name = 'QV'
    #dQV = get_delta_interp_era5(
    #        laffile[var_name],
    #        P_era, var_name_map[var_name], laffile, 
    #        delta_inp_path, delta_time_step, laf_dt)
    #QV_pgw = laffile.QV + dQV
    ##### TODO DEBUG stop

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

        # interpolate variables onto new model levels
        #if np.sum(adj_PS.values != 0):
        var_name = 'T'
        #T_era = interp_nonfixed(laffile.T, 
        #                P_era, P_pgw, 'level', 'level',
        #                extrapolate=extrapolate)
        T_era = laffile[var_name]
        #dT = get_delta_interp_era5(
        #        laffile[var_name],
        #        P_pgw, var_name_map[var_name], laffile, 
        #        delta_inp_path, delta_time_step, laf_dt)
        dT = deltas[var_name]
        T_pgw = T_era + dT
        (T_era - laffile[var_name]).to_netcdf('interp_{}.nc'.format(var_name))
        dT.to_netcdf('delta_{}.nc'.format(var_name))
        #quit()

        if RELHUM_or_QV_delta == 'RELHUM':
            var_name = 'RELHUM'
            #RELHUM_era = interp_nonfixed(laffile.RELHUM, 
            #                P_era, P_pgw, 'level', 'level',
            #                extrapolate=extrapolate)
            RELHUM_era = laffile[var_name]
            #dRELHUM = get_delta_interp_era5(
            #        laffile[var_name],
            #        P_pgw, var_name_map[var_name], laffile, 
            #        delta_inp_path, delta_time_step, laf_dt)
            dRELHUM = deltas[var_name]
            RELHUM_pgw = RELHUM_era + dRELHUM

            # convert relative humidity to speicifc humidity in pgw
            QV_pgw = relative_to_specific_humidity(
                        RELHUM_pgw, P_pgw, T_pgw)

        elif RELHUM_or_QV_delta == 'QV':
            var_name = 'QV'
            #QV_era = interp_nonfixed(laffile.QV, 
            #                P_era, P_pgw, 'level', 'level',
            #                extrapolate=extrapolate)
            QV_era = laffile[var_name]
            #dQV = get_delta_interp_era5(
            #        laffile[var_name],
            #        P_pgw, var_name_map[var_name], laffile, 
            #        delta_inp_path, delta_time_step, laf_dt)
            dQV = deltas[var_name]
            QV_pgw = QV_era + dQV
            (QV_era - laffile[var_name]).to_netcdf('interp_{}.nc'.format(var_name))
            dQV.to_netcdf('delta_{}.nc'.format(var_name))

        else: 
            raise NotImplementedError()

        ## TODO DEBUG start
        #T_pgw = laffile.T
        #QV_pgw = laffile.QV
        ## TODO DEBUG stop

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
            #p_ref.to_netcdf('test.nc')
            #quit()
        else:
            p_ref = p_ref_inp

        # compute updated geopotential at reference pressure
        PHI_hl_pgw, phi_ref_pgw = integ_geopot_era5(P_hl_pgw,
                laffile.FIS, T_pgw, QV_pgw, laffile.level1, p_ref)
        #print('{} phi_ref_pgw'.format(phi_ref_pgw.mean().values/CON_G))


        # recompute original geopotential
        PHI_hl_era, phi_ref_era = integ_geopot_era5(P_hl_era,
                laffile.FIS, laffile.T, laffile.QV, laffile.level1, p_ref)

        delta_phi_ref = phi_ref_pgw - phi_ref_era
        #print('{} delta_phi_ref'.format(delta_phi_ref.mean().values/CON_G))

        ## load climate delta for reference pressure level
        climate_delta_phi_ref = load_delta(delta_inp_path, var_name_map['PHI'],
                            laf_dt, laffile.time, delta_time_step) * CON_G
        climate_delta_phi_ref = climate_delta_phi_ref.assign_coords(
                                    lat=laffile.lat.values)
        climate_delta_phi_ref = climate_delta_phi_ref.sel(plev=p_ref)
        del climate_delta_phi_ref['plev']
        #print('{} gcm delta_phi_ref'.format(climate_delta_phi_ref.mean().values/CON_G))

        # error in future geopotential
        phi_ref_error = delta_phi_ref - climate_delta_phi_ref
        #phi_ref_error.to_netcdf('test_{}.nc'.format(it))
        #quit()
        #print('{} dphi deviation'.format(phi_ref_error.mean().values/CON_G))

        adj_PS = - adj_factor * PS_pgw / (CON_RD * 
                T_pgw.sel(level=np.max(laffile.level))) * phi_ref_error
        del adj_PS['level']
        #print('{} adj_PS'.format(adj_PS.mean().values))

        phi_ref_max_error = np.abs(phi_ref_error).max().values
        print('### iteration {:03d}, phi max error: {}'.
                            format(it, phi_ref_max_error))

        it += 1

        if it > max_n_iter:
            print('DID NOT CONVERGE!')
            break


    ## TODO DEBUG start
    #PHI_hl_pgw, phi_ref_pgw = integ_geopot_era5(P_hl_pgw,
    #        laffile.FIS, T_pgw, QV_pgw, laffile.level1, 50000)
    #phi_ref_pgw.to_netcdf('test.nc')
    #quit()
    dPS = PS_pgw-laffile.PS 
    dPS.to_netcdf('delta_ps_{}_{:%Y%m%d%H}.nc'.format(p_ref_inp, laf_dt))
    dPS.values -= xr.open_dataset('climate_delta_ps.nc').ps.values
    print('mean error PS: {}'.format(dPS.mean().values))
    dPS.to_netcdf('error_delta_ps_{}_{:%Y%m%d%H}.nc'.format(p_ref_inp, laf_dt))
    ## TODO DEBUG stop

    ## replace T and QV in ERA file
    laffile['PS'] = PS_pgw
    laffile['T'] = T_pgw
    laffile['QV'] = QV_pgw
    del laffile['RELHUM']

    P_hl_era.to_netcdf('P_hl_era.nc')
    P_hl_pgw.to_netcdf('P_hl_pgw.nc')

    ## interpolate U and V onto final pgw grid and replace in ERA file
    for var_name in ['U', 'V']:
        print('add {}'.format(var_name))
        #var_era = interp_nonfixed(laffile[var_name], 
        #                P_era, P_pgw, 'level', 'level',
        #                extrapolate=extrapolate)
        var_era = laffile[var_name]
        #dvar = get_delta_interp_era5(
        #        laffile[var_name],
        #        P_pgw, var_name_map[var_name], laffile, 
        #        delta_inp_path, delta_time_step, laf_dt)
        dvar = deltas[var_name]
        #quit()
        var_pgw = var_era + dvar
        (var_era - laffile[var_name]).to_netcdf('interp_{}.nc'.format(var_name))
        dvar.to_netcdf('delta_{}.nc'.format(var_name))
        laffile[var_name] = var_pgw


    # also update T_SKIN 
    var_name = 'T_S'
    print('add {}'.format(var_name))
    delta_T_SKIN = get_delta_era5(laffile.T_SKIN,
                    var_name_map['T_S'], laffile, 
                    delta_inp_path, delta_time_step, laf_dt)
    laffile.T_SKIN.values += delta_T_SKIN.values



    laffile.to_netcdf(out_laf_path, mode='w')
    laffile.close()
    print('saved to file {}.'.format(out_laf_path))
    timer.stop('all')
    timer.print_report()








if __name__ == "__main__":

    ## input arguments
    parser = argparse.ArgumentParser(description =
                    'Perturb COSMO initial condition with PGW climate deltas.')
    # delta hour increments
    parser.add_argument('-d', '--delta_hour_inc', type=int, default=3)
    # sim start date 
    parser.add_argument('-s', '--sim_start_date', type=str, default='2006080100')
    # first bc step to compute 
    parser.add_argument('-f', '--first_bc_step', type=str, default='2006080200')
    # last bc step to compute 
    parser.add_argument('-l', '--last_bc_step', type=str, default='2006080300')
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
    delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded_era5'
    delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded_era5_Emon'

    delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded_delta_era5'
    delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded_delta_era5_test3'
    delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded_delta_era5_Emon'
    #delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded_delta_era5_Emon_test3'

    RELHUM_or_QV_delta = 'RELHUM'
    RELHUM_or_QV_delta = 'QV'

    pgw_sim_name_ending = 'pgw9'


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
    fargs = {'delta_inp_path':delta_inp_path,
             'RELHUM_or_QV_delta':RELHUM_or_QV_delta}
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

        #lafadapt(inp_laf_path, out_laf_path, delta_inp_path,
        #        delta_time_step, new_time_string, laf_dt)

        step_args.append({
            'inp_laf_path':inp_laf_path,
            'out_laf_path':out_laf_path,
            'delta_time_step':delta_time_step,
            'new_time_string':new_time_string,
            'laf_dt':laf_dt}
        )



    IMP.run(lafadapt, fargs, step_args)

