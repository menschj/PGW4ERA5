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
        )
from constants import CON_G, CON_RD
from package.utilities import Timer


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
            delta_time_step, new_time_string, laf_dt):

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

    # update T_SKIN 
    var_name = 'T_S'
    print('add {}'.format(var_name))
    delta_T_SKIN = get_delta_era5(laffile.T_SKIN,
                    var_name_map['T_S'], laffile, 
                    delta_inp_path, delta_time_step, laf_dt)
    laffile.T_SKIN.values += delta_T_SKIN.values


    #handle, = plt.plot(laffile.POTT.mean(dim=['time','lon','lat']),
    #            P_era.mean(dim=['time','lon','lat']),
    #        label='direct')
    #plt.xlim(290,310)
    #plt.ylim(103000,80000)
    #plt.show()
    #quit()


    ## integrate geopotential upward
    #p_ref = 30000
    ##PHI, PHI_hl = integ_geopot_era5(laffile)
    #PHI_hl, phi_ref, p_ref_star, hl_ref_star = integ_geopot_era5(laffile, p_ref)
    #print('geopotential computed')

    ## add climate change delta to geopotential at reference pressure level
    #phi_ref_old = phi_ref.copy()
    #delta_phi_ref = load_delta(delta_inp_path, var_name_map['PHI'],
    #                    laf_dt, delta_time_step).sel(plev=p_ref)
    #phi_ref = phi_ref + delta_phi_ref.values
    ##phi_ref.to_netcdf('test.nc')
    ##print(phi_ref)


    ## integrate geopotential downward
    ##integ_geopot_downward_era5(laffile, p_ref, p_ref_star, phi_ref, hl_ref_star)
    #quit()

    ############################## NEW VERSION START
    p_refs = [100000, 92500, 85000, 70000, 50000, 30000, 10000]
    #p_refs = [50000, 40000, 30000, 20000, 10000]
    #p_refs = [85000, 50000, 10000]
    p_refs = [50000]
    adj_factor = 0.95
    thresh_phi_ref_max_error = 0.10
    i_plot_type = 0
    max_n_iter = 10
    #var_names = ['T', 'U', 'V', 'RELHUM']
    #var_names = ['T', 'RELHUM']




    out_vars = {}
    for p_ref in p_refs:
        print('#################################################################')
        print('#################################################################')
        print('p_ref {}'.format(p_ref))
        print('#################################################################')
        out_vars[p_ref] = {}

        # compute original geopotential
        PHI_hl_era, phi_ref_era, phi_ref_star_era = integ_geopot_era5(P_hl_era,
                laffile.FIS, laffile.T, laffile.QV, laffile.level1, p_ref)
        out_vars[p_ref]['PHI_hl_era'] = PHI_hl_era.copy() / CON_G
        out_vars[p_ref]['phi_ref_era'] = phi_ref_era.copy() / CON_G
        out_vars[p_ref]['P_hl_era'] = P_hl_era.copy().transpose('time','level1','lat','lon')

        # load climate delta for reference pressure level
        climate_delta_phi_ref = load_delta(delta_inp_path, var_name_map['PHI'],
                            laf_dt, laffile.time, delta_time_step).sel(plev=p_ref) * CON_G
        climate_delta_phi_ref = climate_delta_phi_ref.assign_coords(lat=laffile.lat.values)

        ## determine future climate state surface pressure using iterative
        ## procedure
        delta_PS = xr.zeros_like(laffile.PS)
        #adj_PS = xr.zeros_like(laffile.PS)
         
        delta_PS_GCM = get_delta_era5(laffile.PS,
                        var_name_map['PS'], laffile, 
                        delta_inp_path, delta_time_step, laf_dt)
        adj_PS = delta_PS_GCM
        ## TODO debug start
        #climate_delta_phi_ref.values[:] = 0 * CON_G
        ## TODO debug stop

        phi_ref_max_error = np.inf

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
            print('{} delta_PS'.format(delta_PS.mean().values))
            PS_pgw = laffile.PS + delta_PS

            # recompute pressure on half levels
            P_pgw = (laffile.akm + PS_pgw * laffile.bkm).transpose(
                        'time','level','lat','lon')
            P_hl_pgw = (laffile.ak + PS_pgw * laffile.bk).transpose(
                        'time','level1','lat','lon')

            # interpolate variables onto new model levels
            #if np.sum(adj_PS.values != 0):
            var_name = 'T'
            T_era = interp_nonfixed(laffile.T, 
                            P_era, P_pgw, 'level', 'level')
            dT = get_delta_interp_era5(
                    laffile[var_name],
                    P_pgw, var_name_map[var_name], laffile, 
                    delta_inp_path, delta_time_step, laf_dt)
            T_pgw = T_era + dT

            #var_name = 'RELHUM'
            #RELHUM_era = interp_nonfixed(laffile.RELHUM, 
            #                P_era, P_pgw, 'level', 'level')
            #dRELHUM = get_delta_interp_era5(
            #        laffile[var_name],
            #        P_pgw, var_name_map[var_name], laffile, 
            #        delta_inp_path, delta_time_step, laf_dt)
            #RELHUM_pgw = RELHUM_era + dRELHUM

            ## convert relative humidity to speicifc humidity in pgw
            #QV_pgw = relative_to_specific_humidity(
            #            RELHUM_pgw, P_pgw, T_pgw)

            var_name = 'QV'
            QV_era = interp_nonfixed(laffile.QV, 
                            P_era, P_pgw, 'level', 'level')
            dQV = get_delta_interp_era5(
                    laffile[var_name],
                    P_pgw, var_name_map[var_name], laffile, 
                    delta_inp_path, delta_time_step, laf_dt)
            QV_pgw = QV_era + dQV
            


            ## TODO DEBUG start
            #T_pgw = laffile.T
            #QV_pgw = laffile.QV
            ## TODO DEBUG stop

            # compute updated geopotential at reference pressure
            PHI_hl_pgw, phi_ref_pgw, phi_ref_star_pgw = integ_geopot_era5(P_hl_pgw,
                    laffile.FIS, T_pgw, QV_pgw, laffile.level1, p_ref)
            #print('{} phi_ref_pgw'.format(phi_ref_pgw.mean().values/CON_G))
            #print('{} phi_ref_star'.format(phi_ref_star_pgw.mean().values/CON_G))

            delta_phi_ref = phi_ref_pgw - phi_ref_era
            #print('{} delta_phi_ref'.format(delta_phi_ref.mean().values/CON_G))
            #print('{} GCM delta_phi_ref'.format(climate_delta_phi_ref.mean().values/CON_G))

            # error in future geopotential
            phi_ref_error = delta_phi_ref - climate_delta_phi_ref
            #print('{} dphi deviation'.format(phi_ref_error.mean().values/CON_G))

            adj_PS = - adj_factor * PS_pgw / (CON_RD * 
                    T_pgw.sel(level=np.max(laffile.level))) * phi_ref_error
            del adj_PS['level']
            #print('{} adj_PS'.format(adj_PS.mean().values))

            #phi_ref_max_error = np.sqrt(np.square(phi_ref_error).mean()).values
            phi_ref_max_error = phi_ref_error.max().values
            print('####### iteration {:03d}, phi max error: {}'.format(it, phi_ref_max_error))

            it += 1

            if it > max_n_iter:
                print('DID NOT CONVERGE!')
                break


        out_vars[p_ref]['PS_pgw'] = laffile.PS.copy() + delta_PS
        out_vars[p_ref]['P_hl_pgw'] = P_hl_pgw.copy()
        out_vars[p_ref]['PHI_hl_pgw'] = PHI_hl_pgw.copy() / CON_G
        out_vars[p_ref]['phi_ref_pgw'] = phi_ref_pgw.copy() / CON_G

        ## TODO DEBUG start
        dPS = out_vars[p_ref]['PS_pgw']-laffile.PS 
        dPS.to_netcdf('delta_ps_{}_{:%Y%m%d%H}.nc'.format(p_ref, laf_dt))
        dPS.values -= xr.open_dataset('delta_ps.nc').ps.values
        print('mean error PS: {}'.format(dPS.mean().values))
        dPS.to_netcdf('error_delta_ps_{}_{:%Y%m%d%H}.nc'.format(p_ref, laf_dt))
        ## TODO DEBUG stop



    handles = []
    for p_ref in p_refs:
        print(p_ref)

        # interpolate
        PHI_hl_era = debug_interp(out_vars[p_ref]['PHI_hl_era'], 
                                  out_vars[p_ref]['P_hl_era'])
        PHI_hl_pgw = debug_interp(out_vars[p_ref]['PHI_hl_pgw'],
                                  out_vars[p_ref]['P_hl_pgw'])
        dPHI_hl = PHI_hl_pgw - PHI_hl_era
        #print(dPHI_hl)

        phi_ref_era = debug_interp(out_vars[p_ref]['phi_ref_era'])
        phi_ref_pgw = debug_interp(out_vars[p_ref]['phi_ref_pgw'])
        dphi_ref = phi_ref_pgw - phi_ref_era

        PS_pgw = debug_interp(out_vars[p_ref]['PS_pgw'])
        PS_era = debug_interp(laffile.PS)
        dPS = PS_pgw - PS_era

        if i_plot_type == 1:
            plt.scatter(dphi_ref, p_ref)
            handle, = plt.plot(dPHI_hl, dPHI_hl.plev, label=p_ref)
            handles.append(handle)

            plt.xlabel('dphi [m]')

        elif i_plot_type == 2:
            dPS = dPS.where(dPS != 0, 0.001)
            plt.scatter((dphi_ref/dPS), p_ref)
            handle, = plt.plot((dPHI_hl/dPS), dPHI_hl.plev, label=p_ref)
            handles.append(handle)
            plt.xlabel('dphi/dps [m/Pa]')

        elif i_plot_type == 3:
            handle = plt.scatter(dPS, p_ref,
                                label=p_ref)
            handles.append(handle)
            plt.xlabel('dps [Pa]')

    # add GCM
    dPHI_gcm = load_delta(delta_inp_path, var_name_map['PHI'],
                    laf_dt, laffile.time, delta_time_step).sel(plev=slice(100000,10000))
    dPHI_gcm = debug_interp(dPHI_gcm)
    #print(dPHI_gcm)
    dPS_gcm = load_delta(delta_inp_path, var_name_map['PS'],
                        laf_dt, laffile.time, delta_time_step)
    dPS_gcm = debug_interp(dPS_gcm)

    if i_plot_type == 1:
        #pass
        handle, = plt.plot(dPHI_gcm, dPHI_gcm.plev, label='GCM', color='k')
        handles.append(handle)

    elif i_plot_type == 2:
        pass
        #handle, = plt.plot((dPHI_gcm/dPS_gcm),
        #                    dPHI_gcm.plev, label='GCM', color='k')
        #handles.append(handle)

    elif i_plot_type == 3:
        handle = plt.axvline(x=dPS_gcm,
                            label='GCM', color='k')
        handles.append(handle)

    if i_plot_type > 0:
        plt.ylim(103000,10000)
        #plt.ylim(103000,80000)
        #plt.xlim(0,50)
        plt.ylabel('p [Pa]')
        plt.legend(handles=handles)
        plt.show()
        quit()
    ############################## NEW VERSION STOP



    del laffile['RELHUM']

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
    args = parser.parse_args()
    print(args)

    sim_start_date = datetime.strptime(args.sim_start_date, '%Y%m%d%H')
    first_bc_step = datetime.strptime(args.first_bc_step, '%Y%m%d%H')
    last_bc_step = datetime.strptime(args.last_bc_step, '%Y%m%d%H')
    sim_name_base = args.sim_name_base
    wd_path = '/scratch/snx3000/heimc/lmp/wd'
    changeyears = 0
    delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded_era5'
    delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded_delta_era5'
    delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded_delta_era5_test'
    delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded_delta_era5_test2'
    delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded_delta_era5_test3'
    delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded_delta_era5_Emon'
    #delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded_delta_era5_Emon_test3'

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

        lafadapt(inp_laf_path, out_laf_path, delta_inp_path,
                delta_time_step, new_time_string, laf_dt)
