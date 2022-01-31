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
        add_delta_era5,
        add_delta_interp_era5,
        integ_geopot_era5,
        integ_geopot_downward_era5,
        integ_pressure_era5,
        integ_pressure_upward_era5,
        load_delta,
        )
from constants import CON_G, CON_RD
"""
Can be used to add the calculated difference in all necessary variables to the initial condition (laf) file.
This very fast, just run it in the console.
It requires the difference in relative humidity to adapt the specific humidity!

Input:
	inp_laf_path: Path to the original laf-file from the "base" simulation
	(e.g. reanalysis driven or historical simulation). The name of the laf file
	must be as outputted by int2lm (e.g. laf1970010100.nc ).

	newyear: What year to use in the files (change it to the future to adapt CO2 levels)

	delta_time_step: Which timestep within the annual cycle is apropriate to adapt
	the laf file? (0 for beginning of year; otherwise dayofyear*timesteps_per_day)

	new_time_string: What timestamp should be used for the adapted laf file?
	Put the exact time of the new laf file in the format 'seconds since yyyy-mm-dd hh:mm:ss'

	out_laf_dir: In which folder should the adapted laf file be put (
	probably the same as the adapted boudary or lbfd files). Will be created if nonexistent.

	delta_inp_path: Where is the input located, i.e. the single files that have been
	reviously produced by the interpolate.py routie or the regridding routines.
	These are the files called for example T00000.nc

	recompute_pressure: Boolean to indicate whether the pressure of the boundary
	files should be recomputed based on temperature changes (not necessary if a
	difference file for PP already exists (e.g. from previous cosmo simulation).

Output:
	The adapted laf file will be written to the chosen location and should directly be usable for CCLM.
"""


var_name_map = {
    'PHI':      'zg',
    'T':        'ta',
    'T_S':      'tas',
    'RELHUM':   'hur',
    'U':        'ua',
    'V':        'va',
    'PS':       'ps',
}

#var_name_map = {
#    'PHI':      'PHI',
#    'T':        'T',
#    'T_S':      'T_S',
#    'RELHUM':   'RELHUM',
#    'U':        'U',
#    'V':        'V',
#    'PS':       'PS',
#}



def lafadapt(inp_laf_path, out_laf_path, delta_inp_path,
            delta_time_step, new_time_string, laf_dt):

    laffile = xr.open_dataset(inp_laf_path, decode_cf=False)
    print('load done')

    laffile.time.attrs['units'] = new_time_string

    # compute pressure on era5 levels
    P = laffile['PS'].expand_dims(dim={'level':laffile.level})
    P = laffile.akm + P * laffile.bkm
    #laffile['lnP0'] = np.log(P)
    P_hl = laffile['PS'].expand_dims(dim={'level1':laffile.level1})
    P_hl = laffile.ak + P_hl * laffile.bk
    #laffile['P'] = P
    #laffile['P_hl'] = P_hl
    #plt.plot(P.mean(dim=['time','lon','lat']),laffile.level)
    #plt.show()
    #quit()

    # compute relative humidity in ERA climate state
    laffile['RELHUM'] = specific_to_relative_humidity(
            laffile['QV'], P, laffile['T']).transpose(
                                    'time', 'level', 'lat', 'lon')
    print('rh computed')


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







    # reference pressure
    p_ref = 10000
    p_ref = 40000
    ##p_ref = 30000
    #p_ref = 50000
    adj_factor = 0.95
    thresh_phi_ref_rmse = 0.05
    
    ## integrate geopotential upward to obtain reference geopotential
    ## in ERA climate state
    PHI_hl, phi_ref_orig, phi_ref_star = integ_geopot_era5(P_hl,
            laffile.FIS, laffile.T, laffile.QV, laffile.level1, p_ref)
    PHI_orig_100, phi_ref_orig_100, phi_ref_star_orig_100 = integ_geopot_era5(P_hl,
            laffile.FIS, laffile.T, laffile.QV, laffile.level1, 10000)
    PHI_orig_400, phi_ref_orig_400, phi_ref_star_orig_400 = integ_geopot_era5(P_hl,
            laffile.FIS, laffile.T, laffile.QV, laffile.level1, 40000)

    #laffile['PHI_hl'] = PHI_hl
    #hl_ind = None
    ##hl_ind = 80
    ##phi_ref = PHI_hl.sel(level1=hl_ind)
    ##p_ref = laffile['P_hl'].sel(level1=hl_ind)
    #print(phi_ref.mean().values)
    #p_ref_test = integ_pressure_upward_era5(laffile, phi_ref, hl_ind)
    #print(p_ref_test.mean().values)
    #print((p_ref_test-p_ref).mean().values)
    #quit()

    PHI_orig = PHI_hl.copy()
    print('geopotential computed')

    ## add climate change delta to geopotential at reference pressure level
    ## to obtain reference geopotential in future climate state
    ## (the target variable for the iteration below)
    delta_phi_ref = load_delta(delta_inp_path, var_name_map['PHI'],
                        laf_dt, delta_time_step).sel(plev=p_ref) * CON_G
    delta_phi_ref = delta_phi_ref.assign_coords(lat=laffile.lat.values)
    #phi_ref_future = phi_ref_orig + delta_phi_ref.values
    #phi_ref_future = phi_ref

    delta_phi_ref_100 = load_delta(delta_inp_path, var_name_map['PHI'],
                        laf_dt, delta_time_step).sel(plev=10000) * CON_G
    delta_phi_ref_100 = delta_phi_ref_100.assign_coords(lat=laffile.lat.values)
    delta_phi_ref_400 = load_delta(delta_inp_path, var_name_map['PHI'],
                        laf_dt, delta_time_step).sel(plev=40000) * CON_G
    delta_phi_ref_400 = delta_phi_ref_400.assign_coords(lat=laffile.lat.values)

    #print(delta_phi_ref_100.mean().values)
    #print(delta_phi_ref_400.mean().values)
    #quit()

    #phi_ref_future_100 = phi_ref_orig_100 + delta_phi_ref_100.values
    #phi_ref_future_400 = phi_ref_orig_400 + delta_phi_ref_400.values

    #phi_ref_future_100 = phi_ref_star_orig_100 + delta_phi_ref_100.values
    #phi_ref_future_400 = phi_ref_star_orig_400 + delta_phi_ref_400.values

    #laffile['FIS'] += 1000
    #phi_ref_future.to_netcdf('phi_targ.nc')



    ## interpolate climate deltas to ERA5 vertical grid and
    ## shift variables to future climate state
    var_names = ['T', 'U', 'V', 'RELHUM']
    var_names = ['T', 'RELHUM']
    #var_names = []
    for var_name in var_names:
        print('add {}'.format(var_name))
        laffile[var_name] = add_delta_interp_era5(laffile[var_name],
                        P, var_name_map[var_name], laffile, 
                        delta_inp_path, delta_time_step, laf_dt)

    #warming = 10
    #for l in range(130,138):
    #    laffile['T'].loc[dict(level=l)].values += warming


    ### compute QV of future climate (assume pressure of ERA climate)
    #laffile['QV'].values = relative_to_specific_humidity(
    #        laffile['RELHUM'], laffile['P'], laffile['T']).transpose(
    #                                'time', 'level', 'lat', 'lon').values


    ## determine future climate state surface pressure using iterative
    ## procedure
    delta_PS = xr.zeros_like(laffile['PS'])
    delta_PS_100 = xr.zeros_like(laffile['PS'])
    delta_PS_400 = xr.zeros_like(laffile['PS'])

    P_hl_orig = P_hl.copy()

    for it in range(20):

        # recompute pressure on half levels
        P_hl = P_hl_orig + delta_PS * laffile.bk
        P_hl_100 = P_hl_orig + delta_PS_100 * laffile.bk
        P_hl_400 = P_hl_orig + delta_PS_400 * laffile.bk

        # compute updated geopotential at reference pressure
        PHI_hl, phi_ref, phi_ref_star = integ_geopot_era5(P_hl,
                laffile.FIS, laffile.T, laffile.QV, laffile.level1, p_ref)
        PHI_100, phi_ref_100, phi_ref_star_100 = integ_geopot_era5(P_hl_100,
                laffile.FIS, laffile.T, laffile.QV, laffile.level1, 10000)
        PHI_400, phi_ref_400, phi_ref_star_400 = integ_geopot_era5(P_hl_400,
                laffile.FIS, laffile.T, laffile.QV, laffile.level1, 40000)

        #dPHIdPS_100 = (PHI_100-PHI_orig_100)/delta_PS_100
        #dPHIdPS_100 = (PHI_100-PHI_orig_100)
        #dPHIdPS_100.to_netcdf('dPHIdPS_100.nc')

        dphi_ref = phi_ref - phi_ref_orig
        dphi_ref_100 = phi_ref_100 - phi_ref_orig_100
        dphi_ref_400 = phi_ref_400 - phi_ref_orig_400

        #phi_ref_100 = phi_ref_star_100
        #phi_ref_400 = phi_ref_star_400

        #plt.plot((PHI_100-PHI_orig).mean(dim=['time','lon','lat']),laffile.level1)
        #plt.plot((PHI_400-PHI_orig).mean(dim=['time','lon','lat']),laffile.level1)
        #plt.show()

        print('phi')
        print(phi_ref_100.mean().values)
        print(phi_ref_400.mean().values)
        print('dphi')
        print(dphi_ref_100.mean().values)
        print(dphi_ref_400.mean().values)
        print(delta_phi_ref.mean().values)
        print(delta_phi_ref.lat)
        print(dphi_ref.lon - delta_phi_ref.lon)
        #print(delta_phi_ref)
        quit()
        print((dphi_ref - delta_phi_ref).mean().values)
        quit()

        # error in future geopotential
        dphi_ref_error = dphi_ref - delta_phi_ref
        dphi_ref_error_100 = dphi_ref_100 - delta_phi_ref
        dphi_ref_error_400 = dphi_ref_400 - delta_phi_ref

        print('phi ref deviation')
        print(dphi_ref_error_100.mean().values)
        print(dphi_ref_error_400.mean().values)
        quit()

        adj_PS = - adj_factor * (laffile['PS'] + delta_PS) / (CON_RD * 
                laffile.T.sel(level=np.max(laffile.level))) * dphi_ref_error
        del adj_PS['level']
        adj_PS_100 = - adj_factor * (laffile['PS'] + delta_PS_100) / (CON_RD * 
                laffile.T.sel(level=np.max(laffile.level))) * dphi_ref_error_100
        del adj_PS_100['level']
        adj_PS_400 = - adj_factor * (laffile['PS'] + delta_PS_400) / (CON_RD * 
                laffile.T.sel(level=np.max(laffile.level))) * dphi_ref_error_400
        del adj_PS_400['level']
        print('adj')
        print(adj_PS_100.mean().values)
        print(adj_PS_400.mean().values)
        delta_PS += adj_PS
        delta_PS_100 += adj_PS_100
        delta_PS_400 += adj_PS_400
        print('delta')
        print(delta_PS_100.mean().values)
        print(delta_PS_400.mean().values)
        print('DONE')
        #quit()
        ## DEBUG STOP

        #phi_ref.to_netcdf('phi_{:03d}.nc'.format(it).format(it))
        #phi_ref_error.to_netcdf('err_{:03d}.nc'.format(it))
        #delta_PS.to_netcdf('dps_{:03d}.nc'.format(it))
        #adj_PS.to_netcdf('adj_{:03d}.nc'.format(it))

        phi_ref_rmse = np.sqrt(np.square(phi_ref_error).mean()).values
        print('iteration {:03d}, phi rmse: {}'.format(it, phi_ref_rmse))


        if phi_ref_rmse < thresh_phi_ref_rmse:
            break

    #plt.plot((PHI_hl-PHI_orig).mean(dim=['time','lon','lat']),laffile.level1)
    #plt.show()

    laffile['PS'].values += delta_PS
    #delta_PS.to_netcdf(
    #        os.path.join(Path(out_laf_path).parents[0],
    #            'delta_ps_{:%Y%m%d%H}.nc'.format(laf_dt)))
    #delta_PS.to_netcdf(
    #            'delta_ps_{}_{:%Y%m%d%H}.nc'.format(p_ref, laf_dt))
    #print(delta_PS.mean().values)

    delta_PS_100.to_netcdf(
                'delta_ps_{}_{:%Y%m%d%H}.nc'.format(100, laf_dt))
    delta_PS_400.to_netcdf(
                'delta_ps_{}_{:%Y%m%d%H}.nc'.format(400, laf_dt))
    quit()

    ## take surface pressure change from GCM
    #var_name = 'PS'
    #print('add {}'.format(var_name))
    #add_delta_era5(var_name_map[var_name], laffile, delta_inp_path,
    #                delta_time_step, laf_dt,
    #                laf_var_name=var_name, delta_fact=-1)

    ## update surface pressure
    #PS = integ_pressure_era5(laffile)






    # update T_SKIN 
    var_name = 'T_S'
    print('add {}'.format(var_name))
    add_delta_era5(var_name_map[var_name], laffile, delta_inp_path,
                    delta_time_step, laf_dt, laf_var_name='T_SKIN')


    del laffile['RELHUM']
    del laffile['P']
    del laffile['P_hl']

    laffile.to_netcdf(out_laf_path, mode='w')
    laffile.close()
    print('saved to file {}.'.format(out_laf_path))








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
    #delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded_delta_era5_test'

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
