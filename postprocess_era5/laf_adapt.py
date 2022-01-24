import xarray as xr
import sys, argparse
from pathlib import Path
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
import numpy as np
import os
from base.functions import (
        hour_of_year,
        specific_to_relative_humidity,
        relative_to_specific_humidity,
        add_delta_era5,
        add_delta_interp_era5)
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

	terrainpath: Path to a netcdf file containing the height of the terrain in
	the cosmo domain (could be a constant file such as lffd1969120100c.nc)

	recompute_pressure: Boolean to indicate whether the pressure of the boundary
	files should be recomputed based on temperature changes (not necessary if a
	difference file for PP already exists (e.g. from previous cosmo simulation).

Output:
	The adapted laf file will be written to the chosen location and should directly be usable for CCLM.
"""




def lafadapt(inp_laf_path, out_laf_path, delta_inp_path,
            delta_time_step, new_time_string):

    laffile = xr.open_dataset(inp_laf_path, decode_cf=False)
    print('load done')

    laffile.time.attrs['units'] = new_time_string

    # compute pressure on era5 levels
    P0 = laffile['PS'].expand_dims(dim={'level':laffile.level})
    P0 = laffile.akm + P0 * laffile.bkm
    laffile['lnP0'] = np.log(P0)
    #plt.plot(P0.mean(dim=['time','lon','lat']),laffile.level)
    #plt.show()
    #quit()

    # get relative humidity in old laf
    RH0 = specific_to_relative_humidity(
                    laffile['QV'], P0, laffile['T']).transpose(
                                    'time', 'level', 'lat', 'lon')
    laffile['RELHUM'] = RH0
    #plt.plot(RH0.mean(dim=['time','lon','lat']),
    #        P0.mean(dim=['time','lon','lat']))
    #plt.show()
    #quit()
    print('rh done')

    # interpolate deltas to era5 vertical grid and
    # shift variables to future climate
    var_names = ['T', 'U', 'V', 'RELHUM']
    #var_names = ['T', 'U', 'V']
    for var_name in var_names:
        print('add {}'.format(var_name))
        add_delta_interp_era5(var_name, laffile, 
                        delta_inp_path, delta_time_step)

    # update surface pressure
    var_name = 'PS'
    print('add {}'.format(var_name))
    add_delta_era5(var_name, laffile, delta_inp_path, delta_time_step)


    # compute future climate pressure on era5 levels
    P1 = laffile['PS'].expand_dims(dim={'level':laffile.level})
    P1 = laffile.akm + P1 * laffile.bkm
    laffile['P1'] = P1


    #apply moisture function
    #QV0 = laffile['QV'].copy()
    QV1 = relative_to_specific_humidity(
                    laffile['RELHUM'], P0, laffile['T']).transpose(
                                    'time', 'level', 'lat', 'lon')
    laffile['QV'] = QV1

    # update T_SKIN 
    var_name = 'T_S'
    print('add {}'.format(var_name))
    add_delta_era5(var_name, laffile, delta_inp_path,
                    delta_time_step, laf_var_name='T_SKIN')

    #plt.plot(QV0.isel(time=0).mean(dim=['lon','lat']),
    #        np.exp(laffile['lnP0'].isel(time=0)).mean(dim=['lon','lat']))
    #plt.plot(QV1.isel(time=0).mean(dim=['lon','lat']),
    #        np.exp(laffile['lnP0'].isel(time=0)).mean(dim=['lon','lat']))
    #plt.show()

    del laffile['RELHUM']
    del laffile['lnP0']
    del laffile['P1']

    laffile.to_netcdf(out_laf_path, mode='w')
    laffile.close()
    print(f'saved {out_laf_path}')







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
    terrainpath = '/scratch/snx3000/heimc/pgw/constant_era5.nc'

    pgw_sim_name_ending = 'pgw6'


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
                delta_time_step, new_time_string)
