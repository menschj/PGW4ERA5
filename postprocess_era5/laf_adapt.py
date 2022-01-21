import xarray as xr
import sys, argparse
from pathlib import Path
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
import numpy as np
import os
from base.functions import (
        specific_to_relative_humidity,
        add_delta_interp_era5,
        hour_of_year,
        get_alt_half_level, get_alt_full_level,
        pressure_recompute, get_pref_old,
        adjust_pressure_to_new_climate, 
        adjust_pressure_to_new_climate2, 
        get_pref, get_pref_sfc,
        fix_grid_coord_diffs, add_delta)
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




def lafadapt(inp_laf_path, out_laf_path, delta_inp_path, delta_time_step, 
            new_time_string, recompute_pressure):


    #compute new humidity funcion once temperature and pressure were changed
    def computeQVnew(laffile, RH_old, RH_S_old):
        Diffrh = xr.open_dataset(f'{delta_inp_path}/RELHUM{delta_time_step:05d}.nc')['RELHUM']
        Diffrh_s = xr.open_dataset(f'{delta_inp_path}/RELHUM_S{delta_time_step:05d}.nc')['RELHUM_S']
        ### HCH2021 start
        ### There are small grid inconsistencies that have to be fixed... 
        fix_grid_coord_diffs(Diffrh, laffile)
        fix_grid_coord_diffs(Diffrh_s, laffile)
        ### HCH 2021 stop

        newRH = RH_old.data + Diffrh.data.astype('float32')
        newRH_S = RH_S_old.data + Diffrh_s.data.astype('float32')

        p = laffile['PP'] + pref
        T = laffile['T']
        #p_sfc = np.squeeze(laffile['PP'][:,-1,:,:]) + pref_sfc
        p_sfc = np.squeeze(laffile['PP'][:,-1,:,:]) + pref_hl.isel(level1=-1)
        T_S = laffile['T_S']

        newQV = (newRH.data  * np.exp(17.67*(T.data  - 273.15)/(T.data -29.65))) / ( 0.263 * p.data)
        newQV_S = (newRH_S.data  * np.exp(17.67*(T_S.data  - 273.15)/(T_S.data -29.65))) / ( 0.263 * p_sfc.data)

        laffile['QV'].data = newQV.astype('float32')
        laffile['QV_S'].data = newQV_S.astype('float32')

        return laffile


    laffile = xr.open_dataset(inp_laf_path, decode_cf=False)
    print('load done')

    laffile.time.attrs['units'] = new_time_string

    # compute pressure
    P0 = laffile['PS'].expand_dims(dim={'level':laffile.level})
    P0 = laffile.akm + P0 * laffile.bkm
    laffile['lnP0'] = np.log(P0)
    #plt.plot(P0.mean(dim=['time','lon','lat']),laffile.level)
    #plt.show()
    #quit()

    # get relative humidity in old laf
    RH0 = specific_to_relative_humidity(
                    laffile['QV'], P0, laffile['T'])
    #plt.plot(RH0.mean(dim=['time','lon','lat']),
    #        P0.mean(dim=['time','lon','lat']))
    #plt.show()
    #quit()
    print('rh done')

    # interpolate deltas to era5 vertical grid and
    # shift variables to future climate
    var_names = ['T', 'U', 'V']
    for var_name in var_names:
        print('add {}'.format(var_name))
        add_delta_interp_era5(var_name, laffile, 
                        delta_inp_path, delta_time_step)
    quit()

    # Update surface pressure

    #apply moisture function
    laffile = computeQVnew(laffile, RH_old, RH_S_old)
    print('moisture')

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
    parser.add_argument('-s', '--sim_start_date', type=str, default='20060801')
    # sim name excluding ctrl/pgw
    parser.add_argument('-n', '--sim_name_base', type=str, default='SA_3')
    args = parser.parse_args()
    print(args)

    sim_start_date = datetime.strptime(args.sim_start_date, '%Y%m%d')
    sim_name_base = args.sim_name_base
    wd_path = '/scratch/snx3000/heimc/lmp/wd'
    changeyears = 0
    delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded_era5'
    terrainpath = '/scratch/snx3000/heimc/pgw/constant_era5.nc'

    recompute_pressure = 'barometric'
    recompute_pressure = 'hydrostatic'
    #recompute_pressure = False

    pgw_sim_name_ending = 'pgw6'

    pgw_sim_start_date = sim_start_date + relativedelta(years=changeyears)
    new_time_string = 'seconds since {:%Y-%m-%d %H:%M:%S}'.format(pgw_sim_start_date)

    inp_laf_path = os.path.join(wd_path, 
            '{:%y%m%d}00_{}_ctrl/int2lm_in/cas{:%Y%m%d}000000.nc'.format(
                        sim_start_date, sim_name_base, sim_start_date))

    out_laf_dir = os.path.join(wd_path, 
                    '{:%y%m%d}00_{}_{}/int2lm_in/'.format(
                        sim_start_date, sim_name_base, pgw_sim_name_ending))
    #if output directory doesn't exist create it
    Path(out_laf_dir).mkdir(parents=True, exist_ok=True)

    #init_dt = datetime(2006,8,1,0)
    init_dt = sim_start_date
    delta_time_step = int(hour_of_year(init_dt)/args.delta_hour_inc)
    print('use time step {}'.format(delta_time_step))

    out_laf_path = os.path.join(out_laf_dir, 
            'cas{:%Y%m%d}000000.nc'.format(pgw_sim_start_date))
    print(out_laf_path)

    lafadapt(inp_laf_path, out_laf_path, delta_inp_path,
            delta_time_step, new_time_string,
            recompute_pressure)
