import xarray as xr
import sys, argparse
from pathlib import Path
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from my_lbfd_adapt import hour_of_year
import matplotlib.pyplot as plt
import numpy as np
import os
from functions import (get_alt_half_level, get_alt_full_level,
        adjust_pressure_to_new_climate,  get_pref,
        fix_grid_coord_diffs)
"""
Can be used to add the calculated difference in all necessary variables to the initial condition (laf) file.
This very fast, just run it in the console.
It requires the difference in relative humidity to adapt the specific humidity!

Input:
	lafpath: Path to the original laf-file from the "base" simulation
	(e.g. reanalysis driven or historical simulation). The name of the laf file
	must be as outputted by int2lm (e.g. laf1970010100.nc ).

	newyear: What year to use in the files (change it to the future to adapt CO2 levels)

	laftimestep: Which timestep within the annual cycle is apropriate to adapt
	the laf file? (0 for beginning of year; otherwise dayofyear*timesteps_per_day)

	newtimestring: What timestamp should be used for the adapted laf file?
	Put the exact time of the new laf file in the format 'seconds since yyyy-mm-dd hh:mm:ss'

	outputpath: In which folder should the adapted laf file be put (
	probably the same as the adapted boudary or lbfd files). Will be created if nonexistent.

	Diffspath: Where is the input located, i.e. the single files that have been
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
#Diffspath = '/scratch/snx3000/heimc/pgw/vertint_{}_compr'.format(sim_name_base)
Diffspath = '/scratch/snx3000/heimc/pgw/vertint_{}'.format(sim_name_base)
terrainpath = '/scratch/snx3000/heimc/pgw/constant_{}.nc'.format(sim_name_base)
recompute_pressure = True

year = sim_start_date.year
print(year)

lafpath = os.path.join(wd_path, 
        '{:%y%m%d}00_{}_ctrl/int2lm_out/laf{:%Y%m%d}000000.nc'.format(
                    sim_start_date, sim_name_base, sim_start_date))
newyear = year + changeyears
newtimestring = f'seconds since {newyear}-01-01 00:00:00'
outputpath = os.path.join(wd_path, 
                '{:%y%m%d}00_{}_pgw2/int2lm_out/'.format(
                    sim_start_date, sim_name_base))
laftimestring = 'seconds since 2006-01-01 00:00:00'
#init_dt = datetime(2006,8,1,0)
init_dt = sim_start_date
laftimestep = int(hour_of_year(init_dt)/args.delta_hour_inc)
print('use time step {}'.format(laftimestep))

pgw_sim_start_date = sim_start_date + relativedelta(years=changeyears)
output_laf_path = os.path.join(outputpath, 
        'laf{:%Y%m%d}000000.nc'.format(pgw_sim_start_date))
print(output_laf_path)


height_flat_half = xr.open_dataset(lafpath).vcoord #these are half levels
height_flat = xr.open_dataset(lafpath).vcoord[:-1] #to fill with full levels
vcflat=xr.open_dataset(lafpath).vcoord.vcflat

#get the full level height
height_flat.data = height_flat_half.data[1:] + \
(0.5 * (height_flat_half.data[:-1] - height_flat_half.data[1:]) )




def lafadapt(lafpath, output_laf_path, outputpath, Diffspath, laftimestep, 
            newtimestring, pref, pref_sfc, height_array, recompute_pressure):

    laffile = xr.open_dataset(lafpath, decode_cf=False)

    print('load done')
    #print('time' in laffile['T'].dims)
    #print(laffile['T'])
    #quit()

    #if output directory doesn't exist create it
    Path(outputpath).mkdir(parents=True, exist_ok=True)

    #print(lafpath)
    #print(output_laf_path)
    #print(outputpath)
    #print(Diffspath)
    #print(laftimestep)
    #print(newtimestring)
    #print(pref.shape)
    #print(pref_sfc.shape)
    #print(height_array.shape)
    #print(recompute_pressure)


    def comprelhums(laffile, pref, pref_sfc):
        p = laffile['PP'] + pref
        QV = laffile['QV']
        T = laffile['T']

        p_sfc = laffile['PP'][:,-1,:,:] + pref_sfc
        QV_S = laffile['QV_S']
        T_S = laffile['T_S']

        RH = 0.263 * p * QV *(np.exp(17.67*(T - 273.15)/(T-29.65)))**(-1)
        RH_S = 0.263 * p_sfc * QV_S *(np.exp(17.67*(T_S - 273.15)/(T_S-29.65)))**(-1)

        return RH, RH_S


    def diffadd(var, laffile=laffile):
        Diff = xr.open_dataset(
                    f'{Diffspath}/{var}{laftimestep:05d}.nc')[var]
        ## HCH2021 start
        ## There are small grid inconsistencies that have to be fixed... 
        fix_grid_coord_diffs(Diff, laffile)
        fix_grid_coord_diffs(Diff, laffile)
        ## HCH 2021 stop

        if 'time' in laffile[var].dims:
            laffile[var].data[0,::] = (laffile[var].data[0,::] + 
                                Diff.data.astype('float32'))
        else:
            laffile[var].data = (laffile[var].data + 
                                Diff.data.astype('float32'))
        #quit()












    def pressure_recompute(laf_file, pref, height_array, height_flat):
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


    #compute new humidity funcion once temperature and pressure were changed
    def computeQVnew(laffile, RH_old, RH_S_old):
        Diffrh = xr.open_dataset(f'{Diffspath}/RELHUM{laftimestep:05d}.nc')['RELHUM']
        Diffrh_s = xr.open_dataset(f'{Diffspath}/RELHUM_S{laftimestep:05d}.nc')['RELHUM_S']
        ### HCH2021 start
        ### There are small grid inconsistencies that have to be fixed... 
        fix_grid_coord_diffs(Diffrh, laffile)
        fix_grid_coord_diffs(Diffrh_s, laffile)
        ### HCH 2021 stop

        newRH = RH_old.data + Diffrh.data.astype('float32')
        newRH_S = RH_S_old.data + Diffrh_s.data.astype('float32')

        p = laffile['PP'] + pref
        T = laffile['T']
        p_sfc = np.squeeze(laffile['PP'][:,-1,:,:]) + pref_sfc
        T_S = laffile['T_S']

        newQV = (newRH.data  * np.exp(17.67*(T.data  - 273.15)/(T.data -29.65))) / ( 0.263 * p.data)
        newQV_S = (newRH_S.data  * np.exp(17.67*(T_S.data  - 273.15)/(T_S.data -29.65))) / ( 0.263 * p_sfc.data)

        laffile['QV'].data = newQV.astype('float32')
        laffile['QV_S'].data = newQV_S.astype('float32')

        return laffile


    #get relative humidity in old laf
    RH_old, RH_S_old = comprelhums(laffile, pref, pref_sfc)
    print('rh done')

    ## OLD VERSION Roman Brogli
    #if recompute_pressure == True:
    #    laffile = pressure_recompute(laffile, pref, height_array, height_flat)
    #    variables = ['T', 'T_S', 'U', 'V']
    #else:
    #    variables = ['T', 'T_S', 'U', 'V' ,'PP']

    ## HCH2021 Updated pressure computation using hydrostatic equation
    if recompute_pressure == True:
        vcoord = xr.open_dataset(terrainpath).vcoord
        alt_half_level = get_alt_half_level(vcoord, terrainpath)
        alt_full_level = get_alt_full_level(alt_half_level)
        laffile = adjust_pressure_to_new_climate(Diffspath, laftimestep,
                                laffile, 
                                pref, alt_half_level, alt_full_level)
        variables = ['T', 'T_S', 'U', 'V']
    else:
        variables = ['T', 'T_S', 'U', 'V' ,'PP']
    print('pressure done')

    #change other variables
    for var in variables:
        print('add {}'.format(var))
        diffadd(var, laffile)

    laffile.time.attrs['units'] = newtimestring
    laffile['time'].data[0] = 0

    #apply moisture function
    laffile = computeQVnew(laffile, RH_old, RH_S_old)
    print('moisture')

    laffile.to_netcdf(output_laf_path, mode='w')
    laffile.close()
    print(f'saved {output_laf_path}')



pref, pref_sfc, height_array = get_pref(vcflat, terrainpath, height_flat)
lafadapt(lafpath, output_laf_path, outputpath, Diffspath, laftimestep, newtimestring,
        pref, pref_sfc, height_array, recompute_pressure)
