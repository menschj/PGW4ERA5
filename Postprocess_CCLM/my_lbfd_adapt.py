import xarray as xr
import sys, argparse
import os
import glob
from pathlib import Path
import numpy as np
from datetime import datetime,timedelta
from dateutil.relativedelta import relativedelta
from base.functions import (get_alt_half_level, get_alt_full_level,
        adjust_pressure_to_new_climate,  get_pref,
        fix_grid_coord_diffs)

"""
Add the calculated difference in all necessary variables to the boundary condition (lbfd) files for one year.
This should be run on a cluster.

Input:
	year (comandline argument): The year for which the boundary files should be adapted.
	The year refers to the year in the original output from int2lm.

	lbfdpath: path to directory where boundary files are saved

	changeyears: amount of years to be added to timestamp of data (shift to future)

	outputpath: path to put the modifiyed bondary files

	Diffspath: Where to find the changes (climate change deltas) to add to the boundary files.
	This is the ouput of earlier scripts in this repository (i.e. T00000.nc, T00001.nc etc.).

	difftimesteps: Amount of timesteps (or boundary data fields) in one year
	(extrapolate to an entire year even if only a fraction is needed; this depends on the calendar used)

    terrainpath: Path to a netcdf file containing the height of the terrain in
	the cosmo domain (could be a constant file such as lffd1969120100c.nc)

	starttimestep: The files as produced from previous scripts start in january
	and run through the year. If your PGW simulation does not start in january,
	you have to compute the timestep of the start (dayofyear * timesteps_per_day)

	recompute_pressure: Boolean to indicate whether the pressure of the boundary
	files should be recomputed based on temperature changes (not necessary if a
	difference file for PP already exists (e.g. from previous cosmo simulation).

Output:
	For every boundary file in the inputdata, a corresponding output field will
	be written to the path specified as outputpath.
"""

def hour_of_year(dt): 
    beginning_of_year = datetime(dt.year, 1, 1, tzinfo=dt.tzinfo)
    return(int((dt - beginning_of_year).total_seconds() // 3600))




##get reference pressure function
#def getpref(vcflat, terrainpath, height_flat):
#    smoothing = (vcflat - height_flat) / vcflat
#    smoothing = np.where(smoothing > 0, smoothing, 0)
#
#    const = xr.open_dataset(terrainpath)
#
#    hsurf = const['HSURF'].squeeze()
#
#    #the height at which the reference pressure needs to be computed needs to be derived form the terrain   following coordinates:
#    newheights = np.zeros((len(height_flat), hsurf.shape[0], hsurf.shape[1]))
#
#    #avoid forloop
#    newheights = height_flat.values[:,None,None] + hsurf.values[None,:,:] * smoothing[:,None,None]
#
#    #New formulation as researched by Christian Steger (untested)
#    # Constants
#    p0sl = height_flat.p0sl # sea-level pressure [Pa]
#    t0sl = height_flat.t0sl   # sea-level temperature [K]
#    # Source: COSMO description Part I, page 29
#    g = 9.80665     # gravitational acceleration [m s-2]
#    R_d = 287.05    # gas constant for dry air [J K-1 kg-1]
#    # Source: COSMO source code, data_constants.f90
#
#    # irefatm = 2
#    delta_t = height_flat.delta_t
#    h_scal = height_flat.h_scal
#    # Source: COSMO description Part VII, page 66
#    t00 = t0sl - delta_t
#
#    pref = p0sl * np.exp (-g / R_d * h_scal / t00 * \
#               np.log((np.exp(newheights / h_scal) * t00 + delta_t) / \
#                      (t00 + delta_t)) )
#    pref_sfc = p0sl * np.exp (-g / R_d * h_scal / t00 * \
#               np.log((np.exp(hsurf.data / h_scal) * t00 + delta_t) / \
#                      (t00 + delta_t)) )
#
#    return pref, pref_sfc, newheights



#function to adapt all lbfd files:
def lbfdadapt(lbfdpath, outputpath, Diffspath, lbfd_dts, delta_hour_inc, 
                changeyears, pref, pref_sfc, dz, 
                recompute_pressure, height_flat):

    #function to add all variables but humidity to the boundary field (use given timestep)
    def diffadd(var, num, lbfd):
        Diff = xr.open_dataset(f'{Diffspath}/{var}{num:05d}.nc')[var]
        ## HCH2021 start
        ## There are small grid inconsistencies that have to be fixed... 
        fix_grid_coord_diffs(Diff, lbfd)
        fix_grid_coord_diffs(Diff, lbfd)
        ## HCH 2021 stop

        lbfd[var].data = lbfd[var].data + Diff.data.astype('float32')

    #function to calculate relative humidity
    def comprelhums(lbfd, pref, pref_sfc):
        p = lbfd['PP'] + pref
        QV = lbfd['QV']
        T = lbfd['T']

        p_sfc = lbfd['PP'][:,-1,:,:] + pref_sfc
        QV_S = lbfd['QV_S']
        T_S = lbfd['T_S']

        RH = 0.263 * p * QV *(np.exp(17.67*(T - 273.15)/(T-29.65)))**(-1)
        RH_S = 0.263 * p_sfc * QV_S *(np.exp(17.67*(T_S - 273.15)/(T_S-29.65)))**(-1)

        return RH, RH_S


    def pressure_recompute(lbfd, num, pref, dz, height_flat):
        #function to compute pressure field in a differen climate using the barometric
        #formula (maintaining hydrostatic balance)
        #temperature changes
        dT_sfc = xr.open_dataset(f'{Diffspath}/T_S{num:05d}.nc')['T_S']
        dT_atmos = xr.open_dataset(f'{Diffspath}/T{num:05d}.nc')['T']
        ## HCH2021 start
        ## There are small grid inconsistencies that have to be fixed... 
        fix_grid_coord_diffs(dT_sfc, lbfd)
        fix_grid_coord_diffs(dT_atmos, lbfd)
        ## HCH 2021 stop

        #get pressure field
        pressure_original = lbfd['PP'] + pref
        pressure_new = pressure_original.copy()

        temperature = lbfd['T']
        sfc_temperature = lbfd['T_S']

        #define barometric height formula
        def barometric(reference_pressure, reference_temperature, 
                        dz, lapse_rate):
                R = 8.3144598 #universal gas constant
                M = 0.0289644 # molar mass of air
                g = 9.80665
                #lapse_rate = - 0.0065
                exo = - g * M / (R * lapse_rate) #exponent in barometric formula

                pressure = reference_pressure * ( (reference_temperature + (lapse_rate * dz))
                / reference_temperature )**exo

                return pressure

        #compute surface pressure
        surface_press = barometric(pressure_original[:,-1,:,:], temperature[:,-1,:,:], -height_flat[-1], 0.0065)

        #get the lowest model level in warmer climate
        pressure_new[:,-1,:,:] = barometric(surface_press, sfc_temperature+dT_sfc, height_flat[-1], -0.0065)
        #get the rest (loop from ground up)
        for level in range(len(dz)-1, -1, -1):
                pressure_new[:,level,:,:] = barometric(pressure_new[:,level+1,:,:], \
                temperature[:,level+1,:,:]+dT_atmos[level+1,:,:], dz[level,:,:], -0.0065)

        new_pp = pressure_new.data - pref
        #convert to PP
        lbfd['PP'].data = new_pp.astype('float32')

        return lbfd


    #compute new humidity funcion once temperature and pressure were changed
    def computeQVnew(lbfd, num, RH_old, RH_S_old):
        Diffrh = xr.open_dataset(f'{Diffspath}/RELHUM{num:05d}.nc')['RELHUM']
        Diffrh_s = xr.open_dataset(f'{Diffspath}/RELHUM_S{num:05d}.nc')['RELHUM_S']
        ## HCH2021 start
        ## There are small grid inconsistencies that have to be fixed... 
        fix_grid_coord_diffs(Diffrh, lbfd)
        fix_grid_coord_diffs(Diffrh_s, lbfd)
        ## HCH 2021 stop

        newRH = RH_old.data + Diffrh.data.astype('float32')
        newRH_S = RH_S_old.data + Diffrh_s.data.astype('float32')

        p = lbfd['PP'] + pref
        T = lbfd['T']
        p_sfc = lbfd['PP'][:,-1,:,:] + pref_sfc
        T_S = lbfd['T_S']

        newQV = (newRH.data  * np.exp(17.67*(T.data  - 273.15)/(T.data -29.65))) / ( 0.263 * p.data)
        newQV_S = (newRH_S.data  * np.exp(17.67*(T_S.data  - 273.15)/(T_S.data -29.65))) / ( 0.263 * p_sfc.data)

        lbfd['QV'].data = newQV.astype('float32')
        lbfd['QV_S'].data = newQV_S.astype('float32')

        return lbfd



    # calculation part
    #get a list of all lbfd files:
    os.chdir(lbfdpath)
    files = glob.glob('lbfd??????????????.nc')
    files.sort()

    #if output directory doesn't exist create it
    Path(outputpath).mkdir(parents=True, exist_ok=True)


    #loop over all boundary fields:
    for lbfd_dt in lbfd_dts:
        #print and open boundary data
        delta_num = int(hour_of_year(lbfd_dt)/delta_hour_inc)
        #delta_num = 0
        print(delta_num, lbfd_dt)
        inp_lbfd_file = os.path.join(lbfdpath, 
                    'lbfd{:%Y%m%d%H%M%S}.nc'.format(lbfd_dt))
        out_lbfd_dt = lbfd_dt.replace(year=lbfd_dt.year + changeyears)
        out_lbfd_file = os.path.join(outputpath, 
                    'lbfd{:%Y%m%d%H%M%S}.nc'.format(out_lbfd_dt))
        
        print('compute change for lbfd file: {}'.format(inp_lbfd_file))
        print('output file: {}'.format(out_lbfd_file))

        if os.path.exists(inp_lbfd_file):
            lbfd = xr.open_dataset(inp_lbfd_file, decode_cf=False)

            #run the defined functions and change filename & time:
            RH_old, RH_S_old = comprelhums(lbfd, pref, pref_sfc)
            print('RH done')

            ## OLD VERSION Roman Brogli
            #if recompute_pressure == True:
            #    lbfd = pressure_recompute(lbfd, delta_num, pref, dz, height_flat)
            #    variables = ['T', 'T_S', 'U', 'V']
            #else:
            #    variables = ['T', 'T_S', 'U', 'V', 'PP']

            ## HCH2021 Updated pressure computation using hydrostatic equation
            if recompute_pressure == True:
                vcoord = xr.open_dataset(terrainpath).vcoord
                alt_half_level = get_alt_half_level(vcoord, terrainpath)
                alt_full_level = get_alt_full_level(alt_half_level)
                lbfd = adjust_pressure_to_new_climate(Diffspath, delta_num,
                                        lbfd, 
                                        pref, alt_half_level, alt_full_level)
                variables = ['T', 'T_S', 'U', 'V']
            else:
                variables = ['T', 'T_S', 'U', 'V' ,'PP']
            print('pressure done')

            for var in variables:
                print('add {}'.format(var))
                diffadd(var, delta_num, lbfd)

            #change time to future
            endtimestring = lbfd.time.units[-15:]
            old_yyyy_timestamp = int(lbfd.time.units[-19:-15])
            new_yyyy_timestamp = old_yyyy_timestamp + changeyears
            lbfd.time.attrs['units'] = f'seconds since {new_yyyy_timestamp}{endtimestring}'

            lbfd = computeQVnew(lbfd, delta_num, RH_old, RH_S_old)
            print('QV done')

            lbfd.to_netcdf(out_lbfd_file, mode='w')
            lbfd.close()
            print('save done')

            # sanity check
            print('sanity check..')
            ds = xr.open_dataset(out_lbfd_file)
            var_names = ['U','V','T','QV','PP','T_S','QV_S']
            for var_name in var_names:
                if np.sum(np.isnan(ds[var_name])).values > 0:
                    ds.close()
                    os.remove(out_lbfd_file)
                    raise ValueError('File {} var {} broken!'.format(
                                    out_lbfd_file, var_name))
            ds.close()
            print('..completed!')
        else:
            print('#############################################################')
            print('ATTENTION: file {} does not exist. SKIPPING!'.format(inp_lbfd_file))
            print('#############################################################')





		
def workflow_lbfd(lbfdpath, outputpath, Diffspath, 
                terrainpath, lbfd_dts, delta_hour_inc,
                changeyears, recompute_pressure):

    #read height coordinate from file
    os.chdir(lbfdpath)
    files = glob.glob('lbfd??????????????.nc')
    height_flat_half = xr.open_dataset(files[0]).vcoord #these are half levels
    height_flat = xr.open_dataset(files[0]).vcoord[:-1]

    vcflat=xr.open_dataset(files[0]).vcoord.vcflat

    #get the full level height
    height_flat.data = height_flat_half.data[1:] + \
    (0.5 * (height_flat_half.data[:-1] - height_flat_half.data[1:]) )
    
    pref, pref_sfc, height_array = get_pref(vcflat, terrainpath, height_flat)
    dz = height_array[:-1] - height_array[1:] #get height difference between model levels
    lbfdadapt(lbfdpath, outputpath, Diffspath, lbfd_dts, delta_hour_inc,
            changeyears, pref,
            pref_sfc, dz, recompute_pressure, height_flat)

	
if __name__ == "__main__":

    ## input arguments
    parser = argparse.ArgumentParser(description =
                    'Perturb COSMO boundary conditions with PGW climate deltas.')
    ## variable to plot
    # delta hour increments
    parser.add_argument('-d', '--delta_hour_inc', type=int, default=3)
    # first date 
    parser.add_argument('-f', '--first_dt', type=str, default='2006080100')
    # last date 
    parser.add_argument('-l', '--last_dt', type=str, default='2006090100')
    # sim start date 
    parser.add_argument('-s', '--sim_start_date', type=str, default='20060801')
    # sim name excluding ctrl/pgw
    parser.add_argument('-n', '--sim_name_base', type=str, default='SA_3')
    args = parser.parse_args()
    print(args)

    #first_date = datetime.strptime(args.first_date, '%Y%m%d')
    #last_date = datetime.strptime(args.last_date, '%Y%m%d')
    first_dt = datetime.strptime(args.first_dt, '%Y%m%d%H')
    last_dt = datetime.strptime(args.last_dt, '%Y%m%d%H')
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

    lbfdpath = os.path.join(wd_path, 
                    '{:%y%m%d}00_{}_ctrl/int2lm_out/'.format(
                        sim_start_date, sim_name_base))
    outputpath = os.path.join(wd_path, 
                    '{:%y%m%d}00_{}_pgw2/int2lm_out/'.format(
                        sim_start_date, sim_name_base))
    print(lbfdpath)
    print(terrainpath)
    print(outputpath)



    newyear = year + changeyears


    delta_hour_inc = args.delta_hour_inc
    timestep = timedelta(hours=delta_hour_inc)

    #first_dt = first_date
    #last_dt = last_date + timedelta(days=1) - timedelta(hours=delta_hour_inc)


    lbfd_dts = np.arange(first_dt, last_dt+timestep, timestep).tolist()
    print('Run for dates {}'.format(lbfd_dts))
    #quit()

    #for dt in lbfd_dts:
    #    delta_num = int(hour_of_year(dt)/delta_hour_inc)
    #    print(delta_num)
    #quit()

    workflow_lbfd(lbfdpath, outputpath, Diffspath, 
            terrainpath, lbfd_dts, delta_hour_inc,
            changeyears, recompute_pressure)
