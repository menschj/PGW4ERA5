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
from package.utilities import Timer
from cdo import *


delta_var_name_map = {
    'T_CL':      'tas',
}




def extpar_adapt(ext_file_path, delta_inp_path):

    timer = Timer('seconds')
    timer.start('all')

    cdo = Cdo()

    ext_file = xr.open_dataset(ext_file_path, decode_cf=False)

    # update T_C
    var_name = 'T_CL'
    print('update {}'.format(var_name))

    name_base='{}_delta.nc'
    delta_inp_file = os.path.join(delta_inp_path,
                            name_base.format(delta_var_name_map['T_CL']))
    griddes_file = 'grid_extpar'
    ofile = cdo.remapbil(griddes_file, input=delta_inp_file, options='-f nc')
    delta_T_SO_clim = xr.open_dataset(ofile)
    print(delta_T_SO_clim.tas.shape)
    quit()
    delta_T_SO_clim = xr.open_dataset()

    delta_T_SO_clim = delta_T_SO_clim.assign_coords({'lat':ext_file.rlat.values})
    delta_T_SO_clim = delta_T_SO_clim.mean(dim=['time'])

    print(delta_T_SO_clim)
    quit()
    delta_T_SO = (
            delta_T_SO_clim + np.exp(-laffile.soil1/2.8) * 
                    (delta_T_SKIN - delta_T_SO_clim)
    )
    delta_T_SO = delta_T_SO.transpose('time', 'soil1', 'lat', 'lon')
    laffile.T_SO.values += delta_T_SO

    laffile['PS'] = PS_pgw

    laffile.to_netcdf(out_laf_path, mode='w')
    laffile.close()
    timer.stop('all')
    timer.print_report()








if __name__ == "__main__":

    ## input arguments
    parser = argparse.ArgumentParser(description =
                    'Perturb Extpar soil temperature climatology with PGW climate delta.')
    # extpar file to modify
    parser.add_argument('-f', '--extpar_file_path', type=str, default=None)
    args = parser.parse_args()
    print(args)


    delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded/Emon/MPI-ESM1-2-HR'
    #delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded/Amon/MPI-ESM1-2-HR'

    extpar_adapt(args.extpar_file_path, delta_inp_path)

