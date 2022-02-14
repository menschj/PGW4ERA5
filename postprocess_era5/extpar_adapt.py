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


var_name_map = {
    'tas':  'T_CL',
}




def extpar_adapt(ext_file_path, delta_inp_path):

    timer = Timer('seconds')
    timer.start('all')

    ext_file = xr.open_dataset(ext_file_path, decode_cf=False)

    # update T_C
    var_name = 'tas'
    print('update {}'.format(var_name))

    name_base='{}_delta.nc'
    delta_tas = load_delta(delta_inp_path, var_name, None, 
                            None, None)

    # TODO: debug for too small domain 
    delta_tas = delta_tas.where(~np.isnan(delta_tas), 0)

    # Make sure dimensions are exactly the same.
    # There are numerical differences between CDO remapped objects
    # and xarray data...
    delta_tas = delta_tas.assign_coords(
                    {'rlat':ext_file.rlat.values,
                     'rlon':ext_file.rlon.values})

    delta_tas_clim = delta_tas.mean(dim=['time'])
    ext_file[var_name_map[var_name]].values += delta_tas_clim

    ext_file.to_netcdf(ext_file_path, mode='w')
    ext_file.close()
    timer.stop('all')
    timer.print_report()








if __name__ == "__main__":

    ## input arguments
    parser = argparse.ArgumentParser(description =
                    'Perturb Extpar soil temperature climatology with PGW climate delta.')
    # extpar file to modify
    parser.add_argument('extpar_file_path', type=str)
    args = parser.parse_args()
    print(args)

    delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded/Emon/extpar_SA_3'
    #delta_inp_path = '/scratch/snx3000/heimc/pgw/regridded/Amon/MPI-ESM1-2-HR'

    extpar_adapt(args.extpar_file_path, delta_inp_path)

