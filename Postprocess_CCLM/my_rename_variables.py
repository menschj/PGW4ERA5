import xarray as xr
import sys, argparse
import numpy as np
"""
Script to rename variables if necessary; specific to COSMO-CLM
"""

## input arguments
parser = argparse.ArgumentParser(description =
                'Rename 2D fields to COSMO naming convention.')
# variable to plot
parser.add_argument('var_names', type=str,
                    help='separate multiple with ","')
# number of BC files per day
parser.add_argument('-n', '--n_per_day', type=int, default=8)
args = parser.parse_args()
print(args)


n_out_time_steps = 366 * args.n_per_day
#oldpath='/scratch/snx3000/heimc/pgw/regridded_SA_12/'
#newpath='/scratch/snx3000/heimc/pgw/vertint_SA_12/'
oldpath='/scratch/snx3000/heimc/pgw/regridded_SA_3/'
newpath='/scratch/snx3000/heimc/pgw/vertint_SA_3/'

iter_steps = np.arange(n_out_time_steps)
iter_steps = np.arange(1696,1954)

## DEBUG
#iter_steps = iter_steps[2494:]
#iter_steps = [0]

var_name_options = ['hurs', 'tas']
var_dict = {'hurs':'RELHUM_S', 'tas':'T_S'}
var_names = args.var_names.split(',')
for var_name in var_names:
    if var_name not in var_name_options:
            raise ValueError('Input var_name {} is invalid.'.format(
                            var_name))

for var_name in var_names:
    newvariable = var_dict[var_name]
    print('Rename {} to {}'.format(var_name, newvariable))

    for i in iter_steps:
        old = xr.open_dataset(f"{oldpath}/{var_name}{i:05d}.nc")
        old = old.rename({var_name:newvariable})
        old.to_netcdf(f"{newpath}/{newvariable}{i:05d}.nc")
