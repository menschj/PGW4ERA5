# -*- coding: utf-8 -*- 
###############################################################################
import subprocess, os, argparse
from my_regrid_horizontal import regridhorizontal
from interpolate import interpannualcycle
#from era5_vertical import vertinterpol
from package.utilities import cd
from pathlib import Path
###############################################################################

## input arguments
parser = argparse.ArgumentParser(description =
                'Prepare PGW climate deltas for adding to BC files.')
# variable to plot
parser.add_argument('var_names', type=str,
                    help='separate multiple with ","')
# number of parallel processes
parser.add_argument('-p', '--n_par', type=int, default=1)
# number of BC files per day
parser.add_argument('-n', '--n_per_day', type=int, default=8)
args = parser.parse_args()
print(args)

performsmooth   = 0
performinterp   = 1
regridhori      = 1
regridvert      = 0 # not used anymore for era5

# number of output time steps in total
n_out_time_steps = 366 * args.n_per_day

var_name_map = {
        'hur'   :'RELHUM',
        'hurs'  :'RELHUM_S',
        'ta'    :'T', 
        'tas'   :'T_S', 
        #'pa'    :'PP',  # does not work for plev
        'ps'    :'PS', 
        'ua'    :'U', 
        'va'    :'V',
        'zg'    :'PHI',
        }
var_names = args.var_names.split(',')
for var_name in var_names:
    if var_name not in var_name_map:
            raise ValueError('Input var_name {} is invalid.'.format(
                            var_name))


if performsmooth:
#settings for timesmoothing script:

    #list of the files that contain the variables to be smoothed
    samepath = '/project/pr94/robro/inputfiles_for_surrogate_hadgem/input_github/'
    annualcycleraw = [
    samepath+'Diff_HadGEM2-ES_RCP85_PP.nc',
    samepath+'Diff_HadGEM2-ES_RCP85_QV.nc',
    samepath+'Diff_HadGEM2-ES_RCP85_QV_S.nc',
    samepath+'Diff_HadGEM2-ES_RCP85_T.nc',
    samepath+'Diff_HadGEM2-ES_RCP85_T_S.nc',
    samepath+'Diff_HadGEM2-ES_RCP85_T_SO.nc',
    samepath+'Diff_HadGEM2-ES_RCP85_U.nc',
    samepath+'Diff_HadGEM2-ES_RCP85_V.nc',
    samepath+'Diff_HadGEM2-ES_RCP85_RELHUM.nc',
    samepath+'Diff_HadGEM2-ES_RCP85_RELHUM_S.nc'
    ]
    #list of variablenames
    variablename_to_smooth = ['PP', 'QV', 'QV_S', 'T', 'T_S', 'T_SO', 'U', 'V','RELHUM','RELHUM_S']
    #path to put the output netcdf
    outputpath = '/scratch/snx3000/robro/pgwtemp/smoothed/'


    #enter the command to run the script:
    for num,pathi in enumerate(annualcycleraw):
            commandsmooth = f"python timesmoothing.py {pathi} {variablename_to_smooth[num]} {outputpath} > outputfile_smooth.txt &"
            subprocess.run(commandsmooth, shell=True)
            print(commandsmooth)
    

if performinterp:
    ###########################################################################
    ### Namelist
    ###########################################################################
    #see documentation in interpolate.py
    # run e.g. cdo sellonlatbox,-65,35,-45,30
    gcm_data_path='/scratch/snx3000/heimc/pgw/deltas/MPI-ESM1-2-HR/plev'
    gcm_delta_base='plev_{}_delta.nc'
    gcm_data_freq = 'month'

    out_path = '/scratch/snx3000/heimc/pgw/interpolated_plev/'
    ###########################################################################

    for var_name in var_names:  
        print('time interpolation {}'.format(var_name))
        gcm_file_path = os.path.join(gcm_data_path, 
                                gcm_delta_base.format(var_name))
        interpannualcycle(gcm_file_path, var_name, 
                        var_name_map[var_name],
                        n_out_time_steps,
                        gcm_data_freq, out_path, args.n_par)


if regridhori:
    ###########################################################################
    ### Namelist
    ###########################################################################
    infolder = '/scratch/snx3000/heimc/pgw/interpolated_plev/'

    #outputfolder = '/scratch/snx3000/heimc/pgw/regridded_alt2_SA_12/'
    #out_grid_file = 'target_grid_SA_12'

    #outputfolder = '/scratch/snx3000/heimc/pgw/regridded_alt2_SA_3/'
    #out_grid_file = 'target_grid_SA_3'

    outputfolder = '/scratch/snx3000/heimc/pgw/regridded_era5/'
    out_grid_file = 'target_grid_era5'
    ###########################################################################

    #get the python command and write a file to submit to the piz daint machine
    #comandreghor = f"python regrid_horizontal.py {infolder} {variable} {n_out_time_steps} {outputgridfile} {outputfolder} &> out_regrid.txt &" 
    #subprocess.run(comandreghor, shell=True)
    for var_name in var_names:
        print('regrid horizontal {}'.format(var_name))
        regridhorizontal(infolder, var_name_map[var_name], n_out_time_steps,
                        out_grid_file, outputfolder, njobs=args.n_par)



#this part is software/hardware specific for the piz daint supercomputer on CSCS
if regridvert:

    #constant_path = '/scratch/snx3000/heimc/pgw/constant_SA_12.nc'
    #datapath = '/scratch/snx3000/heimc/pgw/regridded_alt2_SA_12/'
    #outputpath = '/scratch/snx3000/heimc/pgw/vertint_alt2_SA_12/'

    #constant_path = '/scratch/snx3000/heimc/pgw/constant_SA_3.nc'
    #datapath = '/scratch/snx3000/heimc/pgw/regridded_alt2_SA_3/'
    #outputpath = '/scratch/snx3000/heimc/pgw/vertint_alt2_SA_3/'

    constant_path = '/scratch/snx3000/heimc/pgw/constant_era5.nc'
    datapath = '/scratch/snx3000/heimc/pgw/regridded_era5/'
    outputpath = '/scratch/snx3000/heimc/pgw/vertint_era5/'

    i_submit = 0
    submit_dir = '/scratch/snx3000/heimc/pgw/submit/'
    #vcflat = 17827 #height where modellevels become flat
    steps_per_job = 10 #split the job into multiple chucks and run in paralell
    lasttime = n_out_time_steps
    #starttime = 0
    starttime = 1696
    lasttime = starttime + 21 * 8

    # copy heights file and script to submission directory
    if i_submit:
        with cd(submit_dir):
            #subprocess.run('cp /project/pr94/heimc/lmp_template/submodules/pgw-python/heights.txt .', shell=True)
            subprocess.run('cp /project/pr94/heimc/lmp_template/submodules/pgw-python/era5_vertical.py .', shell=True)
            subprocess.run('cp /project/pr94/heimc/lmp_template/submodules/pgw-python/functions.py .', shell=True)
            subprocess.run('cp /project/pr94/heimc/lmp_template/submodules/pgw-python/constants.py .', shell=True)

    for var_name in var_names:
        print(var_name)
        #outvar = outvar_dict[var_name]
        for start_job in range(starttime, lasttime, steps_per_job):
            end_job = min(start_job + steps_per_job, lasttime-1)

            print('start at {}'.format(start_job))
            print('end at {}'.format(end_job))
            
            if i_submit == 0:
                ###TODO debug only
                end_job = start_job + 1 # for debugging purpose only run one time step
                vertinterpol(constant_path, datapath, var_name_map[var_name], 
                            outputpath, end_job, start_job)
                break
            else:
                comandregver = f"srun -u python era5_vertical.py {constant_path} {datapath} {var_name_map[var_name]} {outputpath} {end_job} {start_job}"

                #create a run script for afew timesteps and each variable. 
                with open (f'/scratch/snx3000/heimc/pgw/submit/submit_{var_name}.bash', 'w') as rsh:
                    rsh.write(f'''\
#!/bin/bash -l
#SBATCH --job-name="{var_name}_{start_job}"
#SBATCH --time=01:30:00
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=gpu
#SBATCH --hint=nomultithread
#SBATCH --account=pr94
#SBATCH --output=bOut/out_{var_name}_{start_job}
#SBATCH --error=bOut/err_{var_name}_{start_job}

module load daint-gpu

export PATH="/project/pr94/heimc/anaconda3/bin:$PATH"
source /project/pr94/heimc/anaconda3/etc/profile.d/conda.sh

# make sure xarray does not explode with number of threads
export OMP_NUM_THREADS=1

{comandregver}

''')
                #submit the slurm batch job script from the scratch directory (change directory if needed)
                with cd(submit_dir):
                    subprocess.run(f'sbatch submit_{var_name}.bash', shell=True)
