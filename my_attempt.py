import os
import numpy as np
from my_cclm_vertical import vertinterpol




samepath='/scratch/snx3000/heimc/pgw/deltas/GCMdata_SA/'

file_path_int = [
samepath+'diff_hur.nc',
#samepath+'diff_hurs.nc',
#samepath+'diff_ta.nc',
#samepath+'diff_tas.nc',
#samepath+'diff_ua.nc',
#samepath+'diff_va.nc',
]

#variablename = ['hur', 'hurs', 'ta', 'tas', 'ua', 'va']
variablename = ['hur']

output_time_steps = 366 * 2
inputfreq = 'month'
outputpath_int = '/scratch/snx3000/heimc/pgw/interpolated_my/'

outputgridfile = '/project/pr94/heimc/data/cosmo_out/SA_3_long/lm_c/lffd20060801000000c.nc'
outputfolder_regrid = '/scratch/snx3000/heimc/pgw/regridded_debug'



### VERTICAL INTERPOLATION
terrainpath = outputgridfile
datapath = outputfolder_regrid

variablename = ['hur','ta','ua','va'] #only 3D data
out_vars = ['RELHUM', 'T', 'U', 'V'] #rename variable to the correct ones for cosmo
outputpath_vertical = '/scratch/snx3000/heimc/pgw/deltas/GCMdata_SA_vint/'

vcflat = 17827 #height where modellevels become flat (see cosmo namelist)

#steps_per_job = inputtimesteps + 100 #this option is not used in notebook just leave as it is
#starttime = 0 #set to 0; not used in notebook

print('The following height half-levels will be used, should match namelist')
print(np.genfromtxt('heights.txt',skip_header=1)[:,1])

var_ind = 0
inp_var_key = variablename[var_ind]
out_var_key = out_vars[var_ind]
inp_file_path = os.path.join('/scratch/snx3000/heimc/pgw/deltas/GCMdata_SA/','diff_hur.nc')
out_file_path = os.path.join(outputpath_vertical, out_var_key+'.nc')
vertinterpol(terrainpath, datapath, \
            vcflat, out_var_key, inp_var_key, out_file_path, inp_file_path)
