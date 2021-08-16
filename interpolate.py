# -*- coding: utf-8 -*-
###############################################################################
import xarray as xr
import sys
import numpy as np
from pathlib import Path
import gc
import dask
from package.mp import IterMP
###############################################################################

def save_time_step(outfile, out_path, var_name, filenum):
    file_path = f"{out_path}/{var_name}{filenum:05d}.nc"
    outfile[filenum].to_netcdf(file_path, mode='w')       

    # sanity check
    ds = xr.open_dataset(file_path)
    if np.sum(np.isnan(ds[var_name])).values > 0:
        ds.close()
        os.remove(out_lbfd_file)
        raise ValueError('File {} var {} broken!'.format(
                        out_lbfd_file, var_name))

def interpannualcycle(gcm_file_path, var_name, n_out_time_steps,
                    gcm_data_freq, out_path, n_par):
    """
    Interpolate temporally an annual cyle with sparse data points 
    linearly to higher time frequencies.
    This can be used to interpolate a daily or monthly resolution annual 
    cycle to the frequency that is used to update the boundary conditions. 

    Input:
    gcm_file_path:  Path to a netcdf file of the annual cycle of the signal 
                    to be interpolated.  
    var_name:       The name of the variable within the given netcdf file
                    which will be interpolated.
    n_out_time_steps: 
                    Amount of thimesteps needed as output (e.g. 366 * 4 
                    for output every 6 hours in a gregorian calendar)
    gcm_data_freq:  Either 'day' or 'month': Frequency of the input annual 
                    cycle, either daily or monthly averages.
    out_path:       path to a folder to put the output netcdf files.
    
    Output:
    Multiple netcdf files! One netcdf file per output timestep. They are 
    numbered upwards from 0, starting at the beginning of the year and 
    increasing by 1 for each time step. Format var_nameNNNNN.nc, where 
    NNNNN is the number of the timestep 
    (e.g. ta00432.nc; for timestep number 432 and a variable called ta).
    """

    #open the inputfile and variable
    infile = xr.open_dataset(gcm_file_path)[var_name].squeeze()

    #enumerate the time steps for easyier interpolation
    timesteps = len(infile['time'])
    infile['time'] = np.arange(timesteps)

    #create the new time steps if daily inputdata is given 
    # (no special treatment necessary)
    if gcm_data_freq == 'day':
        tnew = np.linspace(0, timesteps-1, num=n_out_time_steps)

    #create new time steps with monthly data (shift the montly means to the 
    # end of the month and add a dublicate the first time step)
    if gcm_data_freq == 'month':
        tnew = np.linspace(0, timesteps, num=n_out_time_steps)

        jan = 0.5*(infile[0] + infile[-1])
        first = infile[:-1, ...]
        second = infile[1:, ...]
        first['time'] = second['time'] #metadata must match for xarray computation. Dangerous thing...

        rest = 0.5*(first + second)
        end = jan.copy()

        jan['time'] = rest['time'][0] - 1
        end['time'] = rest['time'][-1] + 1

        infile = xr.concat([jan,rest,end], dim='time').transpose('time', ...).chunk({'time':13})


    #if len(infile.shape) > 3:
            #infile = infile.chunk({'plev':1})	
    #interpolate new output timesteps
    outfile = infile.interp(time=tnew, method='linear', assume_sorted=True)#.chunk({'time':1})
    infile.close()
    del infile
    gc.collect()

    #numerate n_out_time_steps
    outfile['time'] = np.arange(n_out_time_steps)

    #save each output timestep in a seperate netcdf file for easyier handling later on
    Path(out_path).mkdir(parents=True, exist_ok=True)

    IMP = IterMP(njobs=n_par, run_async=True)
    fargs = {
            'outfile':outfile,
            'out_path':out_path,
            'var_name':var_name,
            }
    step_args = []
    for filenum in range(n_out_time_steps):
        step_args.append({'filenum':filenum})
    IMP.run(save_time_step, fargs, step_args)

    outfile.close()
    del outfile
    gc.collect()



def interpannualcycle_dask(gcm_file_path, var_name, n_out_time_steps, gcm_data_freq, out_path, threads_per_task):
	"""
	Dublication of function above but using dask arrays instead of numpy arrays.
	If memory problems occur with the original function, this can be used instead.
	
	An additional setting is the amount of threads per task, which will be used for paralell saving of 
	netcdf files.
	"""

	#open the inputfile and variable
	infile = xr.open_dataset(gcm_file_path, chunks={'lat':1})[var_name].squeeze()

	#enumerate the time steps for easyier interpolation
	timesteps = len(infile['time'])
	infile['time'] = np.arange(timesteps)

	#create the new time steps if daily inputdata is given (no special treatment necessary)
	if gcm_data_freq == 'day':
		tnew = np.linspace(0, timesteps-1, num=n_out_time_steps)

#create new time steps with monthly data (shift the montly means to the end of the month and add a dublicate the first time step)
	if gcm_data_freq == 'month':
		tnew = np.linspace(0, timesteps, num=n_out_time_steps)

		jan = 0.5*(infile[0] + infile[-1])
		first = infile[:-1, ...]
		second = infile[1:, ...]
		first['time'] = second['time'] #metadata must match for xarray computation. Dangerous thing...
    
		rest = 0.5*(first + second)
		end = jan.copy()
    
		jan['time'] = rest['time'][0] - 1
		end['time'] = rest['time'][-1] + 1
    
		infile = xr.concat([jan,rest,end], dim='time').transpose('time', ...).chunk({'time':13})

	outfile = infile.interp(time=tnew, method='linear', assume_sorted=True)
	infile.close()

	#numerate n_out_time_steps
	outfile['time'] = np.arange(n_out_time_steps)
	outfile = outfile.chunk({'time':1, 'lat':-1})
	
	#save each output timestep in a seperate netcdf file for easyier handling later on
	Path(out_path).mkdir(parents=True, exist_ok=True)
	
	#within processes allow paralell io
	threadsteps = threads_per_task
	for large_file_num in range(0,n_out_time_steps,threadsteps):
		compute_list = []
		for small_file_num in range(large_file_num,large_file_num+threadsteps,1):
			compute_list.append(outfile[small_file_num].to_netcdf(f"{out_path}/{var_name}{small_file_num:05d}.nc", mode='w', compute=False))
		dask.compute(*compute_list, scheduler='threads')
                
if __name__ == "__main__":
	gcm_file_path = str(sys.argv[1])
	var_name = str(sys.argv[2])
	n_out_time_steps = int(sys.argv[3])
	gcm_data_freq = str(sys.argv[4])
	out_path = str(sys.argv[5])
	interpannualcycle(gcm_file_path, var_name, n_out_time_steps, gcm_data_freq, out_path)
