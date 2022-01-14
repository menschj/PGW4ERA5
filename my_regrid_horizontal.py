# -*- coding: utf-8 -*-
###############################################################################
import os
import xarray as xr
import numpy as np
#import xesmf as xe
from package.utilities import cdo_remap
from package.mp import IterMP
import sys
from pathlib import Path
###############################################################################

def run_remap_with_cdo_fix(grid_des_file, inp_file, out_file, method='bil'):
    # call cdo
    cdo_remap(grid_des_file, inp_file, out_file, method='bil')
    print('ATTENTION: fix cdo remapping issue is running!')
    # for strange reason cdo does not fill the rlon values.. (why???)
    # therefore fix this
    ds = xr.open_dataset(out_file)
    ds = ds.assign_coords({'rlon':np.unique(ds.lon)}) 
    ds.to_netcdf(out_file+'.tmp')
    ds.close()
    os.rename(out_file+'.tmp', out_file)

def regridhorizontal(infolder, variablename, n_out_time_steps,
                    out_grid_file, outputfolder, njobs):
    """
    Regrid the output of the interpolate function for one variable to a 
    different domain. Regridding will be performed using xesmf.

    Input:
    infolder:       folder where the data of the target variable is located 
                    (the output of the interpolate.py)
    variablename:   name of the variable to be regridded
    n_out_time_steps: 
                    how many timesteps need to be regridded 
                    (corresponds to the number of files)?
    out_grid_file:  a netcdf file that is in the target grid --> all input 
                    will be regridded to the grid defined in this netcdf file.
    outputfolder:   path to a folder to write output. 
                    Overwriting files seems to cause problems, 
                    so it is recommended to use a new folder.

    Output: One netcdf file per timestep for the selected variable regirdded 
            to the defined target grid.
    """

    Path(outputfolder).mkdir(parents=True, exist_ok=True)

    IMP = IterMP(njobs=njobs, run_async=True)
    fargs = {'grid_des_file':out_grid_file,}
    step_args = []
    for stepnum in range(n_out_time_steps):
    #for stepnum in range(1696,1696+11*8):
        #print(stepnum)
        infile = f"{infolder}/{variablename}{stepnum:05d}.nc"
        outfile = f"{outputfolder}/{variablename}{stepnum:05d}.nc"
        step_args.append({'inp_file':infile,
                          'out_file':outfile})
    IMP.run(run_remap_with_cdo_fix, fargs, step_args)

    #targetgrid = xr.open_dataset(out_grid_file)
    #infile = xr.open_dataset(f"{infolder}/{variablename}00000.nc")


    #regridder = xe.Regridder(infile, targetgrid, 'bilinear', reuse_weights=True, filename='/scratch/snx3000/heimc/'+variablename+'regridder.nc')
#regridder = xe.Regridder(infile, targetgrid, 'bilinear', reuse_weights=True)

    #outfile = regridder(infile)
    #cdo_remap(out_grid_file, infile, outfile, method='bil')
    #quit()
    #infile.close()

    #outfile.to_netcdf(f"{outputfolder}/{variablename}{stepnum:05d}.nc")
    #outfile.close()



if __name__ == "__main__":
    infolder=str(sys.argv[1])
    variablename=str(sys.argv[2])
    n_out_time_steps=int(sys.argv[3])
    out_grid_file=str(sys.argv[4])
    outputfolder=str(sys.argv[5])
    regridhorizontal(infolder, variablename, n_out_time_steps, out_grid_file, outputfolder)
