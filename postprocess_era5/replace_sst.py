import sys, pathlib, os, time, shutil
from netCDF4 import Dataset
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt

i_do_replace = 0

delta_hour_inc = 3
timestep = timedelta(hours=delta_hour_inc)
first_dt = datetime(2006,8,30,18)
last_dt = datetime(2006,9,1,0)
lbfd_dts = np.arange(first_dt, last_dt+timestep, timestep).tolist()
print('Run for dates {}'.format(lbfd_dts))

sst_src_dir = os.path.join('/scratch','snx3000','heimc','lmp',
                           'wd','94080100_SA_3_long','int2lm_out_PGW_full')
year_shift_sst_src = 88
targ_lbfd_dir = os.path.join('/scratch','snx3000','heimc','lmp',
                           'wd','06080100_SA_3_long','int2lm_out_PGW_SST')
const_file_path = os.path.join('/project','pr94','heimc','data',
                           #'cosmo_out','SA_3_long','lm_c', 'FR_LAND.nc')
                           'cosmo_out','SA_3_long','lm_c', 'lffd20060801000000c.nc')

#print(os.listdir(const_file_path))
#quit()

for lbfd_dt in lbfd_dts:
    print(lbfd_dt)

    targ_lbfd_path = os.path.join(targ_lbfd_dir, 
                        'lbfd{:%Y%m%d%H%M%S}.nc'.format(lbfd_dt))
    sst_dt = lbfd_dt.replace(year=lbfd_dt.year + year_shift_sst_src)
    sst_src_path = os.path.join(sst_src_dir, 
                        'lbfd{:%Y%m%d%H%M%S}.nc'.format(sst_dt))
    #print(targ_lbfd_path)
    #print(sst_src_path)

    ## replace W_SO in TMP_laf based on W_SO in src_file
    tar_nc = Dataset(targ_lbfd_path, 'a')
    src_nc = Dataset(sst_src_path, 'r')
    const_nc = Dataset(const_file_path, 'r')

    fr_land = const_nc['FR_LAND'][:]
    hsurf = const_nc['HSURF'][:]
    t_s = src_nc['T_S'][:]
    tar_t_s = tar_nc['T_S'][:]
    tar_t_s[(fr_land == 0) & (hsurf == 0)] = t_s[(fr_land == 0) & (hsurf == 0)]
    tar_nc['T_S'][:] = tar_t_s

    const_nc.close()
    src_nc.close()
    tar_nc.close()

