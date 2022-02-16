import argparse
from netCDF4 import Dataset
from base.functions import load_delta


var_name_map = {
    'tas':  'T_CL',
}




def extpar_adapt(ext_file_path, delta_inp_path):

    ext_file = Dataset(ext_file_path, 'a')

    # update T_C
    var_name = 'tas'
    print('update {}'.format(var_name))

    name_base='{}_delta.nc'
    delta_tas = load_delta(delta_inp_path, var_name, None, 
                            None, None)

    ## TODO: debug for too small domain 
    #delta_tas = delta_tas.where(~np.isnan(delta_tas), 0)

    ## Make sure dimensions are exactly the same.
    ## There are numerical differences between CDO remapped objects
    ## and xarray data...
    #delta_tas = delta_tas.assign_coords(
    #                {'rlat':ext_file.rlat.values,
    #                 'rlon':ext_file.rlon.values})

    delta_tas_clim = delta_tas.mean(dim=['time'])

    ext_file['T_CL'][:] += delta_tas_clim.values.squeeze()

    ext_file.close()








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

