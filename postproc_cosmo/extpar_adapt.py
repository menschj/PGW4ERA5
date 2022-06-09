import argparse
from netCDF4 import Dataset
from functions import load_delta


var_name_map = {
    'tas':  'T_CL',
}




def extpar_adapt(ext_file_path, delta_inp_path):

    ext_file = Dataset(ext_file_path, 'a')

    # update T_C
    var_name = 'tas'
    print('update {}'.format(var_name))

    name_base='{}_delta.nc'
    delta_tas = load_delta(delta_inp_path, var_name, None)
    print(delta_tas)

    ## Make sure dimensions are exactly the same.
    ## There are numerical differences between CDO remapped objects
    ## and xarray data...
    #delta_tas = delta_tas.assign_coords(
    #                {'rlat':ext_file.rlat.values,
    #                 'rlon':ext_file.rlon.values})

    delta_tas_clim = delta_tas.mean(dim=['time'])

    ext_file['T_CL'][:] += delta_tas_clim.values.squeeze()

    ext_file.close()
    print('Done.')








if __name__ == "__main__":

    ## input arguments
    parser = argparse.ArgumentParser(description =
                    'COSMO-specific: Perturb Extpar soil temperature climatology with PGW climate delta.')
    # extpar file to modify
    parser.add_argument('extpar_file_path', type=str,
            help='Path to extpar file to modify T_CL.')

    # climate delta directory (already remapped to ERA5 grid)
    parser.add_argument('-d', '--delta_input_dir', type=str, default=None,
            help='Directory with GCM climate deltas to be used. ' +
            'This directory should have a climate delta for tas ' +
            'already horizontally remapped to the grid of ' +
            'the extpar file (see step_02_preproc_deltas.py).')
    args = parser.parse_args()
    print(args)

    extpar_adapt(args.extpar_file_path, args.delta_input_dir)

