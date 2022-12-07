import argparse
from netCDF4 import Dataset
from functions import load_delta


var_name_map = {
    'ts':   'T_CL',
}




def extpar_adapt(ext_file_path, delta_inp_path):

    ext_file = Dataset(ext_file_path, 'a')

    # update T_C
    print('update deep soil temperature')

    delta_ts = load_delta(delta_inp_path, 'ts', None)

    ## Make sure dimensions are exactly the same.
    ## There are numerical differences between CDO remapped objects
    ## and xarray data...
    #delta_ts = delta_ts.assign_coords(
    #                {'rlat':ext_file.rlat.values,
    #                 'rlon':ext_file.rlon.values})

    delta_ts_clim = delta_ts.mean(dim=['time'])
    print(delta_ts_clim)

    ext_file['T_CL'][:] += delta_ts_clim.values.squeeze()

    ext_file.close()
    print('Done.')








if __name__ == "__main__":

    ## input arguments
    parser = argparse.ArgumentParser(description =
                    'COSMO-specific: Perturb Extpar soil temperature climatology with ts climate delta.')
    # extpar file to modify
    parser.add_argument('extpar_file_path', type=str,
            help='Path to extpar file to modify T_CL.')

    # climate delta directory (already remapped to ERA5 grid)
    parser.add_argument('-d', '--delta_input_dir', type=str, default=None,
            help='Directory with GCM climate deltas to be used. ' +
            'This directory should have a climate delta for ts ' +
            'already horizontally remapped to the grid of ' +
            'the extpar file which can perhaps be done with ' +
            'step_02_preproc_deltas.py or otherwise with CDO.')
    args = parser.parse_args()
    print(args)

    extpar_adapt(args.extpar_file_path, args.delta_input_dir)

