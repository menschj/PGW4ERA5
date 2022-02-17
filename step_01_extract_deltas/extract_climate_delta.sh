#!/bin/bash
##############################################################################
# Template script to extract climate deltas from CMIP GCMs.
# Computes climatologies and the climate deltas for specific CMIP 
# experiments and members.
# The script only serves as inspiration for the extraction of the climate
# deltas and it is not generally valid and has to be adjusted 
# for different GCMs, experiments, table_IDs and climatology periods.

##### IMPORTANT NOTES FOR Emon DATA USERS:
## FIRST
## hur is not available in the Emon output and has to be computed from hus.
## Since this is not excat because it matters if hur is computed from hus
## on the model levels every time the GCM writes an output or whether it is
## computed afterwards. Therefore, the hur is used from the coarse
## resolution Amon output but vertically interpolated to the higher resolution
## of Emon using information from the hur computed on using Emon hus and ta
## data. This helps to do a "better informed" vertical interpolation of the
## coarse Amon hur to the higher resolved Emon grid.
## This steps are done by the 
## SECOND
## After computing hur for Emon, the delta will have to still be computed
## by running the script again only for var_names=(hur).
## THIRD
## The Emon data does not reach up as far as the Amon data
## Thus, after completing this script, run add_Emon_model_top.sh
## to add the model top data from Amon to the Emon fields.
## FOURTH
## These Emon fixes are still a bit experimental and not very user-friendly.
## Sorry for that!
##############################################################################

# USER SETTINGS
##############################################################################
# base directory where cmip6 data is stored
cmip_data_dir=/net/atmos/data/cmip6
# base directory where output should be stored
out_base_dir=/net/o3/hymet_nobackup/heimc/data/pgw/deltas/native

# name of the GCM to extract data for
gcm_name=MPI-ESM1-2-HR

## CMIP experiments to use to compute climate deltas
## --> climate delta = future climatology - ERA climatology
# CMIP experiment to use for ERA climatology 
era_climate_experiment=historical
# CMIP experiment to use for future climatology 
future_climate_experiment=ssp585


## type of CMIP6 model output (e.g. monthly or daily, etc.)
## to use
# standard monthly output
table_ID=Amon
## high-resolution monthly data for only very few GCMs
#table_ID=Emon
## standard daily output
#table_ID=day


## select variables to extract
if [[ "$table_ID" == "Amon" ]]; then
    var_names=(tas hurs ps ua va ta hur zg)
elif [[ "$table_ID" == "day" ]]; then
    var_names=(tas hurs ps ua va ta hur zg)
elif [[ "$table_ID" == "Emon" ]]; then
    var_names=(ua va ta hus zg)
fi
## for Emon, compute hur delta separately after 
## it has been derived from hus
## (turn off i_exctract_vars but keep i_compute_delta)
#var_names=(hur)


# should variables be exctracted for the two climatologies=
i_extract_vars=1
# should climate deltas be computed?
i_compute_delta=1
# for Emon, hur is not available and needs to be approximated
# using the high-resolution hus climatology, as well as the hur climatology
if [[ "$table_ID" == "Emon" ]]; then
    i_convert_hus_to_hur=1
else
    i_convert_hus_to_hur=0
fi

## subdomain for which to extract GCM data
## should be either global (-180,180,-90,90)
## or anything larger than ERA5 subdomain
## except for storage and performance reasons, there is no benefit of
## using a subdomain.
box=-180,180,-90,90
# subdomain
#box=-74,40,-45,35
#box=-90,40,-45,50

# select appropriate cdo time aggregation command
# depending if input data is monthly or daily.
if [[ "$table_ID" == "day" ]]; then
    cdo_agg_command=ydaymean
else 
    cdo_agg_command=ymonmean
fi

# iterate over both experiments to extract data
experiments=($era_climate_experiment $future_climate_experiment)

##############################################################################

out_dir=$out_base_dir/$table_ID/$gcm_name
mkdir -p $out_dir

for var_name in ${var_names[@]}; do
    echo "#################################################################"
    echo $var_name
    echo "#################################################################"

    for experiment in ${experiments[@]}; do
        echo "#######################################"
        echo $experiment
        echo "#######################################"

        if [[ $i_extract_vars == 1 ]]; then

            # data folder hierarchy for CMIP6
            inp_dir=$cmip_data_dir/$experiment/$table_ID/$var_name/$gcm_name/r1i1p1f1/gn

            # start of the CMIP6 file names
            file_name_base=${table_ID}_${gcm_name}_${experiment}_r1i1p1f1_gn

            ## overwrite old data
            #rm $out_dir/${var_name}_${experiment}.nc

            ## compute ERA climatology
            if [[ "$experiment" == "$era_climate_experiment" ]]; then
                cdo $cdo_agg_command \
                    -sellonlatbox,$box \
                    -cat \
                    $inp_dir/${var_name}_${file_name_base}_198[5-9]*.nc \
                    $inp_dir/${var_name}_${file_name_base}_199*.nc \
                    $inp_dir/${var_name}_${file_name_base}_200*.nc \
                    $inp_dir/${var_name}_${file_name_base}_201[0-4]*.nc \
                    $out_dir/${var_name}_${experiment}.nc

            ## compute future experiment climatology
            elif [[ "$experiment" == "$future_climate_experiment" ]]; then
                cdo $cdo_agg_command \
                    -sellonlatbox,$box \
                    -cat \
                    $inp_dir/${var_name}_${file_name_base}_20[7-9]*.nc \
                    $out_dir/${var_name}_${experiment}.nc
            fi
        fi

        ## convert hus to hur if required
        if [[ $i_convert_hus_to_hur == 1 ]]; then
            if [[ "$var_name" == "hus" ]]; then
                echo Convert Emon hus to hur using Amon hur data.

                Amon_out_dir=$out_base_dir/Amon/$gcm_name

                python convert_Emon_hus_to_hur.py  \
                    $out_dir/hus_${experiment}.nc \
                    $out_dir/ta_${experiment}.nc \
                    $out_dir/hur_${experiment}.nc \
                    -a $Amon_out_dir/hur_${experiment}.nc
            fi
        fi

    done

    ## compute delta (future climatology - ERA climatology)
    if [[ $i_compute_delta == 1 ]]; then
        cdo sub $out_dir/${var_name}_$future_climate_experiment.nc \
                $out_dir/${var_name}_$era_climate_experiment.nc \
                $out_dir/${var_name}_delta.nc
    fi

done


# link surface fields from Amon because they are not available
# in Emon
if [[ "$table_ID" == "Emon" ]]; then
    cd $out_dir
    ln -s $out_base_dir/Amon/$gcm_name/*s_*.nc .
    cd -
fi

