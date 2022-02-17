#!/bin/bash
##############################################################################
# Template script to extract climate deltas from CMIP GCMs.
# Computes climatologies for specific CMIP experiments and members
##############################################################################

# USER SETTINGS
##############################################################################
delta_inp_base_dir=/net/atmos/data/cmip6
out_base_dir=/net/o3/hymet_nobackup/heimc/data/pgw

gcm_name=MPI-ESM1-2-HR

table_ID=Amon
#table_ID=Emon
#table_ID=day

cdo_agg_command=ymonmean
cdo_agg_command=ydaymean

model=MPI-ESM1-2-HR
experiments=(historical ssp585)

i_extract_vars=1
i_convert_hus_to_hur=0
i_compute_delta=1

#subsel_command=-sellonlatbox,-74,40,-45,35
box=-180,180,-90,90

# Amon
var_names=(tas hurs ps ua va ta hur zg)
#var_names=(tas)
### Emon
#var_names=(ua va ta hus zg)
#var_names=(hus)
## compute delta separately
#var_names=(hur)

##############################################################################


table_ID_out_base_dir=$out_base_dir/$table_ID

for var_name in ${var_names[@]}; do
    echo "#################################################################"
    echo $var_name
    echo "#################################################################"

    for experiment in ${experiments[@]}; do
        echo "#######################################"
        echo $experiment
        echo "#######################################"

        out_dir=$table_ID_out_base_dir/$model
        mkdir -p $out_dir

        if [[ $i_extract_vars == 1 ]]; then
            inp_dir=$delta_inp_base_dir/$experiment/$table_ID/$var_name/$model/r1i1p1f1/gn
            file_name_base=${table_ID}_${model}_${experiment}_r1i1p1f1_gn
            rm $out_dir/plev_${var_name}_${experiment}.nc

            ## compute historical climatology
            if [[ "$experiment" == "historical" ]]; then
                cdo $cdo_agg_command \
                    -sellonlatbox,$box \
                    -cat \
                    $inp_dir/${var_name}_${file_name_base}_198[5-9]*.nc \
                    $inp_dir/${var_name}_${file_name_base}_199*.nc \
                    $inp_dir/${var_name}_${file_name_base}_200*.nc \
                    $inp_dir/${var_name}_${file_name_base}_201[0-4]*.nc \
                    $out_dir/${var_name}_${experiment}.nc

            ## compute future experiment (here ssp585) climatology
            elif [[ "$experiment" == "ssp585" ]]; then
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

                Amon_out_dir=$out_base_dir/Amon/$model

                python convert_hus_to_hur.py  \
                    $out_dir/hus_${experiment}.nc \
                    $out_dir/ta_${experiment}.nc \
                    $out_dir/hur_${experiment}.nc \
                    -a $Amon_out_dir/hur_${experiment}.nc
            fi
        fi

    done

    ## compute delta (ssp585 - historical)
    if [[ $i_compute_delta == 1 ]]; then
        cdo sub $out_dir/${var_name}_ssp585.nc \
                $out_dir/${var_name}_historical.nc \
                $out_dir/${var_name}_delta.nc
    fi

done
