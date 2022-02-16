#!/bin/bash

out_base_dir=/net/o3/hymet_nobackup/heimc/data/pgw
inp_base_dir=/net/atmos/data/cmip6

table_ID=Emon

model=MPI-ESM1-2-HR
experiments=(historical ssp585 delta)

var_names=(ua va ta zg hur hus)
#var_names=(hur)


table_ID_out_base_dir=$out_base_dir/$table_ID
out_dir=$table_ID_out_base_dir/$model

for var_name in ${var_names[@]}; do
    echo "#################################################################"
    echo $var_name
    echo "#################################################################"


    # add Amon model top values to Emon
    if [[ "$table_ID" == "Emon" ]]; then
        for experiment in ${experiments[@]}
        do
            echo $experiment

            #mv $out_dir/${var_name}_${experiment}.nc \
            #    $out_dir/Emon_model_bottom_${var_name}_${experiment}.nc
            cdo sellevel,100000,97500,95000,92500,90000,87500,85000,82500,80000,77500,75000,70000,65000,60000,55000,50000,45000,40000,35000,30000,25000,22500,20000,17500,15000,12500,10000 $out_dir/${var_name}_${experiment}.nc \
                $out_dir/Emon_model_bottom_${var_name}_${experiment}.nc

            Amon_out_base_dir=$out_base_dir/Amon
            Amon_out_dir=$Amon_out_base_dir/$model
            cdo sellevel,7000,5000,3000,2000,1000,500,100 \
                $Amon_out_dir/${var_name}_${experiment}.nc \
                $out_dir/Amon_model_top_${var_name}_${experiment}.nc
            cdo -O merge \
                $out_dir/Emon_model_bottom_${var_name}_${experiment}.nc \
                $out_dir/Amon_model_top_${var_name}_${experiment}.nc \
                $out_dir/${var_name}_${experiment}.nc
        done
    fi


done

