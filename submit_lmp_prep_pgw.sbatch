#!/bin/bash -l
#
#SBATCH --job-name="lmp pgw"
#SBATCH --time=23:30:00
#SBATCH --cpus-per-task=12
#SBATCH --partition=normal
##SBATCH --output=bOut/%jout
##SBATCH --error=bOut/%jerr
#SBATCH --output=/scratch/snx3000/heimc/lmp/out_prep_pgw/%jout0
#SBATCH --error=/scratch/snx3000/heimc/lmp/out_prep_pgw/%jerr0
#SBATCH --ntasks=1
#SBATCH --account='pr94'
#SBATCH --constraint=gpu

source ~/phd_scripts/load_modules_daint.sh
echo -e "$SLURM_JOB_NAME started on $(date):\n"

## rename 2D fields
#srun -c $SLURM_CPUS_PER_TASK python -u my_rename_variables.py tas
#srun -c $SLURM_CPUS_PER_TASK python -u my_rename_variables.py hurs
#srun -c $SLURM_CPUS_PER_TASK python -u my_rename_variables.py hurs,tas


### MANUALLY
############################################################
############################################################

### perpare laf file
############################################################
#srun -c $SLURM_CPUS_PER_TASK python -u laf_adapt.py

### perpare lbfd files
############################################################

#year=2006
#month=09
#start_day=01
#end_day=30
#
#start_dayhr=${start_day}00
#end_dayhr=${end_day}21
#sim_start_date=${year}${month}01
#
#arg_f=2006090100
#arg_l=2006090100
#
#arg_f=${year}${month}${start_dayhr}
#arg_l=${year}${month}${end_dayhr}
#srun -c $SLURM_CPUS_PER_TASK python -u my_lbfd_adapt.py -d 3 \
#                -f ${arg_f} -l ${arg_l} -s ${sim_start_date}





### AUTOMATICALLY
############################################################
############################################################
sim_start_date=${1}
arg_f=${2}
arg_l=${3}00
busy_file=${4}
sim_name_base=${5}

echo $sim_start_date
echo $arg_f
echo $arg_l

##TODO debug
#arg_l=2006080200

### perpare laf file
############################################################
echo RUN laf_adapt.py
echo "#######################################################"
srun -c $SLURM_CPUS_PER_TASK python -u laf_adapt.py -d 3 \
                -s ${sim_start_date} -n ${sim_name_base}


### perpare lbfd files
############################################################
echo RUN lbfd_adapt.py
echo "#######################################################"
srun -c $SLURM_CPUS_PER_TASK python -u my_lbfd_adapt.py -d 3 \
                -f ${arg_f} -l ${arg_l} -s ${sim_start_date} \
                -n ${sim_name_base}


echo done! delete busy file: $busy_file
rm $busy_file



echo -e "$SLURM_JOB_NAME finished on $(date)\n"


