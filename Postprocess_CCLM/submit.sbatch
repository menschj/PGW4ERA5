#!/bin/bash -l
#
#SBATCH --job-name="post pgw"
#SBATCH --time=23:30:00
#SBATCH --cpus-per-task=12
#SBATCH --partition=normal
#SBATCH --output=bOut/out0
#SBATCH --error=bOut/err0
#SBATCH --ntasks=1
#SBATCH --account='pr94'
#SBATCH --constraint=gpu

source ~/phd_scripts/load_modules_daint.sh
echo -e "$SLURM_JOB_NAME started on $(date):\n $command $1 $2\n"

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
arg_f=${1}00
arg_l=${2}00
sim_start_date=${1}

echo $arg_f
echo $arg_l
echo $sim_start_date

### perpare laf file
############################################################
echo RUN laf_adapt.py
echo "#######################################################"
srun -c $SLURM_CPUS_PER_TASK python -u laf_adapt.py -d 3 \
                -s ${sim_start_date}

### perpare lbfd files
############################################################
echo RUN lbfd_adapt.py
echo "#######################################################"
srun -c $SLURM_CPUS_PER_TASK python -u my_lbfd_adapt.py -d 3 \
                -f ${arg_f} -l ${arg_l} -s ${sim_start_date}






echo -e "$SLURM_JOB_NAME finished on $(date)\n"


