#!/bin/bash -l
#
#SBATCH --job-name="post pgw4"
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=12
#SBATCH --partition=normal
#SBATCH --output=/scratch/snx3000/heimc/pgw/out_scripts/out2
#SBATCH --error=/scratch/snx3000/heimc/pgw/out_scripts/err2
#SBATCH --ntasks=1
#SBATCH --account='pr94'
#SBATCH --constraint=gpu

source ~/phd_scripts/load_modules_daint.sh
echo -e "$SLURM_JOB_NAME started on $(date):\n"

## rename 2D fields
## NOT NECESSARY ANYMORE!
#srun -c $SLURM_CPUS_PER_TASK python -u my_rename_variables.py tas
#srun -c $SLURM_CPUS_PER_TASK python -u my_rename_variables.py hurs
#srun -c $SLURM_CPUS_PER_TASK python -u my_rename_variables.py hurs,tas


### MANUALLY
############################################################
############################################################

### perpare laf file
############################################################
srun -c $SLURM_CPUS_PER_TASK python -u laf_adapt.py

### perpare lbfd files
############################################################
year=2006
start_month=08
end_month=09
start_day=01
end_day=01

start_dayhr=${start_day}00
end_dayhr=${end_day}00
#sim_start_date=${year}${start_month}${start_day}
sim_start_date=${year}${start_month}01 # CH202112

#arg_f=2006080100
#arg_l=2006091100

arg_f=${year}${start_month}${start_dayhr}
arg_l=${year}${end_month}${end_dayhr}

srun -c $SLURM_CPUS_PER_TASK python -u my_lbfd_adapt.py -d 3 \
                -f ${arg_f} -l ${arg_l} -s ${sim_start_date}



echo -e "$SLURM_JOB_NAME finished on $(date)\n"


