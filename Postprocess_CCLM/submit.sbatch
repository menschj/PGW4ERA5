#!/bin/bash -l
#
#SBATCH --job-name="post pgw"
#SBATCH --time=23:30:00
#SBATCH --cpus-per-task=12
#SBATCH --partition=normal
#SBATCH --output=bOut/out8
#SBATCH --error=bOut/err8
#SBATCH --ntasks=1
#SBATCH --account='pr94'
#SBATCH --constraint=gpu

source ~/phd_scripts/load_modules_daint.sh
echo -e "$SLURM_JOB_NAME started on $(date):\n $command $1 $2\n"

## rename 2D fields
#srun -c $SLURM_CPUS_PER_TASK python -u my_rename_variables.py tas
#srun -c $SLURM_CPUS_PER_TASK python -u my_rename_variables.py hurs
#srun -c $SLURM_CPUS_PER_TASK python -u my_rename_variables.py hurs,tas

##### perpare laf file
#srun -c $SLURM_CPUS_PER_TASK python -u laf_adapt.py

year=2006
month=08
start_day=30
end_day=31

start_dayhr=${start_day}00
end_dayhr=${end_day}21
sim_start_date=${year}${month}01

arg_f=${year}${month}${start_dayhr}
arg_l=${year}${month}${end_dayhr}

#arg_f=2006090100
#arg_l=2006090100


### perpare lbfd files
srun -c $SLURM_CPUS_PER_TASK python -u my_lbfd_adapt.py -d 3 \
                -f ${arg_f} -l ${arg_l} -s ${sim_start_date}

echo -e "$SLURM_JOB_NAME finished on $(date)\n"


