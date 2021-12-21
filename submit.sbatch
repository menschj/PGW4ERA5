#!/bin/bash -l
#
#SBATCH --job-name="prep pgw"
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=12
#SBATCH --partition=normal
#SBATCH --output=/scratch/snx3000/heimc/pgw/out_scripts/out0
#SBATCH --error=/scratch/snx3000/heimc/pgw/out_scripts/err0
#SBATCH --ntasks=1
#SBATCH --account='pr94'
#SBATCH --constraint=gpu

### tasks
# time interpolation: 5
# regrid horizontal: 12

source ~/phd_scripts/load_modules_daint.sh
echo -e "$SLURM_JOB_NAME started on $(date):\n $command $1 $2\n"

var_names=hur,hurs,ta,tas,ua,va
#var_names=hur
var_names=ta
#var_names=ua
#var_names=va
#var_names=tas
#var_names=hurs


echo $var_names
srun -c $SLURM_CPUS_PER_TASK python -u my_settings.py \
                -p $SLURM_CPUS_PER_TASK $var_names -n 8

echo -e "$SLURM_JOB_NAME finished on $(date)\n"


