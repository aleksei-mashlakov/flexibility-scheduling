#!/bin/bash
# created: Nov 24, 2019 5:07 PM
# author: mashlakov

#!/bin/bash
#SBATCH --account=Project_2001505
#SBATCH --job-name=online_control
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=23:00:00
#SBATCH --output=job_out_test_abs.txt
#SBATCH --error=job_err_test_abs.txt
#SBATCH --mem-per-cpu=8G
##SBATCH --cpus-per-task=4

# set the number of threads based on --cpus-per-task
#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export MODIN_ENGINE=ray  # Modin will use Ray
#export MODIN_ENGINE=dask  # Modin will use Dask
module load python-data/3.7.6-1
#module load dask
#python -m
#pip install --user pandas
#pip3 install modin --user

#python3 -m pip install --upgrade pip --user
#python3 -m pip install pandas  --force-reinstall --user
#python3 -m pip install modin  --force-reinstall --user
#python3 -m pip install ray  --force-reinstall --user

#pip3 install ray --user
echo "Hola el patron!"

# Train and test the rest of the horizons
echo "running ..."
python3 battery_control_old.py
echo "done"
seff $SLURM_JOBID
