#!/bin/bash
#SBATCH --account=pmg
#SBATCH --job-name=trajectory
#SBATCH --output=logs/traj_%j.out
#SBATCH --error=logs/traj_%j.err
#SBATCH --time=10:00:00
#SBATCH --mem=240G
#SBATCH -c 8
#SBATCH --gres=gpu:1
#SBATCH --exclude=ins092,ins082

source $(conda info --base)/etc/profile.d/conda.sh
conda activate /insomnia001/depts/pmg/users/rc3517/conda_envs/sc_cellrank

# create logs directory if it doesn't exist
mkdir -p logs

# print job information
echo "Starting trajectory analysis (gpu) with cellrank"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "----------------------------------------"

# run the Python script
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
python trajectory.py

# print completion information
echo "----------------------------------------"
echo "End time: $(date)"
echo "Job completed"
