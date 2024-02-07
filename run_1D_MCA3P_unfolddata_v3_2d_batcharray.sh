#!/bin/bash
#SBATCH -p standard
#SBATCH --account=t3
#SBATCH --job-name=job-name
#SBATCH --mem=30G                     # memory 3GB (per job)
#SBATCH --time 12:00:00
#SBATCH --array=0-6
#SBATCH -n 10
obs=("spherocity" "thrust" "transverse_spherocity" "transverse_thrust" "broadening" "isotropy" "mass")

python evaluate_iterations.py --config config/dataset_v7_MCA3P_unfolddata_v3_2d.json   --method multifold --migiter 0 1 --obs nparticle:nparticle_dim1_2D,${obs[${SLURM_ARRAY_TASK_ID}]} --eff-acc
