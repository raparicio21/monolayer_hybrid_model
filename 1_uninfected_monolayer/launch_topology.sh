#!/bin/bash
#SBATCH --job-name=my_topology_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00

#module load matlab
#/opt/matlab/bin/matlab -nosplash -nodesktop -r "run('topology.m'); exit;"
/usr/local/matlab/R2021a/bin/matlab -nosplash -nodesktop -r "run('topology.m'); exit;"

