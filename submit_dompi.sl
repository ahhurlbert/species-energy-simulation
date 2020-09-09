#!/bin/bash

#SBATCH -p 528_queue
#SBATCH -N 2
#SBATCH --ntasks-per-node=44
#SBATCH --cpus-per-task=1
#SBATCH -t 3-
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=Hurlbert@bio.unc.edu

module add r/3.6.0

mpirun Rscript run_sim_on_cluster.r 6501 6504
