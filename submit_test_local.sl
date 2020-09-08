#!/bin/bash

#SBATCH -p debug_queue
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -t 01:00:00
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=Hurlbert@bio.unc.edu

module add r/3.6.0

mpirun Rscript test_fun_on_cluster_nonparallel.r 1 2
