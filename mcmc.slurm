#!/bin/bash

#author:

#Zhen.Cao

#date:

#May 2022

#purpose:

#run data test



#SBATCH --job-name=mcmc

#SBATCH --account=commons

#SBATCH --partition=commons

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=24

#SBATCH --mem-per-cpu=4G

#SBATCH --time=24:00:00

#SBATCH --mail-user=zc36@rice.edu

#SBATCH --mail-type=ALL


module load Java/1.8.0_201
export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH
srun java -jar /scratch/zc36/X2/mcmc/PhyloNet.jar $1 $1/mcmc.nex 10000000 1000000 $2 5000 0 20 $2 null $3 true
