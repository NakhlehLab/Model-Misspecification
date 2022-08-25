#!/bin/bash

# Summarize MCMC_SEQ result
# pypath_mcmc=~/X2/code/mcmc_heter/MCMC_result_summarize.py
# dir_path=/home/zc36/X2/mcmc/net0/long

pypath_mcmc=/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/code/mcmc_heter/ML_result_summarize.py
dir_path=/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/net0/long/

rep=1
sub_dir=$dir_path/$rep/heter/500
burnin=2000000
samplefreq=5000



# python $pypath_mcmc $sub_dir/murate.txt $sub_dir $dir_path/$rep/homo/500/genetrees.txt $burnin $samplefreq
python $pypath_mcmc $sub_dir/murate.txt $sub_dir $dir_path/$rep/homo/500/genetrees.txt $burnin $samplefreq
