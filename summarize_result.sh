#!/bin/bash

dir_path=/home/zc36/X2/taxa5/net0/
num_replicate=100
null_reti_num=0

# summary statistics
/usr/bin/time -v python  collapse_short_gt_branch.py $dir_path $num_replicate short 2>$dir_path/medium/collapse_gt_polytomy.err
/usr/bin/time -v python  collapse_short_gt_branch.py $dir_path $num_replicate medium 2>$dir_path/medium/collapse_gt_polytomy.err
/usr/bin/time -v python  collapse_short_gt_branch.py $dir_path $num_replicate long   2>$dir_path/long/collapse_gt_polytomy.err

# full network
/usr/bin/time -v java -jar PhyloNet.jar $dir_path $num_replicate $null_reti_num 2>$dir_path/phylonet_exp.err 1>phylonet_exp.out
/usr/bin/time -v python summarize_probs.py $dir_path $num_replicate $null_reti_num 2>$dir_path/summarize_probs.err
/usr/bin/time -v python summarize_stats_5taxa.py $dir_path $num_replicate $null_reti_num 2>$dir_path/summarize_stats.err
/usr/bin/time -v python summarize_plot.py $dir_path/statistics_$null_reti_num.csv $num_replicate 2>$dir_path/summarize_plot.err

# triplet
/usr/bin/time -v python summarize_triplet_probs.py $dir_path $num_replicate False $null_reti_num 2>$dir_path/summarize_triplet_probs.err
/usr/bin/time -v python summarize_triplet_stats.py $dir_path $num_replicate False $null_reti_num 2>$dir_path/summarize_triplet_stats.err
/usr/bin/time -v python summarize_triplet_multitest.py  $dir_path $num_replicate False $null_reti_num 2>$dir_path/summarize_triplet_multitest.err
