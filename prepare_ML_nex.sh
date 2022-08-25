#!/bin/bash

nexpath=/home/zc36/X2/code/seq_to_nexus.py

x=0
num_taxa=6
# gtpath=/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5/1/heter/rooted_iqtree.txt
# outpath=/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5/1/heter/ML_$x.nex
#taxon_map="<A:A_0;L:L_0;Q:Q_0;Z:Z_0>"
taxon_map="<C:C_0;G:G_0;L:L_0;R:R_0;Q:Q_0;Z:Z_0>"
#dirpath=/home/zc36/X2/taxa5_tall10/036/
dir=/home/zc36/X2/test_outgroup4/net0/
re_start=1
re_end=100
for scale in short medium long 
do
	dirpath=$dir/$scale/
	for((replicate=$re_start;replicate<=$re_end;++replicate))
	do
		# for locus_length in 500 1000
		for locus_length in 500 2000
		do
			gtpath=$dirpath/$replicate/heter/$locus_length/rooted_iqtree_resolved.txt
			outpath=$dirpath/$replicate/heter/$locus_length/ML_iq_$x.nex
			python $nexpath 2 $num_taxa $locus_length $gtpath $outpath $x $taxon_map
			
			gtpath=$dirpath/$replicate/homo/$locus_length/genetrees.txt
			outpath=$dirpath/$replicate/heter/$locus_length/ML_true_$x.nex
			python $nexpath 2 $num_taxa $locus_length $gtpath $outpath $x $taxon_map
		done
	done
#gtpath=$dirpath/$replicate/heter/rooted_iqtree.txt
#outpath=$dirpath/$replicate/heter/ML_iq_$x.nex
#python $nexpath 2 $num_taxa $locus_length $gtpath $outpath $x $taxon_map


done