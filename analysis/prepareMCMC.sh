#!/bin/bash
nexpath=/home/zc36/X2/code/seq_to_nexus.py

dir_path=/home/zc36/X2/test_outgroup4/net0/
re_start=1
re_end=1
num_taxa=6
taxon_map="<C:C_0;L:L_0;R:R_0;Q:Q_0>"
seed=1
for scale in long
do
	for((replicate=$re_start;replicate<=$re_end;++replicate))
	do
		for locus_length in 500 2000
		do  
            homo_dir=$dir_path/$scale/$replicate/homo/$locus_length/
            heter_dir=$dir_path/$scale/$replicate/heter/$locus_length/
            mkdir $homo_dir/murate
            mkdir $heter_dir/murate

            echo $homo_dir/markers.txt
            python $nexpath 1 $num_taxa $locus_length $homo_dir/markers.txt $homo_dir/murate/mcmc.nex $taxon_map $seed True False
            python $nexpath 1 $num_taxa $locus_length $heter_dir/markers.txt $heter_dir/murate/mcmc.nex $taxon_map $seed True False
            python $nexpath 1 $num_taxa $locus_length $homo_dir/markers.txt $homo_dir/mcmc.nex $taxon_map $seed False False
            python $nexpath 1 $num_taxa $locus_length $heter_dir/markers.txt $heter_dir/mcmc.nex $taxon_map $seed False False

        done
	done
done