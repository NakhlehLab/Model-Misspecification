#!/bin/bash

srcdir=$1
dstdir=$2
num_taxa=$3
NUM_LOCI=$4

mkdir $dstdir

src_gt=$srcdir/genetrees.txt
dst_gt=$dstdir/genetrees.txt

locus_length=500

#taxon_map="<A:A_0;G:G_0;L:L_0;Q:Q_0>"

dst_iq_gt=$dstdir/rooted_iqtree.txt

base_locus_rate=$5
outgroup=$6
LOCUS_HETEROGENEITY=True

# 1. generate heterogenenous gene trees and markers
generatorpath=simulation/model.py
# python $generatorpath $src_gt $dst_gt $heter_marker $dstdir/murate.txt $LINEAGE_HETEROGENEITY $LOCUS_HETEROGENEITY $NUM_LOCI $num_taxa $locus_length $base_locus_rate $RUN_SEQ_GEN

# 1.1 generate per gene lineage heterogeneity gene trees
# python $generatorpath True False $src_gt $NUM_LOCI $dst_gt $num_taxa

# 1.2 produce sequences with seq-gen given the gene trees

rate_dir=$(echo $base_locus_rate  | sed 's/.*\.//')
mkdir $dstdir/$rate_dir
heter_marker=$dstdir/$rate_dir/markers.txt
python $generatorpath False True $dst_gt $NUM_LOCI $LOCUS_HETEROGENEITY $heter_marker $dstdir/$rate_dir/murate.txt $locus_length $base_locus_rate


# # 3. estimate IQTREE
# iqscriptpath=/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/data/simulation/seqgen/code/prepareIQTREE.py
# pruneOutgroup=False
# RUN_IQ=True
# mkdir $dstdir/IQTREE

# python $iqscriptpath $dstdir $pruneOutgroup $RUN_IQ $outgroup $num_taxa $locus_length False

# # 4. compute IQTREE error profile
# proflingPath=/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/data/simulation/seqgen/code/gt_error_profling.py
# python $proflingPath $dst_gt $dstdir/murate.txt $dst_iq_gt $outgroup 1>$dstdir/gt_error.txt
