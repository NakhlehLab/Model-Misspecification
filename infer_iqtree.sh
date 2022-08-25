#!/bin/bash

srcdir=$1
dstdir=$2
num_taxa=$3
NUM_LOCI=$4

mkdir $dstdir

src_gt=$srcdir/genetrees.txt
dst_gt=$dstdir/genetrees.txt

locus_length=$7

#taxon_map="<A:A_0;G:G_0;L:L_0;Q:Q_0>"

base_locus_rate=$5
outgroup=$6
LOCUS_HETEROGENEITY=True
dst_iq_gt=$dstdir/$locus_length/rooted_iqtree.txt

# 3. estimate IQTREE
 iqscriptpath=/home/zc36/X2/code/prepareIQTREE.py
 pruneOutgroup=False
 RUN_IQ=True
mkdir $dstdir/$locus_length/IQTREE
python $iqscriptpath $dstdir/$locus_length/ $pruneOutgroup $RUN_IQ $outgroup $num_taxa $locus_length False $NUM_LOCI

# # 4. compute IQTREE error profile
 proflingPath=/home/zc36/X2/code/gt_error_profling.py
 python $proflingPath $src_gt $dstdir/$locus_length/murate.txt $dst_iq_gt $outgroup $base_locus_rate 1>$dstdir/$locus_length/gt_error.txt

rm -r $dstdir/$locus_length/IQTREE
