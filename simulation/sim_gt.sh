#!/bin/bash

srcdir=$1
dstdir=$2
num_taxa=$3
NUM_LOCI=$4

mkdir $dstdir

locus_length=$7

src_gt=$srcdir/$locus_length/genetrees.txt
dst_gt=$dstdir/$locus_length/genetrees.txt


#taxon_map="<A:A_0;G:G_0;L:L_0;Q:Q_0>"

dst_iq_gt=$dstdir/rooted_iqtree.txt

base_locus_rate=$5
outgroup=$6
LOCUS_HETEROGENEITY=True

# 1. generate heterogenenous gene trees and markers
generatorpath=simulation/model.py

# 1.2 produce sequences with seq-gen given the gene trees

rate_dir=$(echo $base_locus_rate  | sed 's/.*\.//')
mkdir $dstdir/$locus_length
heter_marker=$dstdir/$locus_length/markers.txt
python $generatorpath True $src_gt $NUM_LOCI $LOCUS_HETEROGENEITY $heter_marker $dstdir/$locus_length/murate.txt $locus_length $base_locus_rate



