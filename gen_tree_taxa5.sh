#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
num=0
dir=/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_tall10/01
outfile=$dir/gen_net.nex
phylonet=/Users/zhen/Desktop/phylonet/build/jar/PhyloNet.jar
for numgt in 100
do
    for sitespergt in 500
    do
    for i in 26
#    for i in 1 2
      do
         mkdir $dir/$i 2> /dev/null
         mkdir $dir/$i/homo 2> /dev/null
         rm  $(dirname "$0")seedms 2> /dev/null
         rm $outfile 2> /dev/null
         printf "#NEXUS\n\nBEGIN PHYLONET;\n" >> $outfile
         printf "SimSeqInNetwork -theta 0.005 -base (0.2112,0.2888,0.2896,0.2104) -rate (0.2173,0.9798,0.2575,0.1038,1,0.2070) " >> $outfile
         printf " -numgt %d -sitespergt %d " $numgt $sitespergt >> $outfile
         printf ' -mspath "/Users/zhen/Desktop/Zhen/research/phylogenetics/treeAugment/proj/msdir/ms" ' >> $outfile
         printf ' -seqgenpath "/Users/zhen/Desktop/Zhen/research/phylogenetics/treeAugment/proj/Seq-Gen/source/seq-gen" ' >> $outfile
         printf ' -tm <C:C_0;G:G_0;L:L_0;Q:Q_0;R:R_0;Z:Z_0> -truenet "((((Q:1.5,R:1.5):2.5,L:4):1,(G:2,C:2):3):5,Z:10);" -out "' >> $outfile
         printf "$dir/$i/homo/markers.txt">> $outfile
         printf '"  -gtfile "' >> $outfile
         printf "$dir/$i/homo/genetrees.txt" >> $outfile
         printf '";\n' >> $outfile
         printf 'END;\n' >> $outfile
         java -jar $phylonet $outfile
         printf "Generated numgt: %d sitespergt: %d group: %d num: %d\n" $numgt $sitespergt $i $num
         num=$((num+1))
         cp $(dirname "$0")/seedms $dir/$i/

      done
    done
    cp $(dirname "$0")/seedms $dir/$i/
done


