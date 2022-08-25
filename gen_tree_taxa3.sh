#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
num=0
dir=/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa3_tall/001/
outfile=$dir/gen_net.nex
for numgt in 100
do
    for sitespergt in 500
    do
    for i in 54 55 60 61 62 63 71 81 82 83
#    for i in 1 2
      do
         echo $i
         mkdir $dir/$i 2> /dev/null
         mkdir $dir/$i/homo 2> /dev/null

         rm  $(dirname "$0")seedms 2> /dev/null
         rm $outfile 2> /dev/null
         printf "#NEXUS\n\nBEGIN PHYLONET;\n" >> $outfile
         printf "SimSeqInNetwork -theta 0.0005 -base (0.2112,0.2888,0.2896,0.2104) -rate (0.2173,0.9798,0.2575,0.1038,1,0.2070) " >> $outfile
         printf " -numgt %d -sitespergt %d " $numgt $sitespergt >> $outfile
         printf ' -mspath "/Users/zhen/Desktop/Zhen/research/phylogenetics/treeAugment/proj/msdir/ms" ' >> $outfile
         printf ' -seqgenpath "/Users/zhen/Desktop/Zhen/research/phylogenetics/treeAugment/proj/Seq-Gen/source/seq-gen" ' >> $outfile
         printf ' -tm <A:A_0;Q:Q_0;L:L_0;Z:Z_0> -truenet "(Z:5.0,(L:2.0,(Q:1.0,A:1.0):1.0):3.0);" -out "' >> $outfile
         printf "$dir/$i/homo/markers.txt">> $outfile
         printf '"  -gtfile "' >> $outfile
         printf "$dir/$i/homo/genetrees.txt" >> $outfile
         printf '";\n' >> $outfile
         printf 'END;\n' >> $outfile
         java -jar /Users/zhen/Desktop/phylonet/build/jar/PhyloNet.jar $outfile
         printf "Generated numgt: %d sitespergt: %d group: %d num: %d\n" $numgt $sitespergt $i $num
         num=$((num+1))
         cp $(dirname "$0")/seedms $dir/$i/

      done
    done
    cp $(dirname "$0")/seedms $dir/$i/
done

