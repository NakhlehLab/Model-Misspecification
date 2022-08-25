#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
num=0
dir=/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/test/
outfile=$dir/gen_net.nex
phylonet=/Users/zhen/Desktop/phylonet/build/jar/PhyloNet.jar
for numgt in 10000
do
    for sitespergt in 500
    do
    for i in {1..2}
#    for i in 1 2
      do
         mkdir $dir/$i 2> /dev/null
         mkdir $dir/$i/homo 2> /dev/null
         mkdir $dir/$i/homo/$sitespergt 2> /dev/null
         rm  $(dirname "$0")seedms 2> /dev/null
         rm $outfile 2> /dev/null
         printf "#NEXUS\n\nBEGIN PHYLONET;\n" >> $outfile
         printf "SimSeqInNetwork -theta 0.005 -base (0.2112,0.2888,0.2896,0.2104) -rate (0.2173,0.9798,0.2575,0.1038,1,0.2070) " >> $outfile
         printf " -numgt %d -sitespergt %d " $numgt $sitespergt >> $outfile
         printf ' -mspath "/Users/zhen/Desktop/Zhen/research/phylogenetics/treeAugment/proj/msdir/ms" ' >> $outfile
         printf ' -seqgenpath "/Users/zhen/Desktop/Zhen/research/phylogenetics/treeAugment/proj/Seq-Gen/source/seq-gen" ' >> $outfile
         printf ' -truenet "(((((Q:0.2)#H1:0.1::0.7,R:0.3)I3:0.5,(L:0.4,#H1:0.2::0.3):0.4)I1:0.2,(G:0.4,C:0.4)I2:0.6)I0:4,Z:5);" -out "' >> $outfile
        #  printf ' -truenet "(((B:1.488875,E:1.488875):0.300988,(A:1.489880,(C:0.881394,D:0.881394):0.608486):0.299983):1.7899,Z:3.579763);" -out "' >> $outfile
         printf "$dir/$i/homo/$sitespergt/markers.txt">> $outfile
         printf '"  -gtfile "' >> $outfile
         printf "$dir/$i/homo/$sitespergt/genetrees.txt" >> $outfile
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


