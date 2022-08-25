#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
num=0
dir=$1
outfile=$dir/gen_net.nex
phylonet=/home/zc36/X2/PhyloNet.jar
numgt=$6
sitespergt=$2
poprate=$3
start=$4
end=$5
    for((i=$start;i<=$end;++i))
      do
         mkdir $dir/$i 2> /dev/null
         mkdir $dir/$i/homo 2> /dev/null
	       mkdir $dir/$i/homo/$sitespergt 2> /dev/null
         rm  $(dirname "$0")seedms 2> /dev/null
         rm $outfile 2> /dev/null
         printf "#NEXUS\n\nBEGIN PHYLONET;\n" >> $outfile
         printf "SimSeqInNetwork -theta %.4f -base (0.2112,0.2888,0.2896,0.2104) -rate (0.2173,0.9798,0.2575,0.1038,1,0.2070) " $poprate >> $outfile
         printf " -numgt %d -sitespergt %d " $numgt $sitespergt >> $outfile
         printf ' -mspath "/home/zc36/lib/msdir/ms" ' >> $outfile
         printf ' -seqgenpath "/home/zc36/lib/Seq-Gen-1.3.4/source/seq-gen" ' >> $outfile
         printf ' -tm <C:C_0;G:G_0;L:L_0;Q:Q_0;R:R_0;Z:Z_0> -truenet "(((((Q:0.15,#H1:0.05::0.3):0.15,R:0.3):0.5,(L:0.1)#H1:0.7::0.7):0.2,(C:0.4,G:0.4):0.6):10,Z:11);" -out "' >> $outfile
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

    cp $(dirname "$0")/seedms $dir/$i/


