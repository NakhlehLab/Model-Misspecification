
path=/home/zc36/X2/taxa5
for replicate in {1..10}
do
filename=$path/$replicate/heter/ML_iq_0
/usr/bin/time -v java -jar ~/X2/PhyloNet.jar  $filename.nex 2>$filename.err 1>$filename.out & 
filename=$path/$replicate/heter/ML_true_0
/usr/bin/time -v java -jar ~/X2/PhyloNet.jar  $filename.nex 2>$filename.err 1>$filename.out & 
done