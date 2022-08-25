# path=/scratch/zc36/X2/test_outgroup4/net0/
path=/scratch/zc36/X2/mcmc/net1/

re_start=1
re_end=10
subrun=0

for((replicate=$re_start;replicate<=$re_end;++replicate))
do
    for scale in medium long
    do
        for locus_length in 500
        do
            dir_path=$path/$scale/$replicate/heter/$locus_length/
            filename=$dir_path/mcmc
            cd $dir_path
            # /usr/bin/time -v java -jar ~/X2/PhyloNet.jar $filename.nex 1>$filename.out 2>$filename.err & 
            sbatch /scratch/zc36/X2/code/mcmc.slurm $dir_path $subrun false -o $filename.out -e $filename.err

            dir_path=$path/$scale/$replicate/heter/$locus_length/murate/
            mkdir $dir_path
            filename=$dir_path/mcmc
            cd $dir_path
            # /usr/bin/time -v java -jar ~/X2/PhyloNet.jar $filename.nex 1>$filename.out 2>$filename.err &
            sbatch /scratch/zc36/X2/code/mcmc.slurm $dir_path $subrun true -o $filename.out -e $filename.err

            dir_path=$path/$scale/$replicate/homo/$locus_length/
            filename=$dir_path/mcmc
            cd $dir_path
            # /usr/bin/time -v java -jar ~/X2/PhyloNet.jar $filename.nex 1>$filename.out 2>$filename.err &
            sbatch /scratch/zc36/X2/code/mcmc.slurm $dir_path $subrun false -o $filename.out -e $filename.err

            dir_path=$path/$scale/$replicate/homo/$locus_length/murate/
            mkdir $dir_path
            filename=$dir_path/mcmc
            cd $dir_path
            # /usr/bin/time -v java -jar ~/X2/PhyloNet.jar $filename.nex 1>$filename.out 2>$filename.err &
            sbatch /scratch/zc36/X2/code/mcmc.slurm $dir_path $subrun true -o $filename.out -e $filename.err
        done
    done
    wait
done
