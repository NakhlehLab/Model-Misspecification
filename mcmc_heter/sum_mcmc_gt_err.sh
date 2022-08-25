pypath_mcmc=~/X2/code/mcmc_heter/MCMC_result_summarize.py

mcmc_dir=/home/zc36/X2/mcmc/
locus_length=2000
for net in net0 net1_2
do
    for h in heter homo
    do
        for murate in "murate" ""
        do
            for i in {1..10}
            do  
                path=$mcmc_dir/$net/long/$i/$h/$locus_length/$murate
                true_tree_path=$mcmc_dir/$net/long/$i/homo/$locus_length/genetrees.txt
                echo $path $true_tree_path
                # python $pypath_mcmc $path $true_tree_path 2000000 5000 50

            done
        done
    done
done