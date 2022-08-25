pypath_mcmc=~/X2/code/gt_error_profiling.py

mcmc_dir=/home/zc36/X2/mcmc/
locus_length=2000
for net in net0 net1_2
do
    for h in heter homo
    do
        for i in {1..10}
        do  
            iq_tree_path=$mcmc_dir/$net/long/$i/$h/$locus_length/rooted_iqtree_resolved.txt
            true_tree_path=$mcmc_dir/$net/long/$i/homo/$locus_length/genetrees_pruned.txt
            echo $iq_tree_path $true_tree_path
            python $pypath_mcmc $true_tree_path $iq_tree_path Z $mcmc_dir/$net/long/$i/$h/$locus_length/iqtree_err.csv 0.005 2>iqerr.err &
        done
    done
    wait
done