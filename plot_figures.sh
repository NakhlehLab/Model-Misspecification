
# path="/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_tall10/036/"
path="/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_gt/100/"
num_taxa=5

# path="/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa3_tall/001/"
# num_taxa=3

# python polymorphic.py $path $num_taxa 018
# python gene_tree_error.py $path 018
# python process_gt_error.py $path 018

if [ $num_taxa -eq 5 ];
then
    # java -jar PhyloNet.jar $path $num_taxa
    # python branch_length.py $path True
    # python branch_length.py $path False
    python eval_5taxa.py $path 018
else
    python eval_3taxa.py $path

fi

# 