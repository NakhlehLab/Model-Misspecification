# type=AZ
# ./gen_tree_taxa5.sh
for i in {1..10}
do
# rm -r /Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_tall10/036/$i/heter/IQTREE
# rm -r /Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_100/036/$i/homo/IQTREE

# dir=/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_100/001/$i
# ./simulation/sim_replicate.sh $dir/homo $dir/heter "Z" 6 0.0005

dir=/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_tall10/036/$i
./simulation/sim_replicate.sh $dir/homo $dir/heter "Z" 6 0.018

num_gt=1000
dir=/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_gt/$gt/$i/
./simulation/sim_replicate.sh $dir/homo/ $dir/heter/ "Z" 6 0.018 $num_gt False False False

# dir=/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa3_100/036/$i
# ./simulation/sim_replicate.sh $dir/homo $dir/heter "Z" 4 0.018

done