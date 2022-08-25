# simulate homogeneous gene trees
dir=/home/zc36/X2/taxa5/net1/
re_start=1
re_end=100
for pair in "short 0.025" "medium 0.0125" "long 0.005"
#for pair in "short 0.005" "medium 0.005"
do
	a=( $pair )
	scale=${a[0]}
	poprate=${a[1]}
	for locus_length in 500 2000
	do
		mkdir $dir/$scale/
#		./gen_tree_taxa5_$scale.sh $dir/$scale/ $locus_length $poprate $re_start $re_end
		
		#simulate markers with locus heterogeneous
		for((i=$re_start;i<=$re_end;++i))
		do
			rep_dir=$dir/$scale/$i
			echo rep $i
#			simulation/sim_gt.sh $rep_dir/homo/ $rep_dir/heter/ 6 10000 $poprate "Z" $locus_length
#			rm -r $rep_dir/heter/$locus_length/IQTREE
		done
		#infer iqtree
		for((start=$re_start;start<=$re_end;start=$(($start+10))))
                do
                  	for((i=$start;i<$(($start+10));++i))
                        do 
                           	echo iqtree $i
				rep_dir=$dir/$scale/$i
#				rm -r $rep_dir/heter/$locus_length/IQTREE
                               ./simulation/infer_iqtree.sh $rep_dir/homo/ $rep_dir/heter/ 6 10000 $poprate "Z" $locus_length &
                        done
			wait
                done
	done
done
	
