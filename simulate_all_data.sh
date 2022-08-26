# simulate homogeneous gene trees
dir=/home/zc36/X2/taxa5/net0/
re_start=1
re_end=100
num_loci=10000
for level in level0 level1_2
do
	for pair in "short 0.025" "medium 0.0125" "long 0.005"
	do
		a=( $pair )
		scale=${a[0]}
		poprate=${a[1]}
		for locus_length in 500 2000
		do
			mkdir $dir/$scale/
			simulation_script/$level/gen_tree_taxa5_$scale.sh $dir/$scale/ $locus_length $poprate $re_start $re_end $num_loci
			#simulate markers with locus heterogeneous
			for((i=$re_start;i<=$re_end;++i))
			do
				rep_dir=$dir/$scale/$i
				echo rep $i
				simulation/sim_gt.sh $rep_dir/homo/ $rep_dir/heter/ 6 10000 $poprate "Z" $locus_length
			done
			#infer iqtree
			for((start=$re_start;start<=$re_end;start=$(($start+10))))
			do
				for((i=$start;i<$(($start+10));++i))
				do 
						echo iqtree $i
						rep_dir=$dir/$scale/$i
						./simulation/infer_iqtree.sh $rep_dir/homo/ $rep_dir/heter/ 6 10000 $poprate "Z" $locus_length &
					done
				wait
			done
		done
	done
done

/usr/bin/time -v python  collapse_short_gt_branch.py $dir_path $num_replicate short 2>$dir_path/medium/collapse_gt_polytomy.err
/usr/bin/time -v python  collapse_short_gt_branch.py $dir_path $num_replicate medium 2>$dir_path/medium/collapse_gt_polytomy.err
/usr/bin/time -v python  collapse_short_gt_branch.py $dir_path $num_replicate long   2>$dir_path/long/collapse_gt_polytomy.err

