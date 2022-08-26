dir_path="/home/zc36/X2/test_outgroup4/net1_2"
start=1
end=100

for((rep=$start;rep<=$end;rep=$(($rep+40))))
do

for((i=$rep;i<=$(($rep+40)) && i<=$end;++i))
do
        for pair in "short 0.025" "medium 0.0125" "long 0.005"
do
        a=( $pair )
        scale=${a[0]}
        poprate=${a[1]}
        for locus_length in 500 2000
        do
            sub_dir=$dir_path/$scale/$i/homo/$locus_length
            echo $path
            /usr/bin/time -v python D3_test/D3_gt.py $sub_dir/genetrees.txt $locus_length $poprate 2>$sub_dir/D3_gt.err &
        done
    done
done
    wait
done
~    