dir_path="/home/zc36/X2/test_outgroup4/net0/"
start=1
end=100
remedy=3
for((i=$start;i<=$end;++i))
do  
    for scale in short medium long
    do
        for locus_length in 500 2000
        do
            sub_dir=$dir_path/$scale/$i/heter/$locus_length
            echo $path
            /usr/bin/time -v python D3_test/D3.py $sub_dir/markers.txt $locus_length $remedy 2>$sub_dir/D3.err &
        done
    done
    wait
done
