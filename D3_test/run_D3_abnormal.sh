dir_path="/home/zc36/X2/test_outgroup4/net0/"
start=1
end=100

# for((i=$start;i<=$end;++i))
for i in 93 95 96
do  
    for scale in long
    do
        for locus_length in 500
        do
            sub_dir=$dir_path/$scale/$i/heter/$locus_length
            echo $path
            /usr/bin/time -v python D3_test/D3.py $sub_dir/markers.txt 1 2>$sub_dir/D3.err &
            /usr/bin/time -v python D3_test/D3.py $sub_dir/markers.txt 2 2>$sub_dir/D3.err &
        done
    done
done
