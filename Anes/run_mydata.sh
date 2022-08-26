dir_path="/home/zc36/X2/test_outgroup4/net1_2"
start=1
end=100
ingroup="ingroup"
export PATH="$PATH:/home/zc36/lib/julia-1.7.2/bin"

for((i=$start;i<=$end;i=$(($i+5))))
do  
    for((x=$i;x<$(($i+5));++x))
    do
        for scale in short medium long
        do
            for locus_length in 500 2000
            do
                sub_dir=$dir_path/$scale/$x/heter/$locus_length
                echo $sub_dir
        		if [[ -n $ingroup ]];
        		then
        			/usr/bin/time -v julia --project=/home/zc36/X2/tools/QuartetNetworkGoodnessFit.jl mydata.jl $sub_dir/quartet_probs_true_ingroup.csv $sub_dir/ML_true_0.out $sub_dir/quarnetGoF_true_0_stat_ingroup.csv $sub_dir/quarnetGoF_true_0_res_ingroup.txt true 2>$sub_dir/quarnetGoF_true_0_ingroup.err 1>$sub_dir/quarnetGoF_true_0_ingroup.out  &
                    /usr/bin/time -v julia --project=/home/zc36/X2/tools/QuartetNetworkGoodnessFit.jl mydata.jl $sub_dir/quartet_probs_iq_ingroup.csv $sub_dir/ML_iq_0.out $sub_dir/quarnetGoF_iq_0_stat_ingroup.csv $sub_dir/quarnetGoF_iq_0_res_ingroup.txt true 2>$sub_dir/quarnetGoF_iq_0_ingroup.err 1>$sub_dir/quarnetGoF_iq_0_ingroup.out  &
        		else
                    /usr/bin/time -v julia --project=/home/zc36/X2/tools/QuartetNetworkGoodnessFit.jl mydata.jl $sub_dir/quartet_probs_true.csv $sub_dir/ML_true_0.out $sub_dir/quarnetGoF_true_0_stat.csv $sub_dir/quarnetGoF_true_0_res.txt false 2>$sub_dir/quarnetGoF_true_0.err 1>$sub_dir/quarnetGoF_true_0.out  &
                    /usr/bin/time -v julia --project=/home/zc36/X2/tools/QuartetNetworkGoodnessFit.jl mydata.jl $sub_dir/quartet_probs_iq.csv $sub_dir/ML_iq_0.out $sub_dir/quarnetGoF_iq_0_stat.csv $sub_dir/quarnetGoF_iq_0_res.txt false 2>$sub_dir/quarnetGoF_iq_0.err 1>$sub_dir/quarnetGoF_iq_0.out  &
        		fi
            done

        done
    done
    wait
done