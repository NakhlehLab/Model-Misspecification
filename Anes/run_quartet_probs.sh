for i in {1..100}
do
    for scale in short medium long
    do
        for locus_length in 500 2000
        do
            
            
            /usr/bin/time -v python Anes/summarize_quartet_probs.py $1 $scale $locus_length $i  &
             

        done
    done
    wait
done