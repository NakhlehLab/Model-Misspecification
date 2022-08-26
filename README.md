# 1. simulate data using simulate_all_data.sh
Use simulate_all_data.sh, need to specify your own path to ms and seqgen in simulation/model.py

1. generate homogeneous gene trees and markers
2. generate heterogeneous gene trees and markers
3. infer gene trees using iqtree and obtain estimated gene tree errors
4. randomly resolve short branches 

# 2. Infer species tree for full network using PhyloNet
Use command Infernetwork_ML in Phylonet to infer network
1. prepare_ML_nex.sh
2. run_ML.sh

# 3. Analysis 
Use summarize_result.sh

## 3.1 compute expected gene tree probabilities
You can use CalGTProb Command in PhyloNet for each gene tree topology

## 3.2 compute test statistics on full network

1. write gt probs to a csv: summarize_probs.py
2. compute test statistics- approximate: summarize_stats_5taxa.py
4. plot figures summarize_plot.py

## 3.3 triplets multitest 
    5.1 summarize_triplet_probs.py
    5.2 summarize_triplet_stats.py
    5.3 summarize_triplet_multitest.py
    5.4 plot use summarize_triplet_plot.py


