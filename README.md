# Statistical test
## 1. simulate data using simulate_all_data.sh
Use simulation/simulate_all_data.sh.  Need to specify your own path to ms and seqgen in simulation/model.py
The simulation process is described below:
1. generate homogeneous gene trees and markers
2. generate heterogeneous gene trees and markers
3. infer gene trees using iqtree and obtain estimated gene tree errors
4. randomly resolve short branches 

## 2. Infer species tree for full network PhyloNet
Use command Infernetwork_ML in Phylonet to infer network
prepare_ML_nex.sh, run_ML.sh are two scripts may help

## 3. Analysis 
Use summarize_result.sh

### 3.1 compute expected gene tree probabilities
You can use CalGTProb Command in PhyloNet for each gene tree topology

### 3.2 compute test statistics on full network

1. write gt probs to a csv: summarize_probs.py
2. compute test statistics- approximate: summarize_stats_5taxa.py
3. plot figures with summarize_plot.py

### 3.3 triplets multitest 
1. summarize_triplet_probs.py
2. summarize_triplet_stats.py
3. summarize_triplet_multitest.py
4. plot figures with summarize_triplet_plot.py

## 4. Comparison
### D3
The scripts to run D3 are: D3_test/run_D3.sh and D3_test/run_D3_gt.sh
### QuartetGoodnessFit
1. Use run_quartet_probs.sh to prepare input files.
2. Put mydata.jl and run_mydata.sh in to directory of QuartetGoodnessFit and run 
3. summarize results with summarize_pval.py

# Supporting per-locus rate heterogenety in MCMC_SEQ
This is implemented in the code base of PhyloNet.

Download:
https://github.com/NakhlehLab/PhyloNet/releases

Nexus line:
````
MCMC_SEQ -cl 20000000 -bl 10000000 -sf 5000 -sd 12345 -pl 24 -tm <C:C_0;L:L_0;R:R_0;Q:Q_0>; 
MCMC_SEQ -cl 20000000 -bl 10000000 -sf 5000 -sd 12345 -pl 24 -tm <C:C_0;L:L_0;R:R_0;Q:Q_0> -murate; 
````

The -murate switch is to enable the sampling of varying substitution rates.
