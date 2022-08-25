# 0. sanity check
1. polymorphic.py
2. polymorphic_summary_statistics.py

# 1. simulate data using simulate_all_data.sh
1. generate homogeneous gene trees and markers
2. generate heterogeneous gene trees and markers
3. infer gene trees using iqtree and obtain estimated gene tree errors
# 2. Infer species tree
1. prepare_ML_nex.sh
2. run_ML.sh
# 3. summarize obs and exp gene tree probabilities
## 3.1 compute expected gene tree probabilities
Using PhyloNet.jar calgtprob
## 3.2 write {gt obs exp} csv
# 4. compute test statistics chi2, pvalue-chi2, gtest, pvalue-gtest, exactMultinormial
1. write gt probs to a csv: summarize_probs.py
2. compute test statistics- approximate: summarize_5taxa.py
<!-- 3. compute test statistics- exact: exactMultiNom/true_st_true_gt.R  -->
4. plot figures summarize_plot.py
5. triplets multitest 
    5.1 summarize_triplet_probs.py
    5.2 summarize_triplet_stats.py
    5.3 summarize_triplet_multitest.py
    5.4 plot


