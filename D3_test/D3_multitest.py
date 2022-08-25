import statsmodels.stats.multitest as multitest
import pandas as pd
import numpy as np
import scipy.stats
import sys

from summarize_triplet_multitest import multitest_added
METHODS = ["bonferroni", "sidak", "holm-sidak", "holm", "simes-hochberg", 
"hommel", "fdr_bh", "fdr_by", "fdr_tsbh", "fdr_tsbky"]
ALPHA=0.05

def calc_D3_sig(D3_path, D3_multitest_path):
    re_start = 1
    re_end = 100
    df = pd.read_csv(D3_path)
    res_data = {"scale":[], "locus_length": [], "replicate":[], "pvalue_type":[], "method":[], "significant": []}
    
    for method in ["bonferroni", "simes-hochberg","fdr_bh"]:
        for scale in ["short", "medium", "long"]:
            for locus_length in [500, 2000]:
                for i in range(re_start, re_end+1):
                    df_cur = df[(df["scale"] == scale) & (df["locus_length"] == locus_length) & (df["id"] == i)]
                    if method in ["simes", "fisher"]:
                        significant = multitest_added(list(df_cur["D3_pval"]), alpha=ALPHA, method=method)
                    else:
                        sig_list, pvals_corrected, alphacSidak, alphacBonf = multitest.multipletests(df_cur["D3_pval"], alpha=ALPHA, method=method)
                        significant = any(sig_list)
                    res_data["scale"].append(scale)
                    res_data["locus_length"].append(locus_length)
                    res_data["replicate"].append(i)
                    res_data["pvalue_type"].append("D3_pval")
                    res_data["method"].append(method)
                    res_data["significant"].append(significant)
    read_df = pd.DataFrame(res_data)
    read_df.to_csv(D3_multitest_path)
    print(read_df)
    

if __name__ == "__main__":
    scale = "net1_2"
    # D3_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/D3_res.csv"
    # D3_multitest_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/D3_multitest.csv"
    D3_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/"+scale+"/D3_gt_res.csv"
    D3_multitest_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/"+scale+"/D3_gt_multitest.csv"
    
    # D3_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/D3_res_2stage_faster.csv"
    # D3_multitest_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/D3_2stage_multitest.csv"
    calc_D3_sig(D3_path, D3_multitest_path)