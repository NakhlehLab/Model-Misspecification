
import statsmodels.stats.multitest as multitest
import pandas as pd
import numpy as np
import scipy.stats
import sys
from summarize_triplet_stats import GT_DICT, ST_DICT

METHODS = ["bonferroni", "sidak", "holm-sidak", "holm", "simes-hochberg", 
"hommel", "fdr_bh", "fdr_by", "fdr_tsbh", "fdr_tsbky"]

"""
bonferroni : one-step correction
sidak : one-step correction
holm-sidak : step down method using Sidak adjustments
holm : step-down method using Bonferroni adjustments
simes-hochberg : step-up method (independent)
hommel : closed method based on Simes tests (non-negative)
fdr_bh : Benjamini/Hochberg (non-negative)
fdr_by : Benjamini/Yekutieli (negative)
fdr_tsbh : two stage fdr correction (non-negative)
fdr_tsbky : two stage fdr correction (non-negative)
"""

ALPHA = 0.05

def multitest_added(p_values, alpha, method):
    print(type(p_values[0]))
    n_tests = len(p_values)
    if method.lower() == "simes":
        p_values.sort()
        alpha_multi = alpha / n_tests * np.arange(1, n_tests + 1, 1)
        reject = p_values <= alpha_multi
        return any(reject)
    elif method.lower() == "fisher":
        left = -2 * sum([np.log(pi) for pi in p_values])
        right = scipy.stats.chi2.ppf(1-alpha, 2*n_tests)
        if left >= right:
            return True
        else:
            return False
    else:
        raise ValueError('method not recognized')
    

def multi_test_one_replicate(df, value_type, method):
    if method in ["simes", "fisher"]:
        significant = multitest_added(list(df[value_type]), alpha=ALPHA, method=method)
        print("significant", significant)
    else:
        sig_list, pvals_corrected, alphacSidak, alphacBonf = multitest.multipletests(df[value_type], alpha=ALPHA, method=method)
        significant = any(sig_list)
        print(sig_list)
        print(pvals_corrected, alphacSidak, alphacBonf)
        print(significant)
    return significant

def run_all(dir_path, re_end, exp_type, null_reti_num):
    # dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/"
    re_start = 1
    # re_end = 10
    # true_gt = False

    res_data = {"scale":[], "locus_length": [], "replicate":[], "true_gt":[], "pvalue_type":[], "method":[], "significant": []}
    for scale in ["short", "medium", "long"]:
        for locus_length in [500, 2000]:
            for i in range(re_start, re_end+1):
                for true_gt in [True, False]:
                    if exp_type:
                        if true_gt:
                            stat_path = dir_path + scale+"/"+str(i) + "/heter/"+str(locus_length)+"/triplet_stats_sim.csv"
                        else:
                            stat_path = dir_path + scale+"/"+str(i) + "/heter/"+str(locus_length)+"/triplet_stats_iqgt.csv"
                    else:
                        true_st = False
                        stat_path = dir_path + scale+"/"+str(i) + "/heter/"+str(locus_length)+"/triplet_stats_"+ST_DICT[true_st]+"_"+GT_DICT[true_gt]+"_"+null_reti_num+".csv"
                    df = pd.read_csv(stat_path)
                    for value_type in ["pvalue_chisq", "pvalue_gtest", "pvalue_exact"]:
                        for method in ["bonferroni", "simes-hochberg","fdr_bh"]:
                            significant = multi_test_one_replicate(df, value_type, method)
                            res_data["scale"].append(scale)
                            res_data["locus_length"].append(locus_length)
                            res_data["replicate"].append(i)
                            res_data["true_gt"].append(true_gt)
                            res_data["pvalue_type"].append(value_type)
                            res_data["method"].append(method)
                            res_data["significant"].append(significant)
                            
    res_df = pd.DataFrame(res_data)
    res_df.to_csv(dir_path+"/multitest_"+null_reti_num+"_"+str(exp_type)+".csv", index=False)


if __name__ == "__main__":
    # stat_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/long/1/heter/500/triplet_stats_iqgt.csv"
    # value_type = "pvalue_chisq"
    # method = "fdr_bh"
    # df = pd.read_csv(stat_path)
    # multi_test_one_replicate(df, value_type, method)

    dir_path = sys.argv[1]
    re_end = int(sys.argv[2])
    exp_type = (sys.argv[3].lower() == "true")
    null_reti_num = sys.argv[4]
    run_all(dir_path, re_end, exp_type, null_reti_num)