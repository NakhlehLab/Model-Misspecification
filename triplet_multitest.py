import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
import statsmodels.stats.multitest as multitest

gt_name = {True: "truegt", False:"infergt"}
ALPHA = 0.05
SUBSET_SIZE = 3
DEBUG = False
if DEBUG:
    NUM_REPLICATE = 3
else:
    NUM_REPLICATE = 100
NUM_LOCI = 100
GT_TOPO = 3


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

def added_multitest(p_values, alpha, method):
    n_tests = len(p_values)
    if method.lower() == "simes":
        pvals = p_values.sort()
        alpha_multi = alpha / n_tests * np.arange(1, n_tests + 1, 1)
        reject = pvals <= alpha_multi
    elif method.lower() == "fisher":
        left = -2 * sum([np.log(pi) for pi in p_values])
        right = stats.chi2.ppf(1-alpha, 2*n_tests)
        if left >= right:
            return True
    else:
        raise ValueError('method not recognized')
    return reject

def family_wise_error_rate(obs_dir, true_gt, method):
    csv_path = obs_dir+"/triplet_test_"+gt_name[true_gt]+".csv"
    df = pd.read_csv(csv_path)
    st_list = []
    cnt = 0
    for i in range(1, 1 + NUM_REPLICATE):
        df_replicate = df[df["ID"] == i]
        p_values = list(df_replicate["pvalue"])
        sig_list, pvals_corrected, alphacSidak, alphacBonf = multitest.multipletests(p_values, alpha=ALPHA, method=method)
        significant = False
        # print(i)
        # print(list(df_replicate["pvalue"]))
        # print(sig_list)
        for x in sig_list:
            if x == True:
                significant = True
                # print("significant")
                cnt += 1
                break
        # print(significant)
        st_list.append(significant)
    # print(st_list)
    print(cnt*1.0/NUM_REPLICATE)


def test():
    obs_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_tall10/036/"
    csv_path = obs_dir+"/triplet_test_"+gt_name[False]+".csv"
    df = pd.read_csv(csv_path)
    df_replicate = df[df["ID"] == 24]
    p_values = list(df_replicate["pvalue"])
    print(df_replicate)
    sig_list, pvals_corrected, alphacSidak, alphacBonf = multitest.multipletests(p_values, alpha=ALPHA, method="bonferroni")
    print(sig_list)
    for x in sig_list:
        print(x)
        if x == True:
            significant = True
            # print("significant")
            # cnt += 1
            break

if __name__ == "__main__":
    # obs_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_tall10/036/"
    # true_gt = True
    # family_wise_error_rate(obs_dir, true_gt)

    all_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_tall10/"
    for pop_size in ["036", "01", "001"]:
        obs_dir = all_dir + pop_size + "/"
        for true_gt in [True, False]:
            print(pop_size, true_gt)
            family_wise_error_rate(obs_dir, true_gt, method=METHODS[5])

    # test()