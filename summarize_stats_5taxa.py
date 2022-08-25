import pandas as pd
import seaborn as sns
from dendropy.calculate import treecompare
from scipy.special import rel_entr
import scipy.stats as stats
import sys
from scipy.stats import power_divergence
from pathlib import Path
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from pathlib import Path

from calculator import *
from eval_5taxa import DDOF

N_TAXA = 5
DEBUG=False
if DEBUG:
    NUM_REPLICATE = 3
else:
    NUM_REPLICATE = 100
NUM_GT = 10000
NUM_GT_TOPO = 105

DDOF_dict = {0: 3, 1: 4, 2: 9}
# DDOF = N_TAXA - 2
DOF = NUM_GT_TOPO - 1 


GT_DICT = {True: "truegt", False: "iqgt"}
ST_DICT = {True: "truest", False: "mlst"}

r = ro.r
r['source']('exactMultiNom/exact_mc.R')
exact_function_r = ro.globalenv['exact_multi_nom_1rep']


def test_statistics_one_replicate(df, DDOF):

    obs = [NUM_GT*i for i in df["prob_obs"]]
    exp = [NUM_GT*i for i in df["prob_exp"]]
    
    chisq, pvalue_chisq = power_divergence(obs, exp, DDOF, lambda_='pearson')
    gtest, pvalue_gtest = power_divergence(obs, exp, DDOF, lambda_='log-likelihood')
    print(chisq, pvalue_chisq, gtest, pvalue_gtest)
    return chisq, pvalue_chisq, gtest, pvalue_gtest

def summarize_test_statistics(dir_path, re_end, null_reti_num):
    # dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/"
    re_start = 1
    DDOF = DDOF_dict[int(null_reti_num)]
    # re_end = 10
    res_data = {"scale":[], "locus_length":[], "replicate":[], "chisq":[], "true_st":[], "true_gt":[], "pvalue_chisq":[], "gtest":[], "pvalue_gtest":[], "pvalue_exact_mc":[]}
    for scale in ["short", "medium", "long"]:
        for locus_length in [500, 2000]:
            for i in range(re_start, re_end+1):
                for true_st in [False]:
                    for true_gt in [True, False]:
                        if true_st and true_gt:
                            prob_path = dir_path+scale+"/"+str(i)+"/heter/"+str(locus_length)+"/gt_prob_all_truest_truegt_"+null_reti_num+".csv" 
                        elif true_st and not true_gt:
                            prob_path = dir_path+scale+"/"+str(i)+"/heter/"+str(locus_length)+"/gt_prob_all_truest_iqgt_"+null_reti_num+".csv" 
                        elif true_gt:
                            prob_path = dir_path+scale+"/"+str(i)+"/heter/"+str(locus_length)+"/gt_prob_all_mlst_truegt_"+null_reti_num+".csv" 
                        else:
                            prob_path = dir_path+scale+"/"+str(i)+"/heter/"+str(locus_length)+"/gt_prob_all_mlst_iqgt_"+null_reti_num+".csv" 
                        df = pd.read_csv(prob_path)
                        chisq, pvalue_chisq, gtest, pvalue_gtest = test_statistics_one_replicate(df, DDOF)
                        res_data["scale"].append(scale)
                        res_data["locus_length"].append(locus_length)
                        res_data["replicate"].append(i)
                        res_data["true_gt"].append(true_gt)
                        res_data["true_st"].append(true_st)
                        res_data["chisq"].append(chisq)
                        res_data["pvalue_chisq"].append(pvalue_chisq)
                        res_data["gtest"].append(gtest)
                        res_data["pvalue_gtest"].append(pvalue_gtest)
                        with localconverter(ro.default_converter + pandas2ri.converter):
                            df_r = ro.conversion.py2rpy(df)
                        exact_pvalue = exact_function_r(df_r, NUM_GT, "Monte-Carlo")
                        print(exact_pvalue)
                        res_data["pvalue_exact_mc"].append(exact_pvalue[0])
    df = pd.DataFrame(res_data)
    df.to_csv(dir_path + "/statistics_"+null_reti_num+".csv")
   
                
if __name__ == "__main__":
    
    # obs_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/short/10/heter/500/gt_prob_iqtree.csv"
    # exp_path="/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/short/10/heter/500/gt_prob_exp_iq_0.csv"
    # process_data(obs_path, exp_path)
    # prob_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/short/1/heter/500/gt_prob_all_true_0.csv"
    # test_statistics_one_replicate(prob_path)
    dir_path = sys.argv[1]
    re_end = int(sys.argv[2])
    null_reti_num = sys.argv[3]

    summarize_test_statistics(dir_path, re_end, null_reti_num)
