import pandas as pd
import seaborn as sns
from scipy.special import rel_entr
import sys
from scipy.stats import power_divergence
from summarize_triplet_probs import divide_triplet
from calculator import *

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter


DEBUG=False
if DEBUG:
    NUM_REPLICATE = 3
else:
    NUM_REPLICATE = 100
NUM_GT = 10000

N_TAXA = 3
DDOF = 1

GT_DICT = {True: "truegt", False: "iqgt"}
ST_DICT = {True: "truest", False: "mlst"}


TRIPLETS = divide_triplet()

r = ro.r
r['source']('exactMultiNom/exact_mc.R')
exact_function_r = ro.globalenv['exact_multi_nom_1rep']


def test_statistics_one_replicate_triplets(prob_path, stat_path):
    df = pd.read_csv(prob_path)
    triplet_pvalue = {"triplet":[], "chisq": [], "pvalue_chisq":[], "gtest":[], "pvalue_gtest":[], "pvalue_exact":[]}
    for triplet in TRIPLETS:
        print(triplet)
        
        df_tri = df[df["triplet"]==str(triplet)]
        print(df_tri)
        obs = [NUM_GT*i for i in df_tri["prob_obs"]]
        exp = [NUM_GT*i for i in df_tri["prob_exp"]]
        chisq, pvalue_chisq = power_divergence(obs, exp, DDOF, lambda_='pearson')
        gtest, pvalue_gtest = power_divergence(obs, exp, DDOF, lambda_='log-likelihood')
        triplet_pvalue["triplet"].append(triplet)
        triplet_pvalue["chisq"].append(chisq)
        triplet_pvalue["pvalue_chisq"].append(pvalue_chisq)
        triplet_pvalue["gtest"].append(gtest)
        triplet_pvalue["pvalue_gtest"].append(pvalue_gtest)
        with localconverter(ro.default_converter + pandas2ri.converter):
            df_tri_r = ro.conversion.py2rpy(df_tri)
        exact_pvalue = exact_function_r(df_tri_r, NUM_GT, "exact")
        print(exact_pvalue)
        triplet_pvalue["pvalue_exact"].append(exact_pvalue[0])
    
    df_tri_stat = pd.DataFrame(triplet_pvalue)
    df_tri_stat.to_csv(stat_path)
    return df_tri_stat

def run_all(dir_path, re_end, exp_type, null_reti_num):
    # dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/"
    re_start = 1
    # re_end = 10
    # true_gt = False
    
    df_tri_list = []
    for scale in ["short", "medium", "long"]:
        for locus_length in [500, 2000]:
            for i in range(re_start, re_end+1):
                for true_gt in [True, False]:
                    if exp_type:
                        true_st = False
                        if true_gt:
                            triplet_path = dir_path + scale+"/"+str(i) + "/heter/"+str(locus_length)+"/triplet_probs_sim.csv"
                            stat_path = dir_path + scale+"/"+str(i) + "/heter/"+str(locus_length)+"/triplet_stats_sim.csv"
                        else:
                            triplet_path = dir_path + scale+"/"+str(i) + "/heter/"+str(locus_length)+"/triplet_probs_iqgt.csv"
                            stat_path = dir_path + scale+"/"+str(i) + "/heter/"+str(locus_length)+"/triplet_stats_iqgt.csv"
                    else:
                        true_st = False
                        triplet_path =  dir_path+scale+"/"+str(i)+"/heter/"+str(locus_length)+"/triplet_probs_"+ST_DICT[true_st]+"_"+GT_DICT[true_gt]+"_"+null_reti_num+".csv"
                        stat_path = dir_path + scale+"/"+str(i) + "/heter/"+str(locus_length)+"/triplet_stats_"+ST_DICT[true_st]+"_"+GT_DICT[true_gt]+"_"+null_reti_num+".csv"
                    df_tri_stat = test_statistics_one_replicate_triplets(triplet_path, stat_path)
                    
                    df_tri_stat["scale"] = [scale for x in range(len(df_tri_stat["triplet"]))]
                    df_tri_stat["locus_length"] = [locus_length for x in range(len(df_tri_stat["triplet"]))]
                    df_tri_stat["replicate"] = [i for x in range(len(df_tri_stat["triplet"]))]
                    df_tri_stat["true_gt"] = [true_gt for x in range(len(df_tri_stat["triplet"]))]
                    df_tri_stat["true_st"] = [true_st for x in range(len(df_tri_stat["triplet"]))]
                    df_tri_list.append(df_tri_stat)
    res_df = pd.concat(df_tri_list, ignore_index=True)
    res_df.to_csv(dir_path+"/triplet_stats_all_exp"+str(exp_type)+"_"+null_reti_num+".csv")

if __name__ == "__main__":
    # prob_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/long/1/heter/500/triplet_probs_iqgt.csv"
    # test_statistics_one_replicate_triplets(prob_path)
    dir_path = sys.argv[1]
    re_end = int(sys.argv[2])
    exp_type = (sys.argv[3].lower() == "true")
    null_reti_num = sys.argv[4]
    run_all(dir_path, re_end, exp_type, null_reti_num)

    # triplet_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/long/1/heter/1000/triplet_probs_truest_iqgt_0.csv"
    # stat_path = "./test_stat.csv"
    # test_statistics_one_replicate_triplets(triplet_path, stat_path) 