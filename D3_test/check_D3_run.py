import pandas as pd
import sys


def check_D3_runs():
    dir_path = sys.argv[1]
    for scale in ["short", "medium", "long"]:
        for locus_length in [500, 2000]:
            for i in range(1, 101):
                path = dir_path + scale +"/"+ str(i)+"/heter/"+str(locus_length)+"/D3.txt"
                df = pd.read_csv(path)
                if "D3_pval" not in df.columns:
                    print(path)


def check_D3_abnormal():
    path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/D3_res.csv"
    df = pd.read_csv(path)
    df_abnormal = df[(df["locus_length"] == 2000) & (df["scale"] == "long") &  (df["D3_pval"] < 0.05)]
    print(df)
    print(df_abnormal)
    abs_z_score = abs(df_abnormal["D3_mean"] / df_abnormal["D3_stdev"])
    print(abs_z_score)

check_D3_abnormal()