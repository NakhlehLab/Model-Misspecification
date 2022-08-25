import pandas as pd
import numpy as np

def test_outgroup_num_mixing():
    df = pd.read_csv("/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_all/net0/gt_mixing_outgroup.csv")
    for scale in ["short", "medium", "long"]:                                      
        for locus_length in [200, 1000]:                                               
            df_cur=df[(df["scale"]==scale) & (df["locus_length"] == locus_length)]       
            print(scale, locus_length, sum(df_cur["cnt_mix"])/len(df_cur["cnt_mix"]))


def test_mean_rf():
    df = pd.read_csv()
    for scale in ["short", "medium", "long"]:
        for locus_length in [200, 1000]:                                               
            df_cur=df[(df["scale"]==scale) & (df["locus_length"] == locus_length)]
            print(df["RF-in_unroot"].mean())

if __name__ == "__main__":
    test_outgroup_num_mixing()