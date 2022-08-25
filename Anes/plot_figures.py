import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from pathlib import Path
import numpy as np
import scipy.stats

def plot_outlier_pval(path_iqgt, path_truegt, ignore_outgroup=False, num_reti=0):
    figpath = str(Path(path_iqgt).parent.absolute())+ "/fig_quarnet_outlier_pval_"+str(num_reti)+".pdf"

    df_iqgt = pd.read_csv(path_iqgt)
    df_truegt = pd.read_csv(path_truegt)

    df_iqgt["true_gt"] = [False for x in range(len(df_iqgt.index))]
    df_truegt["true_gt"] = [True for x in range(len(df_truegt.index))]
    df = pd.concat([df_iqgt, df_truegt])

    gt_error_name = {0:"true gt", 1: "inferred gt from 2000 bps", 2: "inferred gt from 500 bps"}
    gt_error_config = {0:{"locus_length": 500, "true_gt": True}, 1:{"locus_length": 2000, "true_gt": False}, 2:{"locus_length": 500, "true_gt": False} }

    if ignore_outgroup:
        # df = df[(df["t1"] != "Z") & (df["t2"] != "Z") & (df["t3"] != "Z") &(df["t4"] != "Z")]
        figpath = str(Path(path_iqgt).parent.absolute())+ "/fig_quarnet_outlier_pval_"+str(num_reti)+"_ingroup.pdf"
    else:
        figpath = str(Path(path_iqgt).parent.absolute())+ "/fig_quarnet_outlier_pval_"+str(num_reti)+".pdf"

    fig, axes = plt.subplots(3, 3, figsize=(18, 10))
    scale_list = ["short", "medium", "long"]
    for scale_idx in [0, 1, 2]:
        scale = scale_list[scale_idx]
        # for gt_error in [0, 1, 2, 3]:
        for gt_error in [0, 1, 2]:
            df_cur = df[(df["scale"] == scale)
            & (df["true_gt"] == gt_error_config[gt_error]["true_gt"])
            & (df["locus_length"] == gt_error_config[gt_error]["locus_length"])]
            
            
            for i, row in df_cur.iterrows():
                df_cur.at[i,"p_value"] = "%.5f" % row["p_value"]
            print(scale, gt_error)
            # print(df_cur["p_value"])
            axes[scale_idx][gt_error].set_xlim(-0.05,1.05)
            sns.histplot(x="p_value",
                    binwidth=0.05,
                    binrange=(0,1),
                    data=df_cur, 
                    ax=axes[scale_idx][gt_error], 
                    color='#607c8e',
                    edgecolor="white",
                    stat="probability"
                    )
            significance_count = df_cur.loc[df_cur["p_value"] < 0.05,"p_value"].count()/df_cur["p_value"].count()
            print(f"{significance_count} -- Proportion of significance:{significance_count}")
            axes[scale_idx][gt_error].set_ylim(0, 1.0)

            axes[scale_idx][gt_error].set_xlabel("", fontsize=16)
            axes[scale_idx][gt_error].set_ylabel("", fontsize=16)
            axes[scale_idx][gt_error].tick_params(axis='y', labelsize=12)
            axes[scale_idx][gt_error].tick_params(axis="x", labelsize=12)
            axes[0][gt_error].set_title(gt_error_name[gt_error], size=20, fontweight="bold")
            
        axes[scale_idx][0].set_ylabel(scale_list[scale_idx], size=16)
    
    plt.tight_layout()
    plt.savefig(figpath, dpi=300)


def plot_summary_pval(path_iqgt, path_truegt, ignore_outgroup = False, num_reti=0):
   
    df_iqgt = pd.read_csv(path_iqgt)
    df_truegt = pd.read_csv(path_truegt)

    df_iqgt["true_gt"] = [False for x in range(len(df_iqgt.index))]
    df_truegt["true_gt"] = [True for x in range(len(df_truegt.index))]
    df = pd.concat([df_iqgt, df_truegt])

    gt_error_name = {0:"true gt", 1: "inferred gt from 2000 bps", 2: "inferred gt from 500 bps"}
    gt_error_config = {0:{"locus_length": 500, "true_gt": True}, 1:{"locus_length": 2000, "true_gt": False}, 2:{"locus_length": 500, "true_gt": False} }

    fig, axes = plt.subplots(3, 3, figsize=(18, 10))
    if ignore_outgroup:
        # df = df[(df["t1"] != "Z") & (df["t2"] != "Z") & (df["t3"] != "Z") &(df["t4"] != "Z")]
        figpath = str(Path(path_iqgt).parent.absolute())+ "/fig_quarnet_summary_pval_"+str(num_reti)+"_ingroup.pdf"
    else:
        figpath = str(Path(path_iqgt).parent.absolute())+ "/fig_quarnet_summary_pval_"+str(num_reti)+".pdf"

    scale_list = ["short", "medium", "long"]
    for scale_idx in [0, 1, 2]:
        scale = scale_list[scale_idx]
        # for gt_error in [0, 1, 2, 3]:
        for gt_error in [0, 1, 2]:
            df_cur = df[(df["scale"] == scale)
            & (df["true_gt"] == gt_error_config[gt_error]["true_gt"])
            & (df["locus_length"] == gt_error_config[gt_error]["locus_length"])]
            
            
            for i, row in df_cur.iterrows():
                df_cur.at[i,"p_value"] = "%.5f" % row["p_value"]
            print(scale, gt_error)
            # print(df_cur["p_value"])
            significance_count = df_cur.loc[df_cur["p_value"] < 0.05,"p_value"].count()/df_cur["p_value"].count()
            print(f"{significance_count} -- Proportion of significance:{significance_count}")
            axes[scale_idx][gt_error].set_xlim(-0.05,1.05)
            sns.histplot(x="p_value",
                    binwidth=0.05,
                    binrange=(0,1),
                    data=df_cur, 
                    ax=axes[scale_idx][gt_error], 
                    color='#607c8e',
                    edgecolor="white",
                    stat="probability"
                    )
            axes[scale_idx][gt_error].set_ylim(0, 1.0)

            axes[scale_idx][gt_error].set_xlabel("", fontsize=16)
            axes[scale_idx][gt_error].set_ylabel("", fontsize=16)
            axes[scale_idx][gt_error].tick_params(axis='y', labelsize=12)
            axes[scale_idx][gt_error].tick_params(axis="x", labelsize=12)
            axes[0][gt_error].set_title(gt_error_name[gt_error], size=20, fontweight="bold")
            
        axes[scale_idx][0].set_ylabel(scale_list[scale_idx], size=16)
    
    plt.tight_layout()
    plt.savefig(figpath, dpi=300)


if __name__ == "__main__":
    ignore_outgroup = True
    num_reti=1
    # net="net0"
    if num_reti == 0:
        if ignore_outgroup:
            outlier_path_iq = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/quarnetGoF_iq_outlierp_ingroup.csv"
            outlier_path_true = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/quarnetGoF_true_outlierp_ingroup.csv"
            summary_path_iq = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/quarnetGoF_iq_summary_ingroup.csv"
            summary_path_true = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/quarnetGoF_true_summary_ingroup.csv"
        else:
            outlier_path_iq = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/new_res_comp/quarnetGoF_iq_outlierp.csv"
            outlier_path_true = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/new_res_comp/quarnetGoF_true_outlierp.csv"
            summary_path_iq = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/new_res_comp/quarnetGoF_iq_summary.csv"
            summary_path_true = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/new_res_comp/quarnetGoF_true_summary.csv"
    elif num_reti == 1:
        if ignore_outgroup:
            outlier_path_iq = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net1_2/quarnetGoF_iq_outlierp_ingroup.csv"
            outlier_path_true = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net1_2/quarnetGoF_true_outlierp_ingroup.csv"
            summary_path_iq = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net1_2/quarnetGoF_iq_summary_ingroup.csv"
            summary_path_true = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net1_2/quarnetGoF_true_summary_ingroup.csv"
        else:
            outlier_path_iq = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net1_2/quarnetGoF_iq_outlierp.csv"
            outlier_path_true = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net1_2/quarnetGoF_true_outlierp.csv"
            summary_path_iq = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net1_2/quarnetGoF_iq_summary.csv"
            summary_path_true = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net1_2/quarnetGoF_true_summary.csv"
    plot_outlier_pval(outlier_path_iq, outlier_path_true, ignore_outgroup)
    plot_summary_pval(summary_path_iq, summary_path_true, ignore_outgroup)