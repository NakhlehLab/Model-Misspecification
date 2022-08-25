import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from pathlib import Path
import numpy as np
import scipy.stats
from summarize_stats_5taxa import DDOF, DOF

def plot(path, stat_type, num_reti):
    figpath = str(Path(path).parent.absolute())+ "/fig_"+stat_type+"_"+str(num_reti)+".pdf"
    df = pd.read_csv(path)
    # fig, axes = plt.subplots(3, 4, figsize=(20, 10))
    fig, axes = plt.subplots(3, 3, figsize=(18, 10))

    
    # gt_error_name = {0:"true st true gt", 1:"true gt", 2: "inferred gt from 1000 bps", 3: "inferred gt from 500 bps"}
    # gt_error_config = {0:{"true_st": True, "true_gt": True, "locus_length": 500}, 1:{"true_st": False, "true_gt": True, "locus_length": 500}, 2:{"true_st": False, "true_gt": False, "locus_length": 1000}, 3:{"true_st": False, "true_gt": False, "locus_length": 500} }
    gt_error_name = {0:"true gt", 1: "inferred gt from 2000 bps", 2: "inferred gt from 500 bps"}
    gt_error_config = {0:{"true_st": False, "true_gt": True, "locus_length": 500}, 1:{"true_st": False, "true_gt": False, "locus_length": 2000}, 2:{"true_st": False, "true_gt": False, "locus_length": 500} }

    scale_list = ["short", "medium", "long"]
    for scale_idx in [0, 1, 2]:
        scale = scale_list[scale_idx]
        # for gt_error in [0, 1, 2, 3]:
        for gt_error in [0, 1, 2]:
            df_cur = df[(df["scale"] == scale) & (df["true_st"] == gt_error_config[gt_error]["true_st"]) 
            & (df["true_gt"] == gt_error_config[gt_error]["true_gt"])
            & (df["locus_length"] == gt_error_config[gt_error]["locus_length"])]
            
            if stat_type.startswith("pvalue"):
                for i, row in df_cur.iterrows():
                    df_cur.at[i,stat_type] = "%.5f" % row[stat_type]
                print(scale, gt_error)
                # print(df_cur[stat_type])
                axes[scale_idx][gt_error].set_xlim(-0.05,1.05)
                sns.histplot(x=stat_type,
                        binwidth=0.05,
                        binrange=(0,1),
                        data=df_cur, 
                        ax=axes[scale_idx][gt_error], 
                        color='#607c8e',
                        edgecolor="white",
                        stat="probability"
                        )
                axes[scale_idx][gt_error].set_ylim(0, 1.0)
                significance_count = df_cur.loc[df_cur[stat_type] < 0.05,stat_type].count()/df_cur[stat_type].count()
                print(f"{stat_type} -- Proportion of significance:{significance_count}")
            else:
                MAX = 210
                x = np.arange(0, MAX, .001)
                axes[scale_idx][gt_error].plot(x, scipy.stats.chi2.pdf(x, df=DOF), color='darkorange', lw=1)
                bins = np.arange(0, MAX+1, 10)
                sns.histplot(x=np.clip(df_cur[stat_type], bins[0], bins[-1]),
                            bins=bins,
                            data=df_cur, 
                            ax=axes[scale_idx][gt_error], 
                            color='#607c8e',
                            edgecolor="white",
                            stat="density"
                            )
                axes[scale_idx][gt_error].set_xlabel("", fontsize=16)
                axes[scale_idx][gt_error].set_ylabel("", fontsize=16)
                axes[scale_idx][gt_error].tick_params(axis='y', labelsize=12)
                axes[scale_idx][gt_error].tick_params(axis="x", labelsize=12)
                axes[scale_idx][gt_error].set_xlim(0, MAX)
                axes[scale_idx][gt_error].set_ylim(0, 0.11)

                # xlabels = bins[1:].astype(str)
                # xlabels[-1] += '+'
                N_labels = len(bins)/2-1
                ticks = (20 * np.arange(N_labels) + 20).astype(int)
                xlabels = ticks.astype(str)
                xlabels[-1] += '+'
                axes[scale_idx][gt_error].set_xticks(ticks)
                axes[scale_idx][gt_error].set_xticklabels(xlabels)
            axes[0][gt_error].set_title(gt_error_name[gt_error], size=20, fontweight="bold")
            
        axes[scale_idx][0].set_ylabel(scale_list[scale_idx], size=16)
    
    plt.tight_layout()
    plt.savefig(figpath, dpi=300)

if __name__ == "__main__":
    path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net1_2/statistics_0.csv"
    # path = sys.argv[1]
    stats_list = ["chisq", "pvalue_chisq", "gtest", "pvalue_gtest", "pvalue_exact_mc"]
    num_reti = 0
    # stats_list = ["chisq", "gtest"]
    # stats_list = ["pvalue_chisq", "pvalue_gtest", "pvalue_exact_mc"]
    for stat_type in stats_list:
        plot(path, stat_type, num_reti)
    # path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/exact_mn.csv"
    # plot(path, "exact_mn_mc")