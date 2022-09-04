import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from pathlib import Path
import numpy as np
import scipy.stats


def plot_figure():
    res_dict = {"scale":["short", "short", "short", "short", "short", "short", "short", "short","short", "short", "short", "short",
                        "medium", "medium", "medium", "medium", "medium", "medium", "medium", "medium", "medium", "medium", "medium", "medium",
                        "long", "long", "long", "long", "long", "long", "long", "long", "long", "long", "long", "long"],
                "gt_error":[0, 0, 0, 0, 1,1,1,1, 2,2,2,2,
                            0, 0, 0, 0, 1,1,1,1, 2,2,2,2,
                            0, 0, 0, 0, 1,1,1,1, 2,2,2,2], 
                "pvalue_type":["pvalue_chisq", "pvalue_gtest", "pvalue_exact", "Cai", "pvalue_chisq", "pvalue_gtest", "pvalue_exact", "Cai","pvalue_chisq", "pvalue_gtest", "pvalue_exact", "Cai",
                               "pvalue_chisq", "pvalue_gtest", "pvalue_exact", "Cai", "pvalue_chisq", "pvalue_gtest", "pvalue_exact", "Cai","pvalue_chisq", "pvalue_gtest", "pvalue_exact", "Cai",
                               "pvalue_chisq", "pvalue_gtest", "pvalue_exact", "Cai", "pvalue_chisq", "pvalue_gtest", "pvalue_exact", "Cai","pvalue_chisq", "pvalue_gtest", "pvalue_exact", "Cai" ],
                "num_sig":[0.05, 0.05, 0.03, 0, 0.59, 0.5, 0.5, 0, 1,1,1,0.02,
                           0.1, 0.14, 0.06, 0.06, 1.0, 1.0, 1, 0.02, 1, 1, 1, 0.08,
                           0.38, 0, 0.02, 0.05, 1.0, 1.0, 1.0, 0.37, 1.0, 1.0, 1.0,1 ]}
    for key in res_dict.keys():
        print(f"{key, len(res_dict[key])}")
    # print(len(res_dict[""]))
    # gt_error_config = {0:{"true_gt": True, "locus_length": 500}, 1:{"true_gt": False, "locus_length": 2000}, 2:{"true_gt": False, "locus_length": 500} }
    scale_list = ["short", "medium", "long"]
    pvalue_list = ["pvalue_chisq", "pvalue_gtest", "pvalue_exact"]
    

    gt_error_name = {0:"True GTs", 1: "Low error (2000 bps)", 2: "High error (500 bps)"}
    res_df = pd.DataFrame(res_dict)
    print(res_df)

    y_label_list = ["High (1x CU)", "Moderate (2x CU)", "Low (5x CU)"]
    # res_df.to_csv(output_path)
    palette=sns.color_palette("colorblind")
    fig, axes = plt.subplots(3, 1, figsize=(6, 10))
    for scale_idx in [0, 1, 2]:
        scale = scale_list[scale_idx]
        # for gt_error in [0, 1, 2]:
        #     df_cur = res_df[(res_df["scale"] == scale) 
        #     & (res_df["gt_error"] == gt_error)]
        #     sns.barplot(x="pvalue_type", y="num_sig", palette=palette,

        #     # sns.barplot(x="method", y="num_sig", hue="pvalue_type", palette="husl",
        #     data=df_cur, ax=axes[scale_idx][gt_error],)

        #     axes[scale_idx][gt_error].set_ylim(0, 1.0)
        #     axes[scale_idx][gt_error].set_xlabel("", fontsize=16)
        #     axes[scale_idx][gt_error].set_ylabel("", fontsize=16)
        #     axes[scale_idx][gt_error].tick_params(axis='y', labelsize=12)
        #     axes[scale_idx][gt_error].tick_params(axis="x", labelsize=12)
        #     # axes[scale_idx][gt_error].get_legend().remove()
        #     axes[scale_idx][gt_error].axhline(y=0.05, color='r', linestyle='-')
        #     axes[0][gt_error].set_title(gt_error_name[gt_error], size=20, fontweight="bold")
        df_cur = res_df[(res_df["scale"] == scale)]
        sns.barplot(x="gt_error", y="num_sig",hue="pvalue_type", palette=palette, edgecolor="white",
        data=df_cur, ax=axes[scale_idx],)
        # sns.barplot(x="method", y="num_sig", hue="pvalue_type", palette="husl",
        

        axes[scale_idx].set_ylim(0, 1.0)
        axes[scale_idx].set_xlabel("", fontsize=16)
        axes[scale_idx].set_ylabel("", fontsize=16)
        axes[scale_idx].tick_params(axis='y', labelsize=12)
        axes[scale_idx].tick_params(axis="x", labelsize=12)
        axes[scale_idx].get_legend().remove()
        axes[scale_idx].axhline(y=0.05, color='r', linestyle='-')
        labels = [item.get_text() for item in axes[scale_idx].get_xticklabels()]
        axes[scale_idx].set_xticklabels(["None (true GTs)", "Low (2000bp)","High (500bp)"])
        # axes[0].set_title(gt_error_name[gt_error], size=20, fontweight="bold")
            
        axes[scale_idx].set_ylabel(y_label_list[scale_idx], size=16)
    
    handles, labels = axes[0].get_legend_handles_labels()
    label_map = {"pvalue_chisq": "Pearson's","pvalue_gtest": "G", "pvalue_exact": "Approximate multinomial", "Cai": "Cai & Ane" }
    new_labels = [label_map[x] for x in labels]
    legend = fig.legend(handles, new_labels, loc='lower center', title="Statistical test", fontsize=16, ncol=1)
    legend.get_title().set_fontsize('16') 
    fig.text(0.02, 0.51, 'Level of ILS', ha='center', va='center', rotation='vertical', fontsize=18, fontweight="bold")
    fig.text(0.35, 0.195, 'Gene tree error', fontsize=18, fontweight="bold")
    fig.subplots_adjust(top=0.95, left=0.2, bottom=0.25, hspace=0.24)
    plt.savefig("/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/teststat.pdf")
if __name__ == "__main__":
    plot_figure()