import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from pathlib import Path
import numpy as np
import scipy.stats

DOF = 1

def plot_stats(path, stat_type, num_reti):
    figpath = str(Path(path).parent.absolute())+ "/fig_triplet_"+stat_type+"_"+str(num_reti)+".pdf"
    df = pd.read_csv(path)
    # fig, axes = plt.subplots(3, 4, figsize=(20, 10))
    fig, axes = plt.subplots(3, 3, figsize=(18, 10))

    
    # gt_error_name = {0:"true st true gt", 1:"true gt", 2: "inferred gt from 1000 bps", 3: "inferred gt from 500 bps"}
    # gt_error_config = {0:{"true_st": True, "true_gt": True, "locus_length": 500}, 1:{"true_st": False, "true_gt": True, "locus_length": 500}, 2:{"true_st": False, "true_gt": False, "locus_length": 1000}, 3:{"true_st": False, "true_gt": False, "locus_length": 500} }
    gt_error_name = {0:"True gt", 1: "Inferred gt from 2000 bps", 2: "Inferred gt from 500 bps"}
    gt_error_config = {0:{"true_st": False, "true_gt": True, "locus_length": 500}, 1:{"true_st": False, "true_gt": False, "locus_length": 2000}, 2:{"true_st": False, "true_gt": False, "locus_length": 500} }

    scale_list = ["short", "medium", "long"]
    for scale_idx in [0, 1, 2]:
        scale = scale_list[scale_idx]
        for gt_error in [0, 1, 2]:
            df_cur = df[(df["scale"] == scale)
            & (df["true_gt"] == gt_error_config[gt_error]["true_gt"])
            & (df["locus_length"] == gt_error_config[gt_error]["locus_length"])]
            
            if stat_type.startswith("pvalue"):
                for i, row in df_cur.iterrows():
                    df_cur.at[i,stat_type] = "%.5f" % row[stat_type]
                print(scale, gt_error)
                print(df_cur[stat_type])
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
            else:
                x = np.arange(0, 10, .01)
                axes[scale_idx][gt_error].plot(x, scipy.stats.chi2.pdf(x, df=DOF), color='darkorange', lw=1)
                sns.histplot(x=stat_type,
                            bins=20,
                            data=df_cur, 
                            ax=axes[scale_idx][gt_error], 
                            color='#607c8e',
                            edgecolor="white",
                            stat="probability"
                            )
            axes[scale_idx][gt_error].set_xlabel("", fontsize=16)
            axes[scale_idx][gt_error].set_ylabel("", fontsize=16)
            axes[scale_idx][gt_error].tick_params(axis='y', labelsize=12)
            axes[scale_idx][gt_error].tick_params(axis="x", labelsize=12)
            axes[0][gt_error].set_title(gt_error_name[gt_error], size=20, fontweight="bold")
        axes[scale_idx][0].set_ylabel(scale_list[scale_idx], size=16)
    
    plt.tight_layout()
    plt.savefig(figpath, dpi=300)


def plot_multitest(path, output_path, output_fig):
    df = pd.read_csv(path)
    method_list = ["bonferroni", "simes", "fisher", "fdr_bh", "fdr_tsbh"]
    pvalue_list = ["pvalue_chisq", "pvalue_gtest", "pvalue_exact"]
    gt_error_name = {0:"True gt", 1: "Inferred gt from 2000 bps", 2: "Inferred gt from 500 bps"}


    res_dict = {"scale":[],"gt_error":[], "pvalue_type":[],"method":[],"num_sig":[]}
    gt_error_config = {0:{"true_gt": True, "locus_length": 500}, 1:{"true_gt": False, "locus_length": 2000}, 2:{"true_gt": False, "locus_length": 500} }
    scale_list = ["short", "medium", "long"]
    for scale_idx in [0, 1, 2]:
        scale = scale_list[scale_idx]
        for gt_error in [0, 1, 2]:
            df_cur = df[(df["scale"] == scale) 
            & (df["true_gt"] == gt_error_config[gt_error]["true_gt"])
            & (df["locus_length"] == gt_error_config[gt_error]["locus_length"])]
            for method in method_list:
                for pvalue_type in pvalue_list:
                    print(method)
                    df_m = df_cur[(df_cur["method"] == method) &(df_cur["pvalue_type"] == pvalue_type)]
                    # print(df_m)
                    num_sig = sum(df_m["significant"])/len(df_m["significant"])
                    res_dict["scale"].append(scale)
                    res_dict["gt_error"].append(gt_error)
                    res_dict["pvalue_type"].append(pvalue_type)
                    res_dict["method"].append(method)
                    res_dict["num_sig"].append(num_sig)
    res_df = pd.DataFrame(res_dict)
    print(res_df)
    res_df.to_csv(output_path)
    
    fig, axes = plt.subplots(3, 3, figsize=(18, 10))
    for scale_idx in [0, 1, 2]:
        scale = scale_list[scale_idx]
        for gt_error in [0, 1, 2]:
            df_cur = res_df[(res_df["scale"] == scale) 
            & (res_df["gt_error"] == gt_error)]
            sns.barplot(x="method", y="num_sig", hue="pvalue_type", palette="husl",
            data=df_cur, ax=axes[scale_idx][gt_error],)

            axes[scale_idx][gt_error].set_ylim(0, 1.0)
            axes[scale_idx][gt_error].set_xlabel("", fontsize=16)
            axes[scale_idx][gt_error].set_ylabel("", fontsize=16)
            axes[scale_idx][gt_error].tick_params(axis='y', labelsize=12)
            axes[scale_idx][gt_error].tick_params(axis="x", labelsize=12)
            axes[scale_idx][gt_error].get_legend().remove()
            axes[scale_idx][gt_error].axhline(y=0.05, color='r', linestyle='-')
            axes[0][gt_error].set_title(gt_error_name[gt_error], size=20, fontweight="bold")
            
        axes[scale_idx][0].set_ylabel(scale_list[scale_idx], size=16)
    
    handles, labels = axes[0][0].get_legend_handles_labels()
    label_map = {"pvalue_chisq": "chisquare pvalue","pvalue_gtest": "G-test pvalue", "pvalue_exact": "Exact pvalue" }
    new_labels = [label_map[x] for x in labels]
    print(handles, labels)
    fig.legend(handles, new_labels, loc='lower center', fontsize=20, ncol=3)
    plt.savefig(output_fig)
        

def plot_multitest_D3(path_multi, path_D3_marker, path_D3_2stage, path_D3_gt, output_path, output_fig):
    df_multi = pd.read_csv(path_multi)
    df_D3_marker = pd.read_csv(path_D3_marker)
    df_D3_2stage = pd.read_csv(path_D3_2stage)
    df_D3_gt = pd.read_csv(path_D3_gt)

    # using markers
    df_D3_marker["true_gt"] = [False for x in range(df_D3_marker.shape[0])]
    df_D3_2stage["true_gt"] = [False for x in range(df_D3_2stage.shape[0])]

    #using gene trees
    df_D3_gt["true_gt"] = [True for x in range(df_D3_gt.shape[0])]


    df_D3_marker['pvalue_type'] = df_D3_marker['pvalue_type'].map({'D3_pval':'D3_pval_1',
                             np.nan:'NY'},
                             na_action=None)
    df_D3_2stage['pvalue_type'] = df_D3_2stage['pvalue_type'].map({'D3_pval':'D3_pval_2',
                             np.nan:'NY'},
                             na_action=None)

    df_D3_gt['pvalue_type'] = df_D3_gt['pvalue_type'].map({'D3_pval':'D3_pval_1',
                             np.nan:'NY'},
                             na_action=None)
    
    df = pd.concat([df_multi, df_D3_marker, df_D3_2stage, df_D3_gt])

    method_list = ["bonferroni", "simes-hochberg", "fdr_bh"]
    pvalue_list = ["pvalue_chisq", "pvalue_gtest", "pvalue_exact", "D3_pval_1", "D3_pval_2"]

    gt_error_name = {0:"None (true GTs)", 1: "Low (2000bp)", 2: "High (500bp)"}
    res_dict = {"scale":[],"gt_error":[], "pvalue_type":[],"method":[],"num_sig":[]}
    gt_error_config = {0:{"true_gt": True}, 1:{"true_gt": False, "locus_length": 2000}, 2:{"true_gt": False, "locus_length": 500} }
    scale_list = ["short", "medium", "long"]
    for scale_idx in range(3):
        scale = scale_list[scale_idx]
        for gt_error in range(3):
            df_cur = df[(df["scale"] == scale) 
            & (df["true_gt"] == gt_error_config[gt_error]["true_gt"])]
            for method in method_list:
                for pvalue_type in pvalue_list:
                    # if pvalue_type == "D3_pval" and gt_error_config[gt_error]["true_gt"] == True:
                    #     continue
                    print(f"{scale, gt_error,method, pvalue_type}")
                    if gt_error != 0:
                        df_m = df_cur[(df_cur["method"] == method) &(df_cur["pvalue_type"] == pvalue_type) & (df_cur["locus_length"] == gt_error_config[gt_error]["locus_length"])]
                    else:
                        df_m = df_cur[(df_cur["method"] == method) &(df_cur["pvalue_type"] == pvalue_type)]
                    if "D3_pval_2" in pvalue_type and gt_error == 0:
                        continue
                    num_sig = sum(df_m["significant"])/len(df_m["significant"])
                    # print(f"{scale, gt_error,method, pvalue_type, num_sig}")
                    # if "D3" in pvalue_type:
                    
                    #     print(sum(df_m["significant"]))
                    #     print(df_m)
                    res_dict["scale"].append(scale)
                    res_dict["gt_error"].append(gt_error)
                    res_dict["pvalue_type"].append(pvalue_type)
                    res_dict["method"].append(method)
                    res_dict["num_sig"].append(num_sig)
    res_df = pd.DataFrame(res_dict)
    # print(res_df)
    res_df.to_csv(output_path)
    
    sns.set_context(rc = {'patch.linewidth': 1.0})
    palette=sns.color_palette("colorblind")

    y_label_list = ["High (1x CU)", "Moderate (2x CU)", "Low (5x CU)"]
    fig, axes = plt.subplots(3, 3, figsize=(18, 10))
    for scale_idx in range(3):
        scale = scale_list[scale_idx]
        for gt_error in range(3):
            df_cur = res_df[(res_df["scale"] == scale) 
            & (res_df["gt_error"] == gt_error)]
            sns.barplot(x="method", y="num_sig", hue="pvalue_type", palette=palette,
            data=df_cur, ax=axes[scale_idx][gt_error], edgecolor="white")
            # print(df_cur[df_cur["pvalue_type"] != "D3_pval"])

            axes[scale_idx][gt_error].set_ylim(0, 1.0)
            axes[scale_idx][gt_error].set_xlabel("", fontsize=16)
            axes[scale_idx][gt_error].set_ylabel("", fontsize=16)
            axes[scale_idx][gt_error].tick_params(axis='y', labelsize=12)
            axes[scale_idx][gt_error].tick_params(axis="x", labelsize=12)
            axes[scale_idx][gt_error].get_legend().remove()
            axes[scale_idx][gt_error].axhline(y=0.05, color='r', linestyle='-')
            axes[0][gt_error].set_title(gt_error_name[gt_error], size=16)
            labels = [item.get_text() for item in axes[scale_idx][gt_error].get_xticklabels()]
            axes[scale_idx][gt_error].set_xticklabels(["Bonferroni\n(FWER)", "Simes-Hochberg\n(FWER)","Benjamini-Hochberg\n(FDR)"])

            print(labels)
        axes[scale_idx][0].set_ylabel(y_label_list[scale_idx], size=16)
        

    handles, labels = axes[1][1].get_legend_handles_labels()
    label_map = {"pvalue_chisq": "Pearson's","pvalue_gtest": "G", "pvalue_exact": "Exact multinomial", "D3_pval_1": "Block bootstrapped D3", "D3_pval_2": "2-stage bootstrapped D3"}
    new_labels = [label_map[x] for x in labels]
    legend = fig.legend(handles, new_labels, loc='lower center', title="Statistical test", fontsize=16, ncol=5)
    legend.get_title().set_fontsize('16') 
    # fig.text(0.055, 0.03, 'Statistical test', ha='center', va='center', fontsize=16)
    fig.text(0.02, 0.51, 'Level of ILS', ha='center', va='center', rotation='vertical', fontsize=18, fontweight="bold")

    
    fig.suptitle('Gene tree error', fontsize=18, fontweight="bold")
    fig.subplots_adjust(top=0.88, )
    fig.tight_layout()
    
    fig.subplots_adjust(left=0.08, bottom=0.15, hspace=0.2)
    
    plt.savefig(output_fig)
    # plt.show()


if __name__ == "__main__":
    exp_type = True
    num_reti = "0"
    dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net1_2/"
    path = dir_path+"multitest_"+num_reti+"_"+str(exp_type)+".csv"

    D3_path_marker = dir_path+"D3_multitest.csv"
    D3_path_gt = dir_path+"D3_gt_multitest.csv"
    D3_2stage = dir_path+"D3_2stage_multitest.csv"

    output_path = dir_path+"multitest_numsig_"+str(exp_type)+"_"+num_reti+"2.csv"
    output_fig = dir_path+"multitest_numsig_"+str(exp_type)+"_"+num_reti+"2.pdf"
    # path = sys.argv[1]
    # stats_list = ["chisq", "pvalue_chisq", "gtest", "pvalue_gtest", "pvalue_exact_mc"]
    # pvalue_list = ["pvalue_chisq", "pvalue_gtest", "pvalue_exact"]
    
    # for stat_type in pvalue_list:
    # plot_multitest(path, output_path, output_fig)
    plot_multitest_D3(path, D3_path_marker, D3_2stage, D3_path_gt, output_path, output_fig)

    # path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/exact_mn.csv"
    # plot(path, "exact_mn_mc")

    # path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net1/triplet_stats_all_expTrue_0.csv"
    # stats_list = ["chisq", "pvalue_chisq", "gtest", "pvalue_gtest", "pvalue_exact"]
    # num_reti = 0
    # stats_list = ["chisq", "gtest"]
    # for stat_type in stats_list:
    #     plot_stats(path, stat_type, num_reti)