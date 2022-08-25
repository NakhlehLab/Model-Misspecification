import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import norm

OUTGROUP = "Z_0"
HYBRID_SPECIES = "Q_0"

def plot_D3(path, fig_dir, D3_type=3):
    if D3_type == 1:
        figpath = fig_dir + "/D3.pdf"
    elif D3_type == 2:
        figpath = fig_dir + "/D3_gt.pdf"
    else:
        figpath = fig_dir + "/D3_2stage_faster.pdf"
    df = pd.read_csv(path)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    df = df[(df["A"] != OUTGROUP) & (df["B"] != OUTGROUP) & (df["C"] != OUTGROUP)]
    
    for j, locus_length in enumerate([2000, 500]):
        sns.violinplot(y="D3", x="scale", data=df[(df["locus_length"] == locus_length)], ax = axes[j])
        # axes[j].set_ylim(-0.08, 0.08)
        # axes[j].set_ylim(-0.011, 0.011)
        # axes[j].yaxis.set_ticks([-0.01, -0.005, 0, 0.005, 0.01])
        axes[j].set_xlabel("", fontsize=16)
        axes[j].set_ylabel("$D_3$", fontsize=16)
        axes[j].tick_params(axis='y', labelsize=12)
        axes[j].tick_params(axis="x", labelsize=12)
        axes[j].set_title(f"locus_length={locus_length}", size=16, fontweight="bold")
        print(f"locus_length={locus_length}")
        print(df[(df["locus_length"] == locus_length)])
    plt.tight_layout()
    plt.savefig(figpath)
    plt.show()

def plot_pvalue(path, fig_dir, D3_type=3):
    if D3_type == 1:
        figpath = fig_dir + "/D3_pval.pdf"
    elif D3_type == 2:
        figpath = fig_dir + "/D3_gt_pval.pdf"
    else:
        figpath = fig_dir + "/D3_2stage.pdf"
    df = pd.read_csv(path)
    df = df[(df["A"] != OUTGROUP) & (df["B"] != OUTGROUP) & (df["C"] != OUTGROUP)]
    fig, axes = plt.subplots(3, 2, figsize=(12, 7))
    # pval = []
    # for i, row in df.iterrows():
    #     pval.append(get_pval(row["D3"], row["D3_stdev"]))
    # df["pval"] = pval    
    
    for j, locus_length in enumerate([2000, 500]):
        for i, scale in enumerate(["short", "medium", "long"]):
        # for i, scale in enumerate(["short"]):
            df_cur = df[(df["locus_length"] == locus_length) & (df["scale"] == scale)]
            sns.histplot(x="D3_pval", data=df_cur, binwidth=0.05, ax = axes[i][j], color='#607c8e',
                    edgecolor="white",
                    stat="probability")
            significance_count = df_cur.loc[df_cur["D3_pval"] < 0.05,"D3_pval"].count()/df_cur["D3_pval"].count()
            print(f"{significance_count}")
            axes[i][j].set_xlim(0, 1)
            axes[i][j].set_ylim(0, 1)
            axes[i][j].set_xlabel("P-value", fontsize=16)
            axes[i][j].set_ylabel("")
            axes[i][0].set_ylabel(scale, fontsize=16, fontweight="bold")
            axes[i][j].tick_params(axis='y', labelsize=12)
            axes[i][j].tick_params(axis="x", labelsize=12)
            axes[i][j].axhline(y=0.05, color='r', linestyle='-')
        axes[0][j].set_title(f"locus_length={locus_length}", size=16, fontweight="bold")
    plt.tight_layout()
    plt.savefig(figpath)
    plt.show()

def plot_sep_pval(path, fig_dir, marker=True):
    if marker:
        figpath = fig_dir + "/D3_pval_sep.pdf"
    else:
        figpath = fig_dir + "/D3_gt_pval_sep.pdf"
    fig, axes = plt.subplots(3, 2, figsize=(12, 7))
    df = pd.read_csv(path)
    df1 = df[(df["A"] != HYBRID_SPECIES) & (df["B"] != HYBRID_SPECIES) & (df["C"] != HYBRID_SPECIES)]
    df2 = df[(df["A"] == HYBRID_SPECIES) | (df["B"] == HYBRID_SPECIES) | (df["C"] == HYBRID_SPECIES)]
    df1["net"] = [False for i in range(len(df1["id"]))]
    df2["net"] = [True for i in range(len(df2["id"]))]
    df = pd.concat([df1, df2])
    # print(df_cur)
    for j, locus_length in enumerate([2000, 500]):
        for i, scale in enumerate(["short", "medium", "long"]):
            df_cur = df[(df["scale"] == scale) & (df["locus_length"] == locus_length)]
            sns.histplot(x="D3_pval", data=df_cur, binwidth=0.05, ax = axes[i][j], color='#607c8e',
                    edgecolor="white",
                    stat="probability",
                    hue = "net",
                    multiple="stack")
            axes[i][j].set_xlim(0, 1)
            axes[i][j].set_ylim(0, 1)
            axes[i][j].set_xlabel("P-value", fontsize=16)
            axes[i][j].set_ylabel("")
            axes[i][0].set_ylabel(scale, fontsize=16, fontweight="bold")
            axes[i][j].tick_params(axis='y', labelsize=12)
            axes[i][j].tick_params(axis="x", labelsize=12)
            df_cur_net = df_cur[df_cur["net"] == True]
            df_cur_tree = df_cur[df_cur["net"] == False]
            significance_count_net = df_cur_net.loc[df_cur_net["D3_pval"] < 0.05,"D3_pval"].count()/df_cur["D3_pval"].count()
            significance_count_tree = df_cur_tree.loc[df_cur_tree["D3_pval"] < 0.05,"D3_pval"].count()/df_cur["D3_pval"].count()
            print(f"scale={scale}, locus_length={locus_length} -- Proportion of significance:{significance_count_net}, {significance_count_tree},")
            if i == 0 and j == 0:
                labels = axes[0][0].get_legend()
                handles = labels.legendHandles
            axes[i][j].get_legend().remove()
            axes[i][j].axhline(y=0.05, color='r', linestyle='-')
        axes[0][j].set_title(f"locus_length={locus_length}", size=16, fontweight="bold")
    # handles, labels = axes[0][0].get_legend_handles_labels()
    
    print(handles, labels)
    # label_map = {True: "net", False: "tree"}
    # new_labels = [label_map[x] for x in labels[0:2]]

    fig.legend(handles, ["tree", "net"], loc='lower center', fontsize=16, ncol=2)
    plt.tight_layout()
    fig.subplots_adjust(bottom=0.15, hspace=0.15)
    plt.savefig(figpath, dpi=300)
    plt.show()


if __name__ == "__main__":
    path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/D3_res_2stage_faster.csv"
    fig_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/"
    D3_type=3
    # path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net1/D3_res.csv"
    # fig_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net1/"
    # path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/test/net1/D3_res.csv"
    # fig_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/test/net1/"
    plot_D3(path, fig_dir, D3_type)
    plot_pvalue(path, fig_dir, D3_type)

    # path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net1/D3_gt_res.csv"
    # fig_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net1/"
    # plot_D3(path, fig_dir, marker)
    # plot_sep_pval(path, fig_dir, marker)

