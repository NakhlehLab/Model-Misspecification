import pandas as pd
import dendropy
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from pathlib import Path


SPECIES = ["C 0", "G 0", "L 0", "Q 0", "R 0", "Z 0"]


def compute_simulated_gt_outgroup_brs(tree_path, murate_path, outgroup = "Z 0"):
    tree_list = dendropy.TreeList.get(path=tree_path, schema="newick", rooting="force-rooted")
    murate_list = []
    if murate_path != None:
        with open(murate_path, "r") as handle:
            for line in handle.readlines():
                murate_list.append(float(line.strip()))
    dist = []
    if len(murate_list) == 0:
        for tree in tree_list:
            tree.encode_bipartitions()
            mrca = tree.mrca(taxon_labels=list(set(SPECIES)-set([outgroup])))
            dist.append(mrca.edge.length)
            
    else:
        for lcs, tree in enumerate(tree_list):
            tree.encode_bipartitions()
            mrca = tree.mrca(taxon_labels=list(set(SPECIES)-set([outgroup])))
            if mrca != tree.seed_node:
                # print(mrca.edge.length)
                dist.append(mrca.edge.length * murate_list[lcs])
                # print(dist[-1])
    return dist


def run_all_branches(dir_path, num_replicate, THETA_DICT, coalescentunit = True, longer=True, log=True):
    df_list = []
    if longer:
        locus_length_list = [500, 2000]
    else:
        locus_length_list = [200, 1000]
    

    for scale in ["short", "medium", "long"]:
        for locus_length in locus_length_list:

            for i in range(1, num_replicate + 1):
                path = dir_path +scale+"/"+ str(i) + "/homo/"+str(locus_length) + "/genetrees.txt"
                murate_path = None
                if coalescentunit == False:
                    murate_path = dir_path +scale+"/"+ str(i) + "/heter/"+str(locus_length) + "/murate.txt"
                
                dist = compute_simulated_gt_outgroup_brs(path, murate_path)
                df = pd.DataFrame()
                if coalescentunit == False:
                    df["branch_lengths"] = [x*THETA_DICT[scale] for x in dist]
                else:
                    df["branch_lengths"] =  dist
                df["scale"] = [scale for x in range(len(df["branch_lengths"]))]
                df["locus_length"] = [locus_length for x in range(len(df["branch_lengths"]))]
                df_list.append(df)
    res_df = pd.concat(df_list)
    if coalescentunit:
        res_df.to_csv(dir_path+"/gt_branches_out2in_distribution_truegt_coal.csv")
    else:
        res_df.to_csv(dir_path+"/gt_branches_out2in_distribution_truegt.csv")
    
    print(res_df)
    plt.clf()
    ax = sns.violinplot(x="scale", y="branch_lengths", data=res_df, hue="locus_length")
    log_dict={True: "_log", False:""}
    if log:
        ax.set_yscale('log')

    plt.grid(axis = 'y')
    if coalescentunit:
        plt.savefig(dir_path+"/gt_branches_out2in_truegt_distribution_coal_"+log_dict[log]+".pdf")
    else:
        plt.savefig(dir_path+"/gt_branches_out2in_truegt_distribution_"+log_dict[log]+".pdf")
    

def plot_compare_branches_cdf(dir_path, longer):
    plt.clf()
    fig, axes = plt.subplots(1, 2, figsize=(15,5))
    if longer:
        locus_length_list = [500, 2000]
    else:
        locus_length_list = [200, 1000]

    df = pd.read_csv(dir_path+"/gt_branches_out2in_distribution_truegt.csv", index_col=False)
    # iqgt_df = pd.read_csv(dir_path+"/gt_branches_distribution_iqgt.csv", index_col=False)

    df["truegt"] = [True for x in range(len(df["branch_lengths"]))]
    # iqgt_df["truegt"] = [False for x in range(len(iqgt_df["branch_lengths"]))]
    
    # df = pd.concat(truegt_df, ignore_index=True)
    

    for idx, locus_length in enumerate(locus_length_list):
        print(df[df["locus_length"] == locus_length])
        hue = df[['scale', 'truegt']].apply(
        lambda row: f"{row.scale}, {row.truegt}", axis=1)
        hue.name = 'scale, truegt'
        sns.ecdfplot(df[df["locus_length"] == locus_length], x="branch_lengths", hue=hue, stat="proportion", ax=axes[idx], linewidth=2, alpha=0.6)
        
        # lss = ['--','--','--','-','-', '-']
        # colors = ["blue", "orange", "green", "blue", "orange", "green"]
        # # labels = axes[idx].get_legend()
        # # handles = labels.legendHandles
        # handles = axes[idx].legend_.legendHandles[::-1]

        # for line, ls, handle, clr in zip(axes[idx].lines, lss, handles, colors):
        #     line.set_linestyle(ls)
        #     line.set_color(clr)
        #     handle.set_ls(ls)
        #     handle.set_color(clr)

        axes[idx].tick_params(axis='y', labelsize=14)
        axes[idx].tick_params(axis='x', labelsize=14)
        axes[idx].set_ylim(0, 1.03)
        # axes[idx].set_xscale('log')
        axes[idx].set_xlim(0.0001, 2.3)
        axes[idx].set_xlabel("Branch lengths in substitution rates", fontsize=18)
        axes[idx].set_ylabel("Proportion of branches", fontsize=18)

    axes[0].set_title("locus length="+str(locus_length_list[0]), size=20, fontweight="bold")
    axes[1].set_title("locus length="+str(locus_length_list[1]), size=20, fontweight="bold")
    
    # labels = axes[0][0].get_legend()
    # handles = labels.legendHandles
    # print(handles, labels)
    # fig.legend(handles, new_labels, loc='lower center', fontsize=20, ncol=3)
    plt.savefig(dir_path+"/gt_branches_out2in_distribution_compare_cdf.pdf")


def plot_out_branches_compare_cdf(dir_path1, dir_path2, coal=True):
    plt.clf()
    fig, axes = plt.subplots(1, 2, figsize=(15,5))
    
    locus_length_list1 = [500, 2000]
    locus_length_list2 = [200, 1000]

    if coal:
        df1 = pd.read_csv(dir_path1+"/gt_branches_out2in_distribution_truegt_coal.csv", index_col=False)
        df2 = pd.read_csv(dir_path2+"/gt_branches_out2in_distribution_truegt_coal.csv", index_col=False)
        
    else:
        df1 = pd.read_csv(dir_path1+"/gt_branches_out2in_distribution_truegt.csv", index_col=False)
        df2 = pd.read_csv(dir_path2+"/gt_branches_out2in_distribution_truegt.csv", index_col=False)

    # df1["truegt"] = [True for x in range(len(df["branch_lengths"]))]
    # iqgt_df["truegt"] = [False for x in range(len(iqgt_df["branch_lengths"]))]
    
    df = pd.concat([df1, df2], ignore_index=True)
    

    for idx in range(2):
        # print(df[df["locus_length"] == locus_length])
        df_cur = df[(df["locus_length"] == locus_length_list1[idx]) | (df["locus_length"] == locus_length_list2[idx])]
        hue = df_cur[['scale', 'locus_length']].apply(
        lambda row: f"{row.scale}, {row.locus_length}", axis=1)
        hue.name = 'scale, locus_length'
        sns.ecdfplot(df_cur, x="branch_lengths", hue=hue, stat="proportion", ax=axes[idx], linewidth=2, alpha=0.6)
        
        lss = [':',':',':','-','-', '-']
        # marker=[""]
        colors = ["blue", "orange", "green", "blue", "orange", "green"]
        # labels = axes[idx].get_legend()
        # handles = labels.legendHandles
        handles = axes[idx].legend_.legendHandles[::-1]

        for line, ls, handle, clr in zip(axes[idx].lines, lss, handles, colors):
            line.set_linestyle(ls)
            line.set_color(clr)
            # line.set_marker("o")
            handle.set_ls(ls)
            handle.set_color(clr)

        axes[idx].tick_params(axis='y', labelsize=14)
        axes[idx].tick_params(axis='x', labelsize=14)
        
        if coal:
            axes[idx].set_xlabel("Branch lengths in coalescent unit", fontsize=18)
            axes[idx].set_ylim(0, 1.03)
            axes[idx].set_xlim(5, 14)
            # axes[idx].set_xscale('log')
        else:
            axes[idx].set_xlabel("Branch lengths in substitution rates", fontsize=18)
            axes[idx].set_ylim(0, 1.03)
            axes[idx].set_xlim(0.0001, 2.3)
        axes[idx].set_ylabel("Proportion of branches", fontsize=18)

    axes[0].set_title("locus length short", size=20, fontweight="bold")
    axes[1].set_title("locus length long", size=20, fontweight="bold")
    
    # labels = axes[0][0].get_legend()
    # handles = labels.legendHandles
    # print(handles, labels)
    # fig.legend(handles, new_labels, loc='lower center', fontsize=20, ncol=3)
    if coal:
        plt.savefig(str(Path(dir_path1).parent.parent.absolute())+"/gt_branches_out2in_distribution_compare_coal_cdf.pdf")
    else:
        plt.savefig(str(Path(dir_path1).parent.parent.absolute())+"/gt_branches_out2in_distribution_compare_cdf.pdf")


if __name__ == "__main__":

    # dir_path = sys.argv[1]
    # num_replicate = int(sys.argv[2])
    # longer = (sys.argv[3].lower() == "true")
    # log = (sys.argv[4].lower() == "true")
    # larger_theta = (sys.argv[5].lower() == "true")
    # coalescentunit = (sys.argv[6].lower() == "true")

    # if larger_theta == False:
    #     THETA_DICT = {"short":0.025, "medium":0.0125, "long":0.005}
    # else:
    #     THETA_DICT = {"short":0.05, "medium":0.025, "long":0.01}

    # run_all_branches(dir_path, num_replicate, THETA_DICT, coalescentunit, longer, log)

    dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/"
    # plot_compare_branches_cdf(dir_path, True)

    # dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup/net0/"
    # plot_compare_branches_cdf(dir_path, False)

    dir_path2 = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/"
    dir_path1 = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup/net0/"

    plot_out_branches_compare_cdf(dir_path1, dir_path2, False)