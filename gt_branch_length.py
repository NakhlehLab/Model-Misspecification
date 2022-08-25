import pandas as pd
import dendropy
from regex import D, F
import seaborn as sns
import matplotlib.pyplot as plt
import sys

SPECIES = ["C 0", "G 0", "L 0", "Q 0", "R 0", "Z 0"]


def compute_simulated_gt_pair_dist(tree_path, dist_csv, theta):
    tree_list = dendropy.TreeList.get(path=tree_path, schema="newick")
    dict_dist = {"species_pair":[], "distances":[]}

    for tree in tree_list:
        pdm = tree.phylogenetic_distance_matrix()
        for idx1 in range(len(SPECIES) - 1):
            idx2 = idx1 + 1
            taxon1 = tree.taxon_namespace.get_taxon(SPECIES[idx1])
            while idx2 < len(tree.taxon_namespace):
                taxon2 = tree.taxon_namespace.get_taxon(SPECIES[idx2])
                weighted_patristic_distance = pdm.patristic_distance(taxon1, taxon2)
                dict_dist["species_pair"].append(taxon1.label+"-"+taxon2.label)
                dict_dist["distances"].append(weighted_patristic_distance*theta)
                idx2 += 1
    res_df = pd.DataFrame(dict_dist)  
    res_df.to_csv(dist_csv)
    return res_df


def compute_simulated_gt_brs(tree_path, murate_path, ingroup, outgroup = "Z 0"):
    tree_list = dendropy.TreeList.get(path=tree_path, schema="newick")
    murate_list = []
    if murate_path != None:
        with open(murate_path, "r") as handle:
            for line in handle.readlines():
                murate_list.append(float(line.strip()))
        
    dist = []
    if len(murate_list) == 0:
        for tree in tree_list:
            if ingroup:
                tree.prune_taxa_with_labels([outgroup])
                tree.is_rooted = False
            tree.encode_bipartitions()
            for node in tree.preorder_node_iter():
                if node.parent_node is not None:
                    dist.append(node.edge.length)
    else:
        for lcs, tree in enumerate(tree_list):
            if ingroup:
                tree.prune_taxa_with_labels([outgroup])
                tree.is_rooted = False
            tree.encode_bipartitions()
            for node in tree.preorder_node_iter():
                if node.parent_node is not None:
                    dist.append(node.edge.length * murate_list[lcs])
    return dist


def run_all_branches(dir_path, num_replicate, THETA_DICT, ingroup = True, longer=True, simulated_gt = True, log=True):
    df_list = []
    if longer:
        locus_length_list = [500, 2000]
    else:
        locus_length_list = [200, 1000]
    

    for scale in ["short", "medium", "long"]:
        for locus_length in locus_length_list:

            for i in range(1, num_replicate + 1):
                if simulated_gt:
                    path = dir_path +scale+"/"+ str(i) + "/homo/"+str(locus_length) + "/genetrees.txt"
                    murate_path = dir_path +scale+"/"+ str(i) + "/heter/"+str(locus_length) + "/murate.txt"
                else:
                    path = dir_path +scale+"/"+ str(i) + "/heter/"+str(locus_length) + "/iqtree.txt"
                    murate_path = None
                
                dist = compute_simulated_gt_brs(path, murate_path, ingroup)
                df = pd.DataFrame()
                if simulated_gt:
                    df["branch_lengths"] = [x*THETA_DICT[scale] for x in dist]
                else:
                    df["branch_lengths"] = dist
                df["scale"] = [scale for x in range(len(df["branch_lengths"]))]
                df["locus_length"] = [locus_length for x in range(len(df["branch_lengths"]))]
                df_list.append(df)
    res_df = pd.concat(df_list)
    if simulated_gt:
        res_df.to_csv(dir_path+"/gt_branches_distribution_truegt.csv")
    else:
        res_df.to_csv(dir_path+"/gt_branches_distribution_iqgt.csv")
    print(res_df)
    plt.clf()
    ax = sns.violinplot(x="scale", y="branch_lengths", data=res_df, hue="locus_length")
    log_dict={True: "_log", False:""}
    if log:
        ax.set_yscale('log')
        ax.set_ylim(sys.float_info.min, )
    plt.grid(axis = 'y')

    if ingroup:
        if simulated_gt:
            plt.savefig(dir_path+"/gt_branches_ingroup_truegt_distribution"+log_dict[log]+".pdf")
        else:
            plt.savefig(dir_path+"/gt_branches_ingroup_iqgt_distribution"+log_dict[log]+".pdf")
    else:
        if simulated_gt:
            plt.savefig(dir_path+"/gt_branches_truegt_distribution"+log_dict[log]+".pdf")
        else:
            plt.savefig(dir_path+"/gt_branches_iqgt_distribution"+log_dict[log]+".pdf")


def run_all_pair(dir_path, num_replicate, THETA_DICT):
    df_list = []
    for scale in ["short", "medium", "long"]:
        for locus_length in [200, 1000]:
            for i in range(1, num_replicate+1):
                sub_dir = dir_path +scale+"/"+ str(i) + "/homo/"+str(locus_length)
                path = sub_dir + "/genetrees.txt"
                csv_path = sub_dir + "gt_dist.csv"
                df = compute_simulated_gt_pair_dist(path, csv_path, THETA_DICT[scale])
                df["scale"] = [scale for x in range(len(df["distances"]))]
                df["locus_length"] = [locus_length for x in range(len(df["distances"]))]
                df_list.append(df)
    res_df = pd.concat(df_list)
    res_df.to_csv(dir_path+"/gt_dist_.csv")
    fig, axes = plt.subplots(1, 2, figsize=(20,5))
    for i, locus_length in enumerate([200, 1000]):
        sns.boxplot(x="species_pair", y="distances", data=res_df[res_df["locus_length"] == locus_length], hue="scale", ax=axes[i], showfliers = False)
        axes[i].set_xticklabels(axes[i].get_xticklabels(),rotation=30)
    plt.tight_layout()
    plt.savefig(dir_path+"/gt_pair_distances.pdf")

def plot_pair_dist(dir_path):
    res_df = pd.read_csv(dir_path+"/gt_dist_.csv")
    plt.clf()
    fig, axes = plt.subplots(1, 2, figsize=(20,5))
    for i, locus_length in enumerate([200, 1000]):
        sns.boxplot(x="species_pair", y="distances", data=res_df[res_df["locus_length"] == locus_length], hue="scale", ax=axes[i], showfliers = False)
        axes[i].set_xticklabels(axes[i].get_xticklabels(),rotation=30)
    plt.tight_layout()
    plt.savefig(dir_path+"/gt_pair_distances.pdf")


def plot_branches(dir_path):
    plt.clf()
    fig = plt.figure(figsize=(10,8))
    res_df = pd.read_csv(dir_path+"/gt_branches_distribution.csv")
    fig = sns.violinplot(x="scale", y="branch_lengths", data=res_df, hue="locus_length")
    plt.savefig(dir_path+"/gt_branches_distribution.pdf")

def plot_compare_branches(dir_path, longer):
    plt.clf()
    fig, axes = plt.subplots(2, 2, figsize=(20,10))
    if longer:
        locus_length_list = [500, 2000]
    else:
        locus_length_list = [200, 1000]

    truegt_df = pd.read_csv(dir_path+"/gt_branches_distribution_truegt.csv")
    iqgt_df = pd.read_csv(dir_path+"/gt_branches_distribution_iqgt.csv")
    
    for idx, locus_length in enumerate(locus_length_list):
        for col, df in enumerate([truegt_df, iqgt_df]):
            print(df[df["locus_length"] == locus_length])
            sns.histplot(df[df["locus_length"] == locus_length], x="branch_lengths", hue="scale", multiple="dodge", hue_order=["short", "medium", "long"], stat="probability", fill=True, binwidth=0.005, ax=axes[idx][col])
            
            axes[idx][col].set_ylim(0, 0.2)
            axes[idx][col].set_xlim(0, 0.4)
            # axes[idx][col].set_xscale('log')
            axes[idx][col].set_xlabel("Branch lengths in substitution rates", fontsize=18)
            axes[idx][col].set_ylabel("Proportion of branches", fontsize=18)
            
    axes[0][0].set_title("true gt", size=20, fontweight="bold")
    axes[0][1].set_title("inferred gt", size=20, fontweight="bold")
    axes[0][0].set_ylabel(f"Proportion (locus length={locus_length_list[0]})", fontsize=20)
    axes[1][0].set_ylabel(f"Proportion (locus length={locus_length_list[1]})", fontsize=20)
    labels = axes[0][0].get_legend()
    handles = labels.legendHandles
    # handles, labels = axes[0][0].get_legend_handles_labels()
    print(handles, labels)
    # fig.legend(handles, labels, loc='lower center', fontsize=16, ncol=3, )

    plt.savefig(dir_path+"/gt_branches_distribution_compare.pdf")

def plot_compare_branches_cdf(dir_path, longer):
    plt.clf()
    fig, axes = plt.subplots(1, 2, figsize=(15,5))
    if longer:
        locus_length_list = [500, 2000]
    else:
        locus_length_list = [200, 1000]

    truegt_df = pd.read_csv(dir_path+"/gt_branches_distribution_truegt.csv", index_col=False)
    iqgt_df = pd.read_csv(dir_path+"/gt_branches_distribution_iqgt.csv", index_col=False)

    truegt_df["truegt"] = [True for x in range(len(truegt_df["branch_lengths"]))]
    iqgt_df["truegt"] = [False for x in range(len(iqgt_df["branch_lengths"]))]
    
    df = pd.concat([truegt_df, iqgt_df], ignore_index=True)
    # df["scale, truegt"] = df[["scale","truegt"]].apply(tuple, axis=1)
    # styles = {True: "solid", False: "dotted"}

    for idx, locus_length in enumerate(locus_length_list):
        print(df[df["locus_length"] == locus_length])
        hue = df[['scale', 'truegt']].apply(
        lambda row: f"{row.scale}, {row.truegt}", axis=1)
        hue.name = 'scale, truegt'
        sns.ecdfplot(df[df["locus_length"] == locus_length], x="branch_lengths", hue=hue, stat="proportion", ax=axes[idx], linewidth=2, alpha=0.6)
        
        lss = ['--','--','--','-','-', '-']
        colors = ["blue", "orange", "green", "blue", "orange", "green"]
        # labels = axes[idx].get_legend()
        # handles = labels.legendHandles
        handles = axes[idx].legend_.legendHandles[::-1]

        for line, ls, handle, clr in zip(axes[idx].lines, lss, handles, colors):
            line.set_linestyle(ls)
            line.set_color(clr)
            handle.set_ls(ls)
            handle.set_color(clr)

        axes[idx].tick_params(axis='y', labelsize=14)
        axes[idx].tick_params(axis='x', labelsize=14)
        axes[idx].set_ylim(0, 1.03)
        axes[idx].set_xscale('log')
        axes[idx].set_xlim(0.0001, 0.3)
        axes[idx].set_xlabel("Branch lengths in substitution rates", fontsize=18)
        axes[idx].set_ylabel("Proportion of branches", fontsize=18)

    axes[0].set_title("locus length="+str(locus_length_list[0]), size=20, fontweight="bold")
    axes[1].set_title("locus length="+str(locus_length_list[1]), size=20, fontweight="bold")
    
    # labels = axes[0][0].get_legend()
    # handles = labels.legendHandles
    # print(handles, labels)
    # fig.legend(handles, new_labels, loc='lower center', fontsize=20, ncol=3)
    plt.savefig(dir_path+"/gt_branches_distribution_compare_cdf.pdf")



if __name__ == "__main__":
    # tree_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/long/1/heter/500/iqtree.txt"
    
    # dist = compute_simulated_gt_brs(tree_path, True)
    # print(len(dist))
    
    
    # dir_path = sys.argv[1]
    # num_replicate = int(sys.argv[2])
    # ingroup = (sys.argv[3].lower() == "true")
    # longer = (sys.argv[4].lower() == "true")
    # truegt = (sys.argv[5].lower() == "true")
    # log = (sys.argv[6].lower() == "true")
    # larger_theta = (sys.argv[7].lower() == "true")
    # # run_all_pair(dir_path, num_replicate)

    # if larger_theta == False:
    #     THETA_DICT = {"short":0.025, "medium":0.0125, "long":0.005}
    # else:
    #     THETA_DICT = {"short":0.05, "medium":0.025, "long":0.01}

    # run_all_branches(dir_path, num_replicate, THETA_DICT, ingroup, longer, truegt, log)

    # plot_pair_dist(dir_path)
    # plot_branches(dir_path)


    dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup/net0/"
    plot_compare_branches_cdf(dir_path, False)

    # dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/"
    # plot_compare_branches_cdf(dir_path, True)
    
