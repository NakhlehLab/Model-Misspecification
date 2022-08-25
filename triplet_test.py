import dendropy
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import scipy.stats as stats
from scipy.special import rel_entr
from itertools import combinations
import numpy as np

# from eval_3taxa import compute_KL
from calculator import *

species = []
taxon_map="<C:C_0;G:G_0;L:L_0;R:R_0;Q:Q_0;Z:Z_0>"
# SPECIES = ["C_0", "G_0", "L_0", "R_0", "Q_0"]
# OUTGROUP="Z_0"
SPECIES = ["C 0", "G 0", "L 0", "R 0", "Q 0"]
OUTGROUP="Z 0"
SUBSET_SIZE = 3
DEBUG = False
if DEBUG:
    NUM_REPLICATE = 3
else:
    NUM_REPLICATE = 100
NUM_LOCI = 100
GT_TOPO = 3

# global NUM_TRIPLET

gt_name = {True: "truegt", False:"infergt"}


def divide_triplet():
    return list(combinations(SPECIES, SUBSET_SIZE))
    
def compute_obs_frequency(tree_list, triplet):
    print(triplet)
    remove_taxa_list = [item for item in SPECIES if item not in triplet]
    remove_taxa_list.append(OUTGROUP)
    # print(remove_taxa_list)
    pruned_tree_list = []
    for tree in tree_list:
        tree.prune_taxa_with_labels(remove_taxa_list)
        # print(tree.as_string(schema="newick"))
        pruned_tree_list.append(tree)

    pruned_tree_list = dendropy.TreeList(pruned_tree_list)
    unique_topologies = pruned_tree_list.as_tree_array().topologies(sort_descending = True, frequency_attr_name='frequency')
    freqs = []
    major_gt = unique_topologies[0].as_string(schema="newick").strip()
    for tree in unique_topologies:
        freqs.append(tree.frequency)
    print(freqs)
    while len(freqs) != GT_TOPO:
        freqs.append(0.0)
    return freqs, major_gt


def get_obs_exp(tree_list, triplet, ):
    obs, major_gt = compute_obs_frequency(tree_list, triplet)
    obs = [x*NUM_LOCI for x in obs]
    t = internal_time(obs[0]/sum(obs))
    expct = compute_exp(t)
    exps = [x*sum(obs) for x in expct]
    # chi2, pvalue = chisquare(obs, exps, ddof=1)
    # print(f"pvalue={pvalue}")
    print(f"obs={obs}, exps={exps}")
    return obs, exps, t, major_gt

def compute_pvalue_species_tree(tree_path, triplets):
    
    tree_list = dendropy.TreeList.get(path=tree_path, schema="newick", rooting="force-rooted")
    res_dict = {"triplet": [], "time":[], "major_gt":[], "pvalue":[], "chi2":[]}

    for triplet in triplets:
        treelist_copy = tree_list.clone(depth=1)
        obs, exps, t, major_gt = get_obs_exp(treelist_copy, triplet)
        # print(obs, exps, t)
        chi2, pvalue = compute_Pvalue(obs, exps, SUBSET_SIZE - 2)
        res_dict["triplet"].append(triplet)
        res_dict["time"].append(t)
        res_dict["pvalue"].append(pvalue)
        res_dict["chi2"].append(chi2)
        res_dict["major_gt"].append("\""+major_gt+"\"")
        print(pvalue)
    print(res_dict)
    res_df = pd.DataFrame.from_dict(res_dict)
    print(res_df)
    return res_df

def compute_all_replicate(obs_dir, true_gt):
    df = None
    triplets = divide_triplet()
    NUM_TRIPLET = len(triplets)

    for i in range(1, NUM_REPLICATE + 1):
        print(i)
        if true_gt:
            path = obs_dir + "/" + str(i) + "/heter/genetrees.txt"
        else:
            path = obs_dir + "/" + str(i) + "/heter/rooted_iqtree.txt"
        df_i = compute_pvalue_species_tree(path, triplets)
        df_name = pd.DataFrame({"ID":[i for j in range(NUM_TRIPLET)]})
        # print(df_name)
        # print([i for i in range(NUM_TRIPLET)])
        df_i = df_i.join(df_name["ID"])
        print(df_i)
        if df is None:
            df = df_i
        else:
            df = pd.concat([df, df_i], ignore_index=True)
    print(df)
    df.to_csv(obs_dir+"/triplet_test_"+gt_name[true_gt]+".csv")
    

def plot_figures(obs_dir, true_gt, plot):
    csv_path = obs_dir+"/triplet_test_"+gt_name[true_gt]+".csv"
    df = pd.read_csv(csv_path)

    plt.clf()
    if plot == "id-p":
        fig = plt.figure(figsize=(18, 6), dpi=80)
        fig = sns.scatterplot(data=df, x="ID", y="pvalue", hue="triplet",style="triplet", s=20)
        fig.axhline(0.05)
        plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
        if not DEBUG:
            plt.xticks(range(1, NUM_REPLICATE+1, 5), fontsize=14)
        else:
            plt.xticks(range(1, NUM_REPLICATE+1), fontsize=14)
        plt.yticks(fontsize=14)
        plt.ylabel("P-value", fontsize=16)
        plt.xlabel("Replicate", fontsize=16)

        
        plt.tight_layout()
        plt.savefig(obs_dir+"/triplet_pvalue_"+gt_name[true_gt]+".pdf")
        plt.show()
    elif plot == "p-hist":
        fig = sns.histplot(data=df, x="pvalue", stat="probability", bins=20)
        
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        plt.xlabel("P-value", fontsize=16)
        plt.ylabel("Proportion", fontsize=16)
        plt.ylim([0, 1])
        plt.xlim(0, )
        
        plt.tight_layout()
        plt.savefig(obs_dir+"/triplet_p_hist_"+gt_name[true_gt]+".pdf")
        plt.show()
    
    elif plot == "chi2-hist":
        fig = sns.histplot(data=df, x="chi2", stat="probability")
        x = np.arange(0, 10, .01)
        plt.plot(x, stats.chi2.pdf(x, df=SUBSET_SIZE-2), color='r', lw=2)
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        plt.xlabel("chi2", fontsize=16)
        plt.ylabel("Proportion", fontsize=16)
        plt.ylim([0, 1])
        plt.xlim(0, )

        
        plt.tight_layout()
        plt.savefig(obs_dir+"/triplet_chi2_hist_"+gt_name[true_gt]+".pdf")
        plt.show()




if __name__ == "__main__":
    # tree_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_tall10/036/2/heter/genetrees.txt"
    # tree_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_tall10/036/1/heter/rooted_iqtree.txt"
    # triplets = divide_triplet()
    # compute_pvalue_species_tree(tree_path, triplets)
    
    
    # obs_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_tall10/036/"
    # # true_gt = True
    # true_gt = False

    # # plot = "id-p"
    # # plot = "p-hist"
    # plot = "chi2-hist"
    # # compute_all_replicate(obs_dir, true_gt)
    # plot_figures(obs_dir, true_gt, plot)

    all_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_tall10/"
    for pop_size in ["01", "001"]:
        obs_dir = all_dir + pop_size + "/"
        for true_gt in [True, False]:
            compute_all_replicate(obs_dir, true_gt)

    for pop_size in ["036", "01", "001"]:
        obs_dir = all_dir + pop_size + "/"
        for true_gt in [True, False]:
            for plot in ["id-p", "p-hist","chi2-hist"]:
                plot_figures(obs_dir, true_gt, plot)
    