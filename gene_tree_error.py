
import scipy.stats
import dendropy
from dendropy.calculate import treecompare
import sys
import os
import matplotlib.pyplot as plt
import seaborn as sns

SCALE_TRUE_GT = 0.018


def tree_ingroup_euclidean_distance(truetree: dendropy.Tree, tree2: dendropy.Tree, OUTGROUP=None):
    truetree.encode_bipartitions()
    tree2.encode_bipartitions()

    if OUTGROUP is not None:
        truetree.prune_taxa_with_labels([OUTGROUP])
        tree2.prune_taxa_with_labels([OUTGROUP])

    distance_ingroup = treecompare.euclidean_distance(truetree, tree2)/truetree.length()
    # print(distance_ingroup)
    return distance_ingroup

def tree_scale(tree, scale):
    for idx, nd in enumerate(tree):
        if nd.edge.length is not None:
            nd.edge.length *= scale
    return tree

def compare_iqtree_truetree(truetree_path, iqtree_path, outgroup):
    tns = dendropy.TaxonNamespace()
    iqtree_list = dendropy.TreeList.get(path=iqtree_path, schema="newick", rooting="force-rooted", taxon_namespace=tns)
    truetree_list = dendropy.TreeList.get(path=truetree_path, schema="newick", rooting="force-rooted",
                                          taxon_namespace=tns)

    true_cnt = 0
    # rf_distance = 0
    nrf_list = []
    nrBS_list = []
    # wrf_distance = 0
    # nrBS = 0
    num_taxa = len(iqtree_list[0].taxon_namespace)
    for iqtree, truetree in zip(iqtree_list, truetree_list):
        truetree = tree_scale(truetree, SCALE_TRUE_GT)
        rf = treecompare.symmetric_difference(iqtree, truetree)
        # print("-------------------------")
        # print(iqtree.as_string("newick"))
        # print(truetree.as_string("newick"))
        # print(rf/(2*(num_taxa-2)))

        # rf_distance += rf/(2*(num_taxa-2))
        nrf_list.append(rf/(2*(num_taxa-2)))
        nrBS_list.append(tree_ingroup_euclidean_distance(truetree, iqtree, outgroup))
        if rf < 0.01:
            true_cnt += 1
        # nrBS += tree_ingroup_euclidean_distance(truetree, iqtree, outgroup)
    # rf_distance /= len(iqtree_list)
    # nrBS /= len(iqtree_list)
    # print("Number of correctly inferred gene trees: ", true_cnt)
    # print("Average normalized RF distance: ", rf_distance)
    # print("Average weighted RF distance: ", wrf_distance)
    # print("Average nrBS of ingroup: ", nrBS)
    return nrf_list, nrBS_list

def plot_RF_hist(rf_list, fig_path):
    sns.histplot(x=rf_list, bins=30, stat="probability")
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylabel("Proportion", fontsize=16)
    plt.xlabel("nrRF", fontsize=16)
    plt.tight_layout()
    plt.savefig(fig_path)
    # plt.show()

def plot_all_replicate(directory, rate_dir=None):
    heter = 1
    
    # directory = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_100/036/"
    # directory = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa3_tall/036/"
    # directory = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa3_100/036/"
    # directory = "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/data/simulation/seqgen/replica/tree/"
    # truegt_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_100/036/1/heter/genetrees.txt"
    # iqtree_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_100/036/1/heter/rooted_iqtree.txt"
    if heter == 1:
        if rate_dir is None:
            fig_path = directory + "gt_rf_heter.pdf"
        else:
            fig_path = directory + "gt_rf_heter_"+rate_dir+".pdf"
    else:
        fig_path = directory + "gt_rf_heter_locus"+".pdf"
    rf_list = []
    outgroup = "Z 0"
    # outgroup = "G 0"
    for i in range(1, 101):
        if heter == 1:
            truegt_path = directory + str(i) + "/heter/genetrees.txt"
            if rate_dir is None:
                
                iqtree_path = directory + str(i) + "/heter/rooted_iqtree.txt"
            else:
                iqtree_path = directory + str(i) + "/heter/"+rate_dir+"/rooted_iqtree.txt"
        else:
            truegt_path = directory + str(i) + "/heter_locus/genetrees.txt"
            iqtree_path = directory + str(i) + "/heter_locus/rooted_iqtree.txt"
        # truegt_path = directory + str(i) + "/heter_locus/genetrees.txt"
        # iqtree_path = directory + str(i) + "/heter_locus/rooted_iqtree.txt"

        rf, nrbs = compare_iqtree_truetree(truegt_path, iqtree_path, outgroup)
        rf_list.extend(rf)
    plot_RF_hist(rf_list, fig_path)

if __name__ == "__main__":
    plot_all_replicate(sys.argv[1], sys.argv[2])