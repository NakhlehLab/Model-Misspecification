import os
from random import seed
import matplotlib.pyplot as plt
from numpy import average
from scipy.stats import bayes_mvs
import dendropy
import sys
import numpy as np
import pandas as pd


# sys.path.append('../')
# from gt_error_profling import tree_ingroup_euclidean_distance
import time
from dendropy.calculate import treecompare

CHAIN_LEN = 20000000
BURN_IN = 2000000
SAMPLE_FREQ = 5000
SEED_INDEX=0
NUM_SEED = 2
NUM_RUN = 8
NUM_REPLICATE = 10
NUM_LOCI = 100
SPECIES = ["C 0", "L 0", "R 0", "Q 0"]

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

# def scale_gts(murate_path, gt_path):
#     murates = []
#     with open(murate_path, "r") as handle:
#         for line in handle.readlines():
#             murates.append(float(line.strip()))
#     treelist = dendropy.TreeList.get(
#         path=gt_path,
#         schema="newick",
#         rooting="force-rooted"
#     )
#     i = 0
#     for tree in treelist:
#         treelist[i] = tree_scale(tree, murates[i])
#         i += 1
#     i = 0

#     for tree in treelist:
#         print("Tree geneTree"+str(i)+"="+tree.as_string(schema="newick")[5:])
#         i += 1
#     return treelist

def read_mcmc_mu_rates(dir_path):
    locus2rates = {}
    for i in range(1, NUM_LOCI+1):
        locus2rates["loci"+str(i)] = []

    for i in range(NUM_RUN):
        sub_dir = dir_path + "/mcmcseq/0_"+str(SEED_INDEX) + "/"+str(CHAIN_LEN)+"_"+str(i)+"/"
        if not os.path.exists(sub_dir):
            break
    
        for file_name in os.listdir(sub_dir):
            if file_name.startswith("mutationRate_"):
                name = file_name.split()[0][13:].split(".")[0]
                # locus2rates[name] = []
                with open(os.path.join(sub_dir, file_name), "r") as handle:
                    for line in handle.readlines():
                        if line.startswith("--"):
                            break
                        locus2rates[name].append(float(line.strip()))
                    handle.close()
    return locus2rates

def test_murate_sum(locus2rates):
    num_iters = 0
    for key in locus2rates.keys():
        num_iters = len(locus2rates[key])
        break
    for i in range(num_iters):
        sum = 0
        for key in locus2rates.keys():
            sum += locus2rates[key]
        print(sum)

def read_true_mu_path(true_mu_path):
    with open(true_mu_path, "r") as handle:
        lines = handle.read().strip().split("\n")
        true_mu_rates = [float(x.strip()) for x in lines]
    return true_mu_rates

def compare_true_mcmc_inferred_murates(true_mu_path, mcmc_dir, burnin=10):
    true_mu_rates = read_true_mu_path(true_mu_path)
    locus2rates = read_mcmc_mu_rates(mcmc_dir)
    inferred_mu_rates = []
    start = int(burnin / 100 * len(locus2rates["loci1"]))
    for i in range(1, len(true_mu_rates)+1):
        inferred_mu_rates.append(sum(locus2rates["loci"+str(i)][start:])/len(locus2rates["loci"+str(i)][start:]))
    fig, ax = plt.subplots()
    ax.scatter(true_mu_rates, inferred_mu_rates, c="black", marker=".")
    plt.xlabel("true mutation rates")
    plt.ylabel("mcmc inferred mutation rates")
    plt.plot(true_mu_rates, true_mu_rates, '--', color="blue", )
    plt.xlim(0.0, )
    plt.ylim(0.0, )
    plt.savefig(os.path.join(mcmc_dir,"mcmc_murates.pdf"))
    # plt.show()

def trace_mcmc_inferred_murates_each_locus(mcmc_dir, burnin=10):
    locus2rates = read_mcmc_mu_rates(mcmc_dir)
    inferred_mu_rates = []

    i = 0
    k = "loci1"
    length = len(locus2rates[k])
    for k, v in locus2rates.items():
        if length > len(locus2rates[k]):
            length = len(locus2rates[k])

    for j in range(length):
        sum = 0
        for k, v in locus2rates.items():
            sum += locus2rates[k][j]
            # inferred_mu_rates.append(sum(locus2rates["locus"+str(i)])/len(locus2rates["locus"+str(i)]))
            i += 1
        print(sum)
    # plt.ylim(0.0, )

    start = int(burnin * length/100)
    for k, v in locus2rates.items():
        fig, ax = plt.subplots()
        ax.plot(range(len(locus2rates[k][start:])), locus2rates[k][start:], c="black")
        plt.savefig(os.path.join(mcmc_dir, "trace_murates_" + k + ".pdf"))
    # plt.show()


# def read_beast_mu_rates(path):
#     with open(path, "r") as handle:
#         begin = False
#         index2key = {}
#         locus2rates = {}
#         for line in handle.readlines():
#             if line.startswith("Sample"):
#                 begin = True
#                 arr = line.strip().split("\t")
#                 index = 0
#                 for x in arr:
#                     if x.startswith("mutationRate.s"):
#                         locus_name = x.split(":")[1]
#                         index2key[index] = locus_name
#                         locus2rates[locus_name] = []
#                     index += 1
#             elif begin:
#                 values = line.strip().split("\t")
#                 for i in index2key.keys():
#                     locus2rates[index2key[i]].append(float(values[i]))
#     return locus2rates


# def compare_true_beast_inferred_murates(true_mu_path, beast_log_path):
#     with open(true_mu_path, "r") as handle:
#         lines = handle.read().strip().split("\n")
#         true_mu_rates = [float(x.strip()) for x in lines]
#     locus2rates = read_beast_mu_rates(beast_log_path)
#     inferred_mu_rates = []
#     lows = []
#     highs = []
#     for i in range(len(true_mu_rates)):
#         mean, var, std = bayes_mvs(locus2rates["locus" + str(i)], alpha=0.95)
#         inferred_mu_rates.append(mean[0])
#         lows.append(mean[1][0])
#         highs.append(mean[1][1])
#         # inferred_mu_rates.append(sum(locus2rates["locus" + str(i)]) / len(locus2rates["locus" + str(i)]))
#     print(true_mu_rates)
#     print(inferred_mu_rates)

#     fig, ax = plt.subplots()
#     ax.scatter(true_mu_rates, inferred_mu_rates, c="black", marker=".")
#     plt.xlabel("true mutation rates")
#     plt.ylabel("beast inferred mutation rates")
#     plt.plot(true_mu_rates, true_mu_rates, '--', color="blue", )
#     plt.xlim(0.0, 4.0)
#     plt.ylim(0.0, 4.0)
#     plt.savefig("beast.pdf")


def read_mcmc_gene_trees(dir_path):
    locus2gt = {}
    for i in range(1, NUM_LOCI+1):
        locus2gt["loci"+str(i)] = []

    for i in range(NUM_RUN):
        sub_dir = dir_path + "/mcmcseq/0_"+str(SEED_INDEX) + "/"+str(CHAIN_LEN)+"_"+str(i)+"/"
        if not os.path.exists(sub_dir):
            break
        for file_name in os.listdir(sub_dir):
            if file_name.startswith("tree_"):
                name = file_name.split()[0][5:].split(".")[0]
                # locus2gt[name] = []
                # print(name)
                with open(os.path.join(sub_dir, file_name)) as handle:
                    for line in handle.readlines():

                        if line.startswith("--"):
                            break
                        locus2gt[name].append(line.strip())
    return locus2gt


def read_mcmc_popsize(dir_path):
    popsize_list = []
    for i in range(NUM_RUN):
        sub_dir = dir_path + "/mcmcseq/0_"+str(SEED_INDEX) + "/"+str(CHAIN_LEN)+"_"+str(i)+"/"
        if not os.path.exists(sub_dir):
            break
        for filename in os.listdir(sub_dir):
            if filename.startswith("slurm"):
                with open(os.path.join(sub_dir, filename)) as handle:
                    begin = False
                    for line in handle.readlines():
                        if begin:
                            if line.startswith("["):
                                popsize_list.append(float(line[1:line.index("]")]))
                        elif line.startswith("Iteration"):
                            begin = True
    return popsize_list


"""
plot trace of gene tree errors
Todo: test
"""
def trace_mcmc_gene_tree_errors(dir_path, truetree_path, bl, sf, outgroup=None, burnin=0, plot=False):
    starttime = time.time()
    locus2gt = read_mcmc_gene_trees(dir_path)
    popsize_list = read_mcmc_popsize(dir_path)
    tns = dendropy.TaxonNamespace()
    truetree_list = dendropy.TreeList.get(path=truetree_path, schema="newick", rooting="force-rooted",
                                          taxon_namespace=tns)

    print("Finished reading...")
    print(len(locus2gt["loci1"]), len(popsize_list))
    length = len(locus2gt["loci1"])

    for k in locus2gt.keys():
        if length > len(locus2gt[k]):
            length = len(locus2gt[k])

    print(length)
    start = max(int(burnin * length // 100), int(bl//sf + 1))
    theta_mean = np.mean(popsize_list[start:])

    taxon_list = []
    for taxon in truetree_list[0].taxon_namespace:
        taxon_list.append(taxon.label)
    taxon_to_remove = list(set(taxon_list)-set(SPECIES))
    for j in range(len(truetree_list)):
        truetree_list[j].prune_taxa_with_labels(taxon_to_remove)
        tree_scale(truetree_list[j], theta_mean / 2)

    rf_list = []
    nrbs_list = []

    for i in range(start, int(length/100)*100):
        rf_distance = 0
        nrbs = 0

        for k in locus2gt.keys():
            try:
                inferred_tree = dendropy.Tree.get(data=locus2gt[k][i], schema="newick", rooting="force-rooted",
                                          taxon_namespace=tns)
            except:
                print(f"An exception occurred:{k, i, locus2gt[k][i]}")
            # true_tree = truetree_list[int(k[5:])]
            true_tree = truetree_list[int(k[4:])-1]

            rf_distance += treecompare.symmetric_difference(inferred_tree, true_tree)
            nrbs += tree_ingroup_euclidean_distance(true_tree, inferred_tree, outgroup)

        rf_distance /= (2 * (len(tns) - 2))
        rf_distance /= len(locus2gt)
        nrbs /= len(locus2gt)
        if i % 100 == 0:
            print(f"iteration={i}, rf_distance={rf_distance}, nrbs={nrbs}")
        rf_list.append(rf_distance)
        nrbs_list.append(nrbs)
    if plot:
        fig, ax = plt.subplots()
        ax.plot(range(length), rf_list, c="green")
        plt.savefig(os.path.join(dir_path, "trace_gt_RF.pdf"))

        fig, ax = plt.subplots()
        ax.plot(range(length), nrbs_list, c="green")
        plt.savefig(os.path.join(dir_path, "trace_gt_nrBS.pdf"))
    with open(os.path.join(dir_path, "mcmc_gt_err.txt"), "w") as handle:
        handle.write(f"average RF distance: {sum(rf_list)/len(rf_list)}\n")
        handle.write(f"average nrBS: {sum(nrbs_list) / len(nrbs_list)}\n")
        handle.write(f"Run Time:{time.time()-starttime}\n")
        handle.close()
    # print("Run Time:", time.time()-starttime)

def make_canonical(node):
	min_labels = []
	node_children = node.child_nodes()
	for child in node_children:
		if child.is_leaf():
			min_labels.append(child.taxon.label)
		else:
			min_labels.append(make_canonical(child))
	if min_labels[0] > min_labels[1]:
		node.set_child_nodes(reversed(node_children))
	return min(min_labels)

def sum_tree_list(treelist):
    unique_topologies = treelist.as_tree_array().topologies(sort_descending = True, frequency_attr_name='frequency')
    most_freq_trees = []
    unique_topologies[0].prune_taxa_with_labels(["Z 0"])
    res_tree = None
    # mcct = treelist.maximum_product_of_split_support_tree()
    # con_tree = treelist.consensus(min_freq=0.95)
    brl_list_internal = []
    # brl_list_leaf = []
    for tree in treelist:
        if treecompare.symmetric_difference(unique_topologies[0], tree) < 0.01:
            make_canonical(tree.seed_node)
            if res_tree is None:
                res_tree = tree.clone()
            # print(tree.as_string(schema="newick"))
            dist = []
            for nd in tree.preorder_node_iter():
                if nd.edge_length is not None:
                    dist.append(nd.edge_length)
            most_freq_trees.append(tree)
            brl_list_internal.append(dist)
            # dist = []
            # for nd in tree.leaf_node_iter():
            #     dist.append(nd.edge_length)
            # brl_list_leaf.append(dist)
            # pdc = tree.phylogenetic_distance_matrix()
    # dist_list = np.array(dist_list)
    # x = dist_list.sum(axis=0)
    # res_tree =  most_freq_trees[0]
    i = 0
    internal_brls = np.array(brl_list_internal).sum(axis=0)/len(brl_list_internal)
    for nd in res_tree.preorder_node_iter():
        if nd.edge_length is not None:
            nd.edge_length = internal_brls[i]
            i += 1
    
    # leaf_brls = np.array(brl_list_leaf).sum(axis=0)/len(brl_list_leaf)
    return res_tree
        

def trace_mcmc_gene_tree_errors_2(dir_path, truetree_path, bl, sf, outgroup=None, burnin=0, plot=False):
    starttime = time.time()
    locus2gt = read_mcmc_gene_trees(dir_path)
    popsize_list = read_mcmc_popsize(dir_path)
    tns = dendropy.TaxonNamespace()
    truetree_list = dendropy.TreeList.get(path=truetree_path, schema="newick", rooting="force-rooted",
                                          taxon_namespace=tns)

    print("Finished reading...")
    print(len(locus2gt["loci1"]), len(popsize_list))
    length = len(locus2gt["loci1"])

    for k in locus2gt.keys():
        if length > len(locus2gt[k]):
            length = len(locus2gt[k])

    end = length//100 * 100

    for k in locus2gt.keys():
        locus2gt[k] = locus2gt[k][:end]

    print(length)
    start = max(int(burnin * length // 100), int(bl//sf + 1))
    theta_mean = np.mean(popsize_list[start:])

    taxon_list = []
    for taxon in truetree_list[0].taxon_namespace:
        taxon_list.append(taxon.label)

    taxon_to_remove = list(set(taxon_list)-set(SPECIES))
    for j in range(len(truetree_list)):
        truetree_list[j].prune_taxa_with_labels(taxon_to_remove)
        tree_scale(truetree_list[j], theta_mean / 2)

    rf_list = []
    nrbs_list = []
    gt_df = pd.DataFrame(locus2gt)
    end = length//100 * 100
    num_taxa = len(truetree_list[0].taxon_namespace)
    
    for k in locus2gt.keys():
        print(k)
        # try:
        # print(list(gt_df[k].loc[start:end]))
        inferred_tree_list = dendropy.TreeList.get(data="\n".join(list(gt_df[k].loc[start:end])), schema="newick", rooting="force-rooted",
                                    taxon_namespace=tns)
        inferred_tree = sum_tree_list(inferred_tree_list)
        true_tree = truetree_list[int(k[4:])-1]
        rf_distance = treecompare.symmetric_difference(inferred_tree, true_tree)/(2*(num_taxa-2))
        nrbs = tree_ingroup_euclidean_distance(true_tree, inferred_tree, outgroup)
        rf_list.append(rf_distance)
        nrbs_list.append(nrbs)
        print(rf_distance, nrbs)
        # except:
        #     print(f"An exception occurred:{k}")

    # for i in range(start, int(length/100)*100):
    #     rf_distance = 0
    #     nrbs = 0

    #     for k in locus2gt.keys():
    #         try:
    #             inferred_tree = dendropy.Tree.get(data=locus2gt[k][i], schema="newick", rooting="force-rooted",
    #                                       taxon_namespace=tns)
    #         except:
    #             print(f"An exception occurred:{k, i, locus2gt[k][i]}")
    #         # true_tree = truetree_list[int(k[5:])]
    #         true_tree = truetree_list[int(k[4:])-1]

    #         rf_distance += treecompare.symmetric_difference(inferred_tree, true_tree)
    #         nrbs += tree_ingroup_euclidean_distance(true_tree, inferred_tree, outgroup)

    #     rf_distance /= (2 * (len(tns) - 2))
    #     rf_distance /= len(locus2gt)
    #     nrbs /= len(locus2gt)
    #     if i % 100 == 0:
    #         print(f"iteration={i}, rf_distance={rf_distance}, nrbs={nrbs}")
    #     rf_list.append(rf_distance)
    #     nrbs_list.append(nrbs)
    if plot:
        fig, ax = plt.subplots()
        ax.plot(range(length), rf_list, c="green")
        plt.savefig(os.path.join(dir_path, "trace_gt_RF.pdf"))

        fig, ax = plt.subplots()
        ax.plot(range(length), nrbs_list, c="green")
        plt.savefig(os.path.join(dir_path, "trace_gt_nrBS.pdf"))
    with open(os.path.join(dir_path, "mcmc_gt_err.txt"), "w") as handle:
        handle.write(f"average RF distance: {sum(rf_list)/len(rf_list)}\n")
        handle.write(f"average nrBS: {sum(nrbs_list) / len(nrbs_list)}\n")
        handle.write(f"Run Time:{time.time()-starttime}\n")
        handle.close()

def mcmc_murates_error(mcmc_dir, true_mu_path, burnin=0, plot=False):
    locus2rates = read_mcmc_mu_rates(mcmc_dir)
    true_mu_rates = read_true_mu_path(true_mu_path)

    # length = len(locus2rates["locus0"])
    length = len(locus2rates["loci1"])

    for k in locus2rates.keys():
        if length > len(locus2rates[k]):
            length = len(locus2rates[k])

    start = int(burnin * length / 100)

    error_list = []
    for i in range(start, length):
        err = 0
        for k in range(len(true_mu_rates)):
            err += abs(locus2rates["loci" + str(k+1)][i] - true_mu_rates[k])
        error_list.append(err / len(true_mu_rates))
    if plot:
        fig, ax = plt.subplots()
        ax.plot(range(len(error_list)), error_list, c="black")
        plt.savefig(os.path.join(mcmc_dir, "trace_murates_error.pdf"))
    with open(os.path.join(mcmc_dir, "mcmc_murate_err.txt"), "w") as handle:
        handle.write(f"average mutation rate error: {sum(error_list) / len(error_list)}\n")
        handle.close()

if __name__ == '__main__':
    # dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/net/mcmc/heter_locus/jc"
    # true_gt_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/data/simulation/seqgen/net/heter_locus/genetrees.txt"
    # bl = 5000000
    # sf = 5000
    # trace_mcmc_gene_tree_errors(dir,true_gt_path, bl, sf)
    # # trace_mcmc_gene_tree_errors(dir,true_gt_path, bl, sf, "G 0")
    # # trace_mcmc_gene_tree_errors(dir,true_gt_path, bl, sf)
    #
    # trace_mcmc_murates_error(dir, true_mu_path)

    SEED_INDEX = 0
    # dir_path="/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/net0/long/"
    # rep=1
    # sub_dir=dir_path+str(rep)+"/heter/500/"
    # burn_in_len = 2000000
    # samplefreq = 500
    # burnin_prop=80
    # compare_true_mcmc_inferred_murates(sub_dir+"/murate.txt", sub_dir+"murate/", burnin=70)
    # mcmc_murates_error(sub_dir+"murate/", sub_dir+"/murate.txt", burnin=70)
    # trace_mcmc_inferred_murates_each_locus(sub_dir+"murate/")
    mcmc_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/net0/long/7/heter/2000/"
    true_gt_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/net0/long/7/homo/2000/genetrees_pruned.txt"
    # trace_mcmc_gene_tree_errors_2(mcmc_dir, true_gt_path, 2000000, 5000, burnin=50)
    trace_mcmc_gene_tree_errors_2(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), burnin=int(sys.argv[5]))


    # true mu path
    # mcmc dir
    # true gt path
    # bl
    # sf
    # burn in : optional
    
    
    
    # if len(sys.argv) == 6:
    #     compare_true_mcmc_inferred_murates(sys.argv[1], sys.argv[2], burnin)
    #     mcmc_murates_error(sys.argv[2], sys.argv[1], burnin=70)
    #     trace_mcmc_inferred_murates_each_locus(sys.argv[2], burnin=70)
    #     trace_mcmc_gene_tree_errors(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))
    # elif len(sys.argv) == 7:
    #     # compare_true_mcmc_inferred_murates(sys.argv[1], sys.argv[2])
    #     # trace_mcmc_inferred_murates(sys.argv[2])
    #     # mcmc_trace(sys.argv[2], int(sys.argv[4]), int(sys.argv[5]))
    #     trace_mcmc_gene_tree_errors(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]+" 0")
    #     mcmc_murates_error(sys.argv[2], sys.argv[1])
    # elif len(sys.argv) == 8:
    #     compare_true_mcmc_inferred_murates(sys.argv[1], sys.argv[2], int(sys.argv[6]))
    #     # trace_mcmc_inferred_murates(sys.argv[2], int(sys.argv[6]))
    #     # mcmc_trace(sys.argv[2], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]))
    #     trace_mcmc_gene_tree_errors(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]))
    #     mcmc_murates_error(sys.argv[2], sys.argv[1], int(sys.argv[6]))
    # trace_mcmc_gene_tree_errors(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
