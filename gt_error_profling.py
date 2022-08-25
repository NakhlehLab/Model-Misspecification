import scipy.stats
import dendropy
from dendropy.calculate import treecompare
import sys
import os
import pandas as pd
from sympy import true


def get_trees_heights(gene_tree_path, murate_path):
    gt_list = dendropy.TreeList.get(path=gene_tree_path, schema="newick", rooting="force-rooted")
    murate_list = get_murates(murate_path)
    distances = []
    for i in range(len(gt_list)):
        tree = gt_list[i]
        scale = murate_list[i]
        distances.append(tree.calc_node_root_distances(return_leaf_distances_only=True)[0] * scale)
    # print(distances)
    return distances

def get_murates(murate_path):
    murate_list = []
    with open(murate_path, "r") as handle:
        for line in handle.readlines():
            murate_list.append(float(line.strip()))
    return murate_list

def heterogeneity_rate(gene_tree_path, murate_path):

    root_to_tip_distances = get_trees_heights(gene_tree_path, murate_path)
    cov = scipy.stats.variation(root_to_tip_distances)
    print("The CoV if mutation rate * gene tree height is： ", cov)
    return cov

def murate_level(murate_path):
    murates_list = get_murates(murate_path)
    cov = scipy.stats.variation(murates_list)
    print("The CoV of mutation rate is: ", cov)
    return cov

def tree_ingroup_rf_distance(truetree: dendropy.Tree, tree2: dendropy.Tree, OUTGROUP=None):
    truetree.encode_bipartitions()
    tree2.encode_bipartitions()

    if OUTGROUP is not None:
        truetree.prune_taxa_with_labels([OUTGROUP])
        tree2.prune_taxa_with_labels([OUTGROUP])
    distance_ingroup = treecompare.symmetric_difference(truetree, tree2)
    # print(distance_ingroup)
    return distance_ingroup

def tree_ingroup_unroot_rf_distance(truetree: dendropy.Tree, tree2: dendropy.Tree, OUTGROUP=None):
    if OUTGROUP is not None:
        truetree.prune_taxa_with_labels([OUTGROUP])
        tree2.prune_taxa_with_labels([OUTGROUP])
    truetree.encode_bipartitions()
    tree2.encode_bipartitions()
    unroot_truetree = dendropy.Tree.get(data=truetree.as_string(schema="newick"), schema="newick", rooting="force-unrooted", taxon_namespace=truetree.taxon_namespace)
    unroot_iqtree = dendropy.Tree.get(data=tree2.as_string(schema="newick"), schema="newick", rooting="force-unrooted", taxon_namespace=truetree.taxon_namespace)
    distance_ingroup = treecompare.symmetric_difference(unroot_truetree, unroot_iqtree)
    return distance_ingroup


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

def compare_iqtree_truetree(truetree_path, iqtree_path, outgroup, scale_true_gt=0.01):
    tns = dendropy.TaxonNamespace()
    iqtree_list = dendropy.TreeList.get(path=iqtree_path, schema="newick", rooting="force-rooted", taxon_namespace=tns)
    truetree_list = dendropy.TreeList.get(path=truetree_path, schema="newick", rooting="force-rooted",
                                          taxon_namespace=tns)
    true_cnt = 0
    rf_ingroup_distance = 0
    rf_outgroup_distance = 0
    wrf_distance = 0
    nrBS = 0
    num_taxa = len(iqtree_list[0].taxon_namespace)
    for iqtree, truetree in zip(iqtree_list, truetree_list):
        truetree = tree_scale(truetree, scale_true_gt)
        
        rf_o = treecompare.symmetric_difference(iqtree, truetree)
        rf_outgroup_distance += rf_o/(2*(num_taxa-2))
       
        rf = tree_ingroup_rf_distance(truetree, iqtree, outgroup)
        rf_ingroup_distance += rf/(2*(num_taxa-2))
        if rf < 0.01:
            true_cnt += 1
        nrBS += tree_ingroup_euclidean_distance(truetree, iqtree, outgroup)
    rf_ingroup_distance /= len(iqtree_list)
    rf_outgroup_distance /= len(iqtree_list)
    
    wrf_distance /= len(iqtree_list)
    nrBS /= len(iqtree_list)
    print("Number of correctly inferred gene trees: ", true_cnt)
    print("Average normalized ingroup RF distance: ", rf_ingroup_distance)
    print("Average normalized outgroup RF distance: ", rf_outgroup_distance)
    # print("Average weighted RF distance: ", wrf_distance)
    print("Average nrBS of ingroup: ", nrBS)

def write_all_rfs(truetree_path, iqtree_path, outgroup, csv_path, scale_true_gt=0.01):
    tns = dendropy.TaxonNamespace()
    iqtree_list = dendropy.TreeList.get(path=iqtree_path, schema="newick", rooting="force-rooted", taxon_namespace=tns)
    truetree_list = dendropy.TreeList.get(path=truetree_path, schema="newick", rooting="force-rooted",
                                          taxon_namespace=tns)
    rf_outgroup_list = []
    rf_ingroup_list = []
    rf_ingroup_unroot_list = []
    nrBS_list = []

    num_taxa = len(iqtree_list[0].taxon_namespace)
    
    df = pd.DataFrame()
    for iqtree, truetree in zip(iqtree_list, truetree_list):
        truetree = tree_scale(truetree, scale_true_gt)

        rf_o = treecompare.symmetric_difference(truetree, iqtree)
        rf_outgroup_list.append(rf_o/(2*(num_taxa-2)))
    
        rf = tree_ingroup_rf_distance(truetree, iqtree, outgroup)
        rf_ingroup_list.append(rf/(2*(num_taxa-3)))

        nrBS_list.append(tree_ingroup_euclidean_distance(truetree, iqtree, outgroup))

        rf_unroot_in = tree_ingroup_unroot_rf_distance(truetree, iqtree, outgroup)
        rf_ingroup_unroot_list.append(rf_unroot_in/(2*(num_taxa-1-3)))
    
    df["RF-out"] = rf_outgroup_list
    df["RF-in"] = rf_ingroup_list
    df["nrBS"] = nrBS_list
    df["RF-in_unroot"] = rf_ingroup_unroot_list
    print(df)
    df.to_csv(csv_path, index=False)
    
def test_unroot_ingroup():
    tns = dendropy.TaxonNamespace()
    tree1 = dendropy.Tree.get(data="((((A,B),C),D),O);", schema="newick", taxon_namespace=tns) 
    tree2 = dendropy.Tree.get(data="((((A,C),B),D),O);", schema="newick", taxon_namespace=tns) 
    dist = tree_ingroup_rf_distance(tree1, tree2, "O")
    print(dist)

if __name__ == '__main__':
    true_tree_path = sys.argv[1]
    # murate_path = sys.argv[2]
    iq_tree_path = sys.argv[2]
    OUTGROUP = sys.argv[3]
    csv_path = sys.argv[4]
    poprate = float(sys.argv[5])

    # if os.path.exists(murate_path):
    #     heterogeneity_rate(true_tree_path, murate_path)
    #     murate_level(murate_path)
    # compare_iqtree_truetree(true_tree_path, iq_tree_path, OUTGROUP+" 0")
    write_all_rfs(true_tree_path, iq_tree_path, OUTGROUP+" 0", csv_path, poprate)

    # test_unroot_ingroup()