import pandas as pd
import dendropy
from itertools import combinations
from summarize_triplet_probs import make_canonical
from dendropy.calculate import treecompare

SPECIES = ["C 0", "G 0", "L 0", "Q 0", "R 0"]
OUTGROUP = ["Z 0"]
SUBSET_SIZE = 3
DEBUG = False
if DEBUG:
    NUM_REPLICATE = 3
else:
    NUM_REPLICATE = 100


def divide_triplet():
    return list(combinations(SPECIES, 3))

def divide_quartet():
    return list(combinations(SPECIES, 4))

TRIPLETS = divide_triplet()
QUARTETS = divide_quartet()

def check_iqtree_resolved(path):
    tree_list = dendropy.TreeList.get(path = path, schema="newick", rooting="force-rooted")
    for tree in tree_list:
        tree.prune_taxa_with_labels(OUTGROUP)
    triplet_gt_obs_freq_dict = {}
    gt_triplet_dict = {}
    for triplet in TRIPLETS:
        remove_taxa_list = [item for item in SPECIES if item not in triplet]
        for tree in tree_list:
            tree_cp = tree.clone()
            tree_cp.prune_taxa_with_labels(remove_taxa_list)
            find = False
            for tri_tree in triplet_gt_obs_freq_dict.keys():
                if treecompare.symmetric_difference(tri_tree, tree_cp) <= 0.001:
                    triplet_gt_obs_freq_dict[tri_tree] += 1
                    find = True
                    break
            if find == False:
                triplet_gt_obs_freq_dict[tree_cp] = 1
                gt_triplet_dict[tree_cp] = triplet


    df_triplet_probs = pd.DataFrame(triplet_gt_obs_freq_dict.items(), columns=['gt', "prob_obs"])
    df_gt_to_trip = pd.DataFrame(gt_triplet_dict.items(), columns=['gt', "triplet"])
    df_triplet_probs = df_triplet_probs.sort_values(by="prob_obs", ascending=False)    
    triplet_result = pd.merge(df_triplet_probs, df_gt_to_trip, on="gt")
    triplet_result["newick"] = [x.as_string(schema="newick", suppress_rooting = True).rstrip() for x in triplet_result["gt"]]   
    triplet_result = triplet_result.drop(["gt"], axis=1)
    print(triplet_result)

    


def check_quartet(iq_path, quartet_path):
    df = pd.read_csv(quartet_path)
    tree_list = dendropy.TreeList.get(path = iq_path, schema="newick", rooting="force-unrooted")

    quartet_gt_obs_freq_dict = {}
    gt_quartet_dict = {}
    for tree in tree_list:
        tree.is_rooted=False
        tree.encode_bipartitions()
    for quartet in QUARTETS:
        remove_taxa_list = [item for item in SPECIES if item not in quartet]
        for tree in tree_list:
            tree_cp = tree.clone()
            tree_cp.prune_taxa_with_labels(remove_taxa_list)
            
            find = False
            for tri_tree in quartet_gt_obs_freq_dict.keys():
                if treecompare.symmetric_difference(tri_tree, tree_cp) <= 0.001:
                    quartet_gt_obs_freq_dict[tri_tree] += 1
                    find = True
                    gt_quartet_dict[tri_tree] = quartet
                    break
            if find == False:
                quartet_gt_obs_freq_dict[tree_cp] = 1

    df_quartet_probs = pd.DataFrame(quartet_gt_obs_freq_dict.items(), columns=['gt', "prob_obs"])
    df_gt_to_quartet = pd.DataFrame(gt_quartet_dict.items(), columns=['gt', "quartet"])
    df_quartet_probs = df_quartet_probs.sort_values(by="prob_obs", ascending=False)    
    quartet_result = pd.merge(df_quartet_probs, df_gt_to_quartet, on="gt")
    quartet_result["newick"] = [x.as_string(schema="newick", suppress_rooting = True).rstrip() for x in quartet_result["gt"]]   
    quartet_result = quartet_result.drop(["gt"], axis=1)
    print(quartet_result)


    return quartet_result



def check_triplet(path):
    df = pd.read_csv(path)


def check_iq_probs(path_tree, path_df):
    tree_list = dendropy.TreeList.get(path = path_tree, schema="newick", rooting="force-rooted")
    tns = tree_list[0].taxon_namespace
    df = pd.read_csv(path_df)
    topology_list = list(df["gt"])
    
    dict_freqs = {}
    for topology in topology_list:
        dict_freqs[topology] = 0
    
    for tree1 in tree_list:
        tree1.prune_taxa_with_labels(OUTGROUP)
        for tree2_newick in topology_list:
            tree2 = dendropy.Tree.get(data=tree2_newick, schema="newick", rooting="force-rooted", taxon_namespace=tns)
            if treecompare.symmetric_difference(tree1, tree2) <= 0.001:
                dict_freqs[tree2_newick] += 1
                break
    for topology in topology_list:
        dict_freqs[topology] /= 10000

    df_freqs = pd.DataFrame(dict_freqs.items(), columns=['gt', "freqs"])
    df_res = pd.merge(df, df_freqs, on="gt")
    print(df_res)
    for index, row in df_res.iterrows():
        if abs(row["freqs"] - row["prob_obs"]) > 0.0001:
            print('wrong')
        

if __name__ == "__main__":
    triplet_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/long/3/heter/500/triplet_probs_iqgt.csv"
    quartet_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/long/3/heter/500/quartet_probs_iq.csv"
    resolved_iq_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/long/3/heter/500/rooted_iqtree_resolved.txt"
    # check_triplet(triplet_path)
    iq_probs_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/long/3/heter/500/gt_prob_obs_iqtree.csv"
    # check_iqtree_resolved(resolved_iq_path)
    # check_iq_probs(resolved_iq_path, iq_probs_path)
    check_quartet(resolved_iq_path, quartet_path)

   