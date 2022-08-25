import dendropy
import pandas as pd
from summarize_stats_5taxa import GT_DICT, ST_DICT
from summarize_probs import make_canonical
from itertools import combinations
from calculator import *
import sys
from dendropy.calculate import treecompare



SPECIES = ["C 0", "G 0", "L 0", "Q 0", "R 0"]
SUBSET_SIZE = 3
DEBUG = False
if DEBUG:
    NUM_REPLICATE = 3
else:
    NUM_REPLICATE = 100


def divide_triplet():
    return list(combinations(SPECIES, SUBSET_SIZE))

TRIPLETS = divide_triplet()
TRUE_ST = "(((Q_0:0.3,R_0:0.3)I3:0.5,L_0:0.8)I1:0.2,(G_0:0.4,C_0:0.4)I2:0.6)I0;"

# use full gt topologies as both observed and expected
def subtree_freq(gt_freq_df, triplet):
    print(triplet)
    remove_taxa_list = [item for item in SPECIES if item not in triplet]
    triplet_gt_obs_freq_dict = {}
    triplet_gt_exp_freq_dict = {}

    for i, row in gt_freq_df.iterrows():
        tree = dendropy.Tree.get(data=row["gt"], schema="newick", rooting="force-rooted")
        tree.prune_taxa_with_labels(remove_taxa_list)
        make_canonical(tree.seed_node)
        trip_str = tree.as_string("newick", suppress_rooting = True).strip()

        if trip_str in triplet_gt_obs_freq_dict:
            triplet_gt_obs_freq_dict[trip_str] += row["prob_obs"]
        else:
            triplet_gt_obs_freq_dict[trip_str] = row["prob_obs"]
        
        if trip_str in triplet_gt_exp_freq_dict:
            triplet_gt_exp_freq_dict[trip_str] += row["prob_exp"]
        else:
            triplet_gt_exp_freq_dict[trip_str] = row["prob_exp"]
        
    df_obs = pd.DataFrame(triplet_gt_obs_freq_dict.items(), columns=['gt', "prob_obs"])
    df_exp = pd.DataFrame(triplet_gt_exp_freq_dict.items(), columns=['gt', "prob_exp"])
    prob_df = df_obs.merge(df_exp, on='gt', how="outer")
    prob_df['prob_obs'] = prob_df['prob_obs'].fillna(0)
    prob_df = prob_df.sort_values(by="prob_obs", ascending=False)

    major_tree = dendropy.Tree.get(data=prob_df.at[0, "gt"].strip(),schema="newick", rooting="force-rooted")
    sub_leaf_list = [node.taxon.label for node in major_tree.leaf_node_iter()]
    true_tree = dendropy.Tree.get(data=TRUE_ST,schema="newick", rooting="force-rooted")
    full_leaf_list = [node.taxon.label for node in true_tree.leaf_node_iter()]
    leaf_to_remove = list(set(full_leaf_list)-set(sub_leaf_list))
    true_tree.prune_taxa_with_labels(leaf_to_remove)
    true_sub_tree = dendropy.Tree.get(data=true_tree.as_string(schema="newick"), schema="newick", taxon_namespace=major_tree.taxon_namespace)
    if treecompare.symmetric_difference(true_sub_tree, major_tree) == 0:
        st = True
    else:
        st = False

    return prob_df, st


# use full gt observations but computed expectations by equation
def subtree_obs_freq(gt_obs_freq_df, triplet):
    print(triplet)
    remove_taxa_list = [item for item in SPECIES if item not in triplet]
    triplet_gt_obs_freq_dict = {}

    for i, row in gt_obs_freq_df.iterrows():
        tree = dendropy.Tree.get(data=row["gt"], schema="newick", rooting="force-rooted")
        tree.prune_taxa_with_labels(remove_taxa_list)
        make_canonical(tree.seed_node)
        trip_str = tree.as_string("newick", suppress_rooting = True).strip()

        if trip_str in triplet_gt_obs_freq_dict:
            triplet_gt_obs_freq_dict[trip_str] += row["prob_obs"]
        else:
            triplet_gt_obs_freq_dict[trip_str] = row["prob_obs"] 
        
    df_probs = pd.DataFrame(triplet_gt_obs_freq_dict.items(), columns=['gt', "prob_obs"])
    df_probs = df_probs.sort_values(by="prob_obs", ascending=False)
    print(df_probs)
    # check major_tree topology
    major_tree = dendropy.Tree.get(data=df_probs.at[0, "gt"].strip(),schema="newick", rooting="force-rooted")
    sub_leaf_list = [node.taxon.label for node in major_tree.leaf_node_iter()]
    true_tree = dendropy.Tree.get(data=TRUE_ST,schema="newick", rooting="force-rooted")
    full_leaf_list = [node.taxon.label for node in true_tree.leaf_node_iter()]
    leaf_to_remove = list(set(full_leaf_list)-set(sub_leaf_list))
    true_tree.prune_taxa_with_labels(leaf_to_remove)
    true_sub_tree = dendropy.Tree.get(data=true_tree.as_string(schema="newick"), schema="newick", taxon_namespace=major_tree.taxon_namespace)
    if treecompare.symmetric_difference(true_sub_tree, major_tree) == 0:
        st = True
    else:
        st = False
    print(df_probs.at[0, "prob_obs"])
    t = internal_time(df_probs.at[0, "prob_obs"])
    exp = compute_exp(t)
    df_probs["prob_exp"] = exp
    return df_probs, st

def compute_triplet_probs(prob_csv_path, triplet_csv_path, exp_type=True):
    # triplets = divide_triplet()
    df_prob = pd.read_csv(prob_csv_path)
    df_list = []
    st_dict = {"triplet":[], "st_infer":[]}
    for triplet in TRIPLETS:
        if exp_type:
            tri_df, st = subtree_obs_freq(df_prob, triplet)
        else:
            tri_df, st = subtree_freq(df_prob, triplet)
        
        tri_df["triplet"] = [triplet for x in range(len(tri_df["gt"]))]
        df_list.append(tri_df)
        st_dict["triplet"].append(triplet)
        st_dict["st_infer"].append(st)
    res_df = pd.concat(df_list, ignore_index=True)
    res_df.to_csv(triplet_csv_path, index=False)
    # print(res_df)
    # print(pd.DataFrame(st_dict))
    return pd.DataFrame(st_dict)

def run_all(dir_path, re_end, exp_type, null_reti_num):
    # dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/"
    re_start = 1
    # re_end = 1
    # exp_type = True
    # divide_triplet()
    st_list = []
    for scale in ["short", "medium", "long"]:
        for locus_length in [500, 2000]:
            for i in range(re_start, re_end+1):
                for true_gt in [True, False]:
                    if exp_type:
                        if true_gt:
                            prob_path = dir_path+scale+"/"+str(i)+"/homo/"+str(locus_length)+"/gt_prob_obs_sim.csv"
                            triplet_path =  dir_path+scale+"/"+str(i)+"/heter/"+str(locus_length)+"/triplet_probs_sim.csv"
                        else:
                            prob_path = dir_path+scale+"/"+str(i)+"/heter/"+str(locus_length)+"/gt_prob_obs_iqtree.csv"
                            triplet_path =  dir_path+scale+"/"+str(i)+"/heter/"+str(locus_length)+"/triplet_probs_iqgt.csv"
                    else:
                        true_st = False
                        prob_path = dir_path+scale+"/"+str(i)+"/heter/"+str(locus_length)+"/gt_prob_all_"+ST_DICT[true_st]+"_"+GT_DICT[true_gt]+"_"+null_reti_num+".csv"
                        triplet_path =  dir_path+scale+"/"+str(i)+"/heter/"+str(locus_length)+"/triplet_probs_"+ST_DICT[true_st]+"_"+GT_DICT[true_gt]+"_"+null_reti_num+".csv"
                    st_infer_df = compute_triplet_probs(prob_path, triplet_path, exp_type)
                    st_infer_df["scale"] = [scale for x in range(len(st_infer_df["st_infer"]))]
                    st_infer_df["locus_length"] = [locus_length for x in range(len(st_infer_df["st_infer"]))]
                    st_infer_df["replicate"] = [i for x in range(len(st_infer_df["st_infer"]))]
                    st_infer_df["true_gt"] = [true_gt for x in range(len(st_infer_df["st_infer"]))]
                    st_list.append(st_infer_df)
    res_df = pd.concat(st_list, ignore_index=True)
    res_df.to_csv(dir_path+"/triplet_st_infer_exp"+str(exp_type)+".csv")


if __name__ == "__main__":
    # obs_csv_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/long/1/heter/500/gt_prob_all_truest_iqgt_0.csv"
    # compute_triplet_probs(obs_csv_path, "./test.csv", True)
    dir_path = sys.argv[1]
    re_end = int(sys.argv[2])
    exp_type = (sys.argv[3].lower() == "true")
    null_reti_num = sys.argv[4]
    # dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/"
    # re_end = 2
    # exp_type = False
    run_all(dir_path, re_end, exp_type, null_reti_num)