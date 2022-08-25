import dendropy
import pandas as pd
import tqdm
from itertools import combinations
import sys
from dendropy.calculate import treecompare

# with outgroup Z 0 as unrooted tree
SPECIES = ["C 0", "G 0", "L 0", "Q 0", "R 0", "Z 0"]
OUTGROUP = ["Z 0"]
INGROUP = list(set(SPECIES) - set(OUTGROUP))
SUBSET_SIZE = 4
DEBUG = False
if DEBUG:
    NUM_REPLICATE = 3
else:
    NUM_REPLICATE = 100
NUM_LOCI = 10000

def divide_quartet(species_list):
    return list(combinations(species_list, SUBSET_SIZE))

# QUARTET = divide_quartet()
TRUE_ST = "(((Q_0:0.3,R_0:0.3)I3:0.5,L_0:0.8)I1:0.2,(G_0:0.4,C_0:0.4)I2:0.6)I0;"

def check_tree_in_list(tree0, tree_list):
    for tree in tree_list:
        if treecompare.symmetric_difference(tree0, tree) <= 0.00001:
            return tree
    return None


# use full gt topologies as both observed and expected
def quartet_obs_freq(gt_freq_df, quartet):
    
    remove_taxa_list = [item for item in SPECIES if item not in quartet]
    triplet_gt_obs_freq_dict = {}
    
    tns = None
    for i, row in gt_freq_df.iterrows():
        if tns is not None:
            tree = dendropy.Tree.get(data=row["gt"], schema="newick", rooting="force-unrooted", taxon_namespace=tns)
        else:
            tree = dendropy.Tree.get(data=row["gt"], schema="newick", rooting="force-unrooted")
            tns = tree.taxon_namespace
        tree.prune_taxa_with_labels(remove_taxa_list)
        tree.update_bipartitions()
        trip_str = tree.as_string("newick", suppress_rooting = True).strip()
            
        idendical_tree = check_tree_in_list(tree, triplet_gt_obs_freq_dict.keys())
        if idendical_tree is not None:
            triplet_gt_obs_freq_dict[idendical_tree] += row["prob_obs"]
        else:
            triplet_gt_obs_freq_dict[tree] = row["prob_obs"]
    return triplet_gt_obs_freq_dict


def compute_quartet_probs(tree_path, quartet_csv_path, quartet_list):
    tree_list = dendropy.TreeList.get(path=tree_path, schema="newick", rooting="force-unrooted")
    for tree in tree_list:
        tree.encode_bipartitions()
    unique_topologies = tree_list.as_tree_array().topologies(sort_descending = True, frequency_attr_name='frequency')
    observed_frequencies = {"gt":[], "prob_obs":[]}
    for topology in unique_topologies:
        topology.update_bipartitions()
        # make_canonical(topology.seed_node)
        topology_newick = topology.as_string(schema = "newick", suppress_rooting = True).rstrip()
        observed_frequencies["gt"].append(topology_newick)
        observed_frequencies["prob_obs"].append(topology.frequency)
    gt_obs_freq_df = pd.DataFrame(observed_frequencies)
    print(gt_obs_freq_df)
    # df_quartet = {"t1":[],"t2":[],"t3":[], "t4":[], "CF12_34":[], "CF13_24":[], "CF14_23":[], "ngenes":[]}
    quartet_obs_dict = {}
    df_list = []
    
    for quartet in quartet_list:
        print(quartet)
        topo_2_freq = quartet_obs_freq(gt_obs_freq_df, quartet)
        print(topo_2_freq)
        for i in range(4):
            label = quartet[i].split()[0]
            quartet_obs_dict["t"+str(i+1)] = label

        for k, v in topo_2_freq.items():
            pair = []
            for c in k.seed_node.child_nodes():
                if c.taxon is not None:
                    label = c.taxon.label.split()[0]
                    pair.append(label)
            print(pair)
            if (quartet_obs_dict["t1"] in pair and quartet_obs_dict["t2"] in pair) or (quartet_obs_dict["t3"] in pair and quartet_obs_dict["t4"] in pair):
                quartet_obs_dict["CF12_34"] = v
            elif (quartet_obs_dict["t1"] in pair and quartet_obs_dict["t3"] in pair) or (quartet_obs_dict["t2"] in pair and quartet_obs_dict["t4"] in pair):
                quartet_obs_dict["CF13_24"] = v
            elif (quartet_obs_dict["t1"] in pair and quartet_obs_dict["t4"] in pair) or (quartet_obs_dict["t2"] in pair and quartet_obs_dict["t3"] in pair):
                quartet_obs_dict["CF14_23"] = v
            else:
                raise Exception("wrong quartet splits")
            print(k.as_string(schema="newick"))
        
        print(quartet_obs_dict)
        df = pd.DataFrame(quartet_obs_dict, index=['i',])
        df_list.append(df)
    df_res = pd.concat(df_list, ignore_index=True)
    df_res["ngenes"] = [NUM_LOCI for i in range(len(df_res["t1"]))]
    print(df_res)
    df_res.to_csv(quartet_csv_path, index = False)

    
# def run_all(dir_path, re_end):
#     re_start = 1
#     for scale in ["short", "medium", "long"]:
#         for locus_length in [500, 2000]:
#             for i in range(re_start, re_end+1):
#                 for true_gt in [True, False]:
#                     if true_gt == True:
#                         gt_path = dir_path + scale +"/"+ str(i)+ "/homo/"+str(locus_length) + "/genetrees.txt"
#                         quartet_path = dir_path + scale +"/"+ str(i)+ "/heter/"+str(locus_length) + "/quartet_probs_true.csv"
#                     else:
#                         gt_path = dir_path + scale +"/"+ str(i)+ "/heter/"+str(locus_length) + "/rooted_iqtree_resolved.txt"
#                         quartet_path = dir_path + scale +"/"+ str(i)+ "/heter/"+str(locus_length) + "/quartet_probs_iq.csv"
#                     compute_quartet_probs(gt_path, quartet_path)

def run(dir_path, scale, locus_length, i, true_gt=False, in_group = False):
    if true_gt == True:
        gt_path = dir_path + scale +"/"+ str(i)+ "/homo/"+str(locus_length) + "/genetrees.txt"
        if in_group == True:
            quartet_path = dir_path + scale +"/"+ str(i)+ "/heter/"+str(locus_length) + "/quartet_probs_true_ingroup.csv"
        else:
            quartet_path = dir_path + scale +"/"+ str(i)+ "/heter/"+str(locus_length) + "/quartet_probs_true.csv"
    else:
        gt_path = dir_path + scale +"/"+ str(i)+ "/heter/"+str(locus_length) + "/rooted_iqtree_resolved.txt"
        if in_group == True:
            quartet_path = dir_path + scale +"/"+ str(i)+ "/heter/"+str(locus_length) + "/quartet_probs_iq_ingroup.csv"
        else:
            quartet_path = dir_path + scale +"/"+ str(i)+ "/heter/"+str(locus_length) + "/quartet_probs_iq.csv"
    if in_group == True:
        quartet_list = divide_quartet(INGROUP)
    else:
        quartet_list = divide_quartet(SPECIES)
    compute_quartet_probs(gt_path, quartet_path, quartet_list)

if __name__ == "__main__":
    # gt_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net1/long/1/heter/1000/iqtree.txt"
    # quartet_csv_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net1/long/1/heter/1000/quartet_probs.csv"
    # compute_quartet_probs(gt_path, quartet_csv_path)
    # dir_path = sys.argv[1]
    # re_end = int(sys.argv[2])
    # run_all(dir_path, re_end)

    run(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), (sys.argv[5].lower() == "true"), (sys.argv[6].lower() == "true"))