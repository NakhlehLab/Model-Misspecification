import pandas as pd
import dendropy
from pathlib import Path
import sys


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


def write_obs(tree_path, true_gt, outgroup=None):
    tree_list = dendropy.TreeList.get(path=tree_path, schema="newick", rooting="force-rooted")
    pruned_tree_list = []
    for tree in tree_list:
        if outgroup is not None:
            tree.prune_taxa_with_labels([outgroup])
        pruned_tree_list.append(tree)
    pruned_tree_list = dendropy.TreeList(pruned_tree_list)
    unique_topologies = pruned_tree_list.as_tree_array().topologies(sort_descending = True, frequency_attr_name='frequency')
    observed_frequencies = {"gt":[], "prob_obs":[]}
    for topology in unique_topologies:
        make_canonical(topology.seed_node)
        topology_newick = topology.as_string(schema = "newick", suppress_rooting = True).rstrip()
        observed_frequencies["gt"].append(topology_newick)
        observed_frequencies["prob_obs"].append(topology.frequency)
    df = pd.DataFrame(observed_frequencies)
    df.drop(df.columns[df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
    if true_gt:
        df.to_csv(str(Path(tree_path).parent.absolute())+"/gt_prob_obs_sim.csv", index=False)
    else:
        df.to_csv(str(Path(tree_path).parent.absolute())+"/gt_prob_obs_iqtree.csv", index=False)


def get_obs(csv_path):
    df = pd.read_csv(csv_path)
    return df

def get_exp(csv_path):
    df = pd.read_csv(csv_path)
    for i, row in df.iterrows():
        topology = dendropy.Tree.get(data=row["gt"], schema="newick", rooting="force-rooted")
        make_canonical(topology.seed_node)
        topology_newick = topology.as_string(schema = "newick", suppress_rooting = True).rstrip()
        df.at[i,'gt'] = topology_newick
    print(df)
    return df

def write_obs_exp_probs(obs_path, exp_path, prob_exp):
    obs_df = get_obs(obs_path)
    exp_df = get_exp(exp_path)
    prob_df = obs_df.merge(exp_df, on='gt', how="outer")
    print(prob_df)
    prob_df['prob_obs'] = prob_df['prob_obs'].fillna(0)
    prob_df.drop(prob_df.columns[prob_df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
    prob_df.to_csv(prob_exp, index=False)
    return prob_df


def summarize_obs_probs(dir_path, re_end):
    # dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/"
    # true_gt = True
    re_start = 1
    # re_end = 1
    outgroup = "Z 0"
    for scale in ["short", "medium", "long"]:
        for locus_length in [500, 2000]:
            for i in range(re_start, re_end+1):
                for true_gt in [True, False]:
                    if true_gt:
                        obs_path = dir_path+scale+"/"+str(i)+"/homo/"+str(locus_length)+"/genetrees.txt" 
                    else:
                        obs_path = dir_path+scale+"/"+str(i)+"/heter/"+str(locus_length)+"/rooted_iqtree_resolved.txt" 
                    write_obs(obs_path, true_gt, outgroup)
                    print(obs_path)

def summarize_probs(dir_path, re_end, null_num_reti):
    # dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/"
    re_start = 1
    # re_end = 10
    # true_st = False
    # true_gt = False
    for scale in ["short", "medium", "long"]:
        for locus_length in [500, 2000]:
            for i in range(re_start, re_end+1):
                for true_st in [False]:
                    for true_gt in [True, False]:
                        if true_gt:
                            obs_path = dir_path+scale+"/"+str(i)+"/homo/"+str(locus_length)+"/gt_prob_obs_sim.csv" 
                        else:
                            obs_path = dir_path+scale+"/"+str(i)+"/heter/"+str(locus_length)+"/gt_prob_obs_iqtree.csv" 
                        
                        if true_st:
                            exp_path = dir_path+scale+"/gt_prob_exp_true_st.csv"
                            if true_gt:
                                prob_path = dir_path+scale+"/"+str(i)+"/heter/"+str(locus_length)+"/gt_prob_all_truest_truegt_"+null_num_reti+".csv" 
                            else:
                                prob_path = dir_path+scale+"/"+str(i)+"/heter/"+str(locus_length)+"/gt_prob_all_truest_iqgt_"+null_num_reti+".csv" 
                        else:
                            if true_gt:
                                exp_path = dir_path+scale+"/"+str(i)+"/heter/"+str(locus_length)+"/gt_prob_exp_true_"+null_num_reti+".csv" 
                                prob_path = dir_path+scale+"/"+str(i)+"/heter/"+str(locus_length)+"/gt_prob_all_mlst_truegt_"+null_num_reti+".csv" 
                            else:
                                exp_path = dir_path+scale+"/"+str(i)+"/heter/"+str(locus_length)+"/gt_prob_exp_iq_"+null_num_reti+".csv" 
                                prob_path = dir_path+scale+"/"+str(i)+"/heter/"+str(locus_length)+"/gt_prob_all_mlst_iqgt_"+null_num_reti+".csv" 
                        write_obs_exp_probs(obs_path, exp_path, prob_path)


if __name__ == "__main__":
    dir_path = sys.argv[1]
    num_replicate = int(sys.argv[2])
    null_reti_num = sys.argv[3]

    summarize_obs_probs(dir_path, num_replicate)
    summarize_probs(dir_path, num_replicate, null_reti_num)