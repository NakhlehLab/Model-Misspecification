
import sys
import pandas as pd
import dendropy
from dendropy.calculate import treecompare
import matplotlib.pyplot as plt
import seaborn as sns

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

def gt_parital_info_distribution(gt_dir, num_replicate, locus_length):
    iqtree_list = []
    resolvetree_list = []
    partial_informative_list_iq = dendropy.TreeList()
    partial_informative_list_resolve = dendropy.TreeList()
    
    tns = None
    threshold =1e-5
    for i in range(1, num_replicate+1):
        iq_gt_path = gt_dir + str(i) + "/heter/"+str(locus_length)+"/rooted_iqtree.txt"
        resolve_gt_path = gt_dir + str(i) + "/heter/"+str(locus_length)+"/rooted_iqtree_resolved.txt"

        if tns is None:
            replicate_iqtree_list = dendropy.TreeList.get(path=iq_gt_path, schema="newick", rooting="force-rooted")            
            tns = replicate_iqtree_list[0].taxon_namespace
        else:
            replicate_iqtree_list = dendropy.TreeList.get(path=iq_gt_path, schema="newick", rooting="force-rooted", taxon_namespace=tns)
            
        iqtree_list.extend(replicate_iqtree_list)    
        replicate_resolve_list = dendropy.TreeList.get(path=resolve_gt_path, schema="newick", rooting="force-rooted", taxon_namespace=tns)
        resolvetree_list.extend(replicate_resolve_list)
    for tree1, tree2 in zip(iqtree_list, resolvetree_list):
        collapsed = False
        for e in tree1.postorder_edge_iter():
            if e.length is not None and ((e.length <= threshold) and e.is_internal()):
               collapsed = True
               break
        if collapsed:
            partial_informative_list_iq.append(tree1)
            partial_informative_list_resolve.append(tree2)

    topologies_iq = partial_informative_list_iq.as_tree_array().topologies(sort_descending = True, frequency_attr_name='frequency')
    topologies_resolve = partial_informative_list_resolve.as_tree_array().topologies(sort_descending = True, frequency_attr_name='frequency')

    # freqs1 = []
    # freqs2 = []
    # topo_freq_dict = {"gt": [], "iq_freq":[], "resolved_freq":[]}
    topo_freq_dict = {}
    for tree in topologies_iq:
        make_canonical(tree.seed_node)
        tree_newick = tree.as_string(schema = "newick", suppress_rooting = True).rstrip()
        topo_freq_dict[tree_newick]=[tree.frequency,0]

    for tree in topologies_resolve:
        make_canonical(tree.seed_node)
        tree_newick = tree.as_string(schema = "newick", suppress_rooting = True).rstrip()
        if tree_newick in topo_freq_dict.keys():
            topo_freq_dict[tree_newick][1] = tree.frequency
        else:
            topo_freq_dict[tree_newick] = [0, tree.frequency]
    
    df_topo_freq = pd.DataFrame.from_dict(topo_freq_dict, orient="index",
                       columns=['iq_freq', 'resolve_freq'])
    df_topo_freq.reset_index(inplace=True)
    df_topo_freq = df_topo_freq.rename(columns = {'index':'gt'})
    df_topo_freq.to_csv(gt_dir+"gt_info_dist_"+str(locus_length)+".csv")
    return df_topo_freq
    # for tree1 in topologies_iq:
    #     for tree2 in topologies_resolve:
    #         if treecompare.symmetric_difference(tree1, tree2) <= 0.0001:
    #             tree1.append
    # for tree in topologies_iq:
    #     make_canonical(tree.seed_node)
    #     tree_newick = tree.as_string(schema = "newick", suppress_rooting = True).rstrip()
    #     topo_freq_dict["gt"].append(tree_newick)
    #     topo_freq_dict["iq_freq"].append(tree.frequency)
    #     freqs1.append(tree.frequency)

    # for tree in topologies_resolve:
    #     make_canonical(tree.seed_node)
    #     tree_newick = tree.as_string(schema = "newick", suppress_rooting = True).rstrip()
    #     if tree_newick in topo_freq_dict["gt"]:
            
    #         topo_freq_dict["gt"].append(tree_newick)
    #     else:
    #         topo_freq_dict["gt"].append(tree_newick)
    #         topo_freq_dict["iq_freq"].append(0)
    #     freqs2.append(tree.frequency)

    # for tree1 in topologies_iq:
    #     for tree2 in topologies_resolve:
    #         if treecompare.symmetric_difference(tree1, tree2) <= 0.0001:
    #             tree1.append
    # return freqs1, freqs2

        
def run_all_gt_info(dir_path, num_replicate):
    number_of_genes = 105
    for locus_length in [500, 2000]:
        plt.clf()
        fig, axes = plt.subplots(2, 3, figsize=(20, 5))
        for i, scale in enumerate(["short", "medium", "long"]):
            path = dir_path + scale + "/"
            df_freqs = gt_parital_info_distribution(path, num_replicate, locus_length)
            # dfm = df_freqs.melt('gt', var_name='gt_type', value_name='freqs')
            # sns.barplot(x="gt", y="freqs",hue = "gt_type", data=dfm, ax=axes[i])

            sns.barplot(x="gt", y="iq_freq", data=df_freqs, ax=axes[0][i])
            sns.barplot(x="gt", y="resolve_freq", data=df_freqs, ax=axes[1][i])
            for x in [0,1]:
                axes[x][i].set_ylim(0, 0.3)
                axes[x][i].set_ylabel("Probability", fontsize=16)
                axes[x][i].set_xlabel("Gene Tree", fontsize=16)
                axes[x][i].set_title(scale, fontsize=16, fontweight="bold")
                axes[x][i].set_xticks([])
        plt.tight_layout()
        # plt.tick_params(
        #     axis='x',          # changes apply to the x-axis
        #     which='both',      # both major and minor ticks are affected
        #     bottom=False,      # ticks along the bottom edge are off
        #     top=False,         # ticks along the top edge are off
        #     labelbottom=False)
        plt.savefig(dir_path+"/gt_info_distribution_"+str(locus_length) + ".pdf")


if __name__ == "__main__":
    # plot_polymorphism()
    run_all_gt_info(sys.argv[1], int(sys.argv[2]))
    # plot_gt_non_info_distribution()
    # plot_non_info_sites_cnt()
