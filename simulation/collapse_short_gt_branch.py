import dendropy
import random
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from dendropy.calculate import treecompare


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


def collapse_short_gt_branch(tree_path, out_tree_path, rng, prune_taxa_list=['G 0']):
    # fig, axes = plt.subplots(1, 1, figsize=(7, 5))
    tree_list = dendropy.TreeList.get(path=tree_path, schema="newick", rooting="force-rooted")

    # resolved_list = dendropy.TreeList()
    
    # tree_0 = dendropy.Tree.get(data= "(Z_0,(C_0:0, G_0:0, L_0:0, R_0:0, Q_0:0));", schema="newick", rooting="force-rooted", taxon_namespace=tree_list[0].taxon_namespace)

    
    cnt = 0
    threshold=1e-5
    for i, tree in enumerate(tree_list):
        # print(tree.as_string(schema='newick'))
        
        collapsed = False
        tree.prune_taxa_with_labels(prune_taxa_list)
        for e in tree.postorder_edge_iter():
            if e.length is not None and ((e.length <= threshold) and e.is_internal()):
               e.collapse()
               collapsed = True
        
        # tree.collapse_unweighted_edges(threshold=1e-5)
        if collapsed:
            # print(tree.as_string(schema='newick'))
            # if treecompare.symmetric_difference(tree, tree_0) < 0.001:
            #     cnt += 1
            #     resolved_list.append(tree)
            make_canonical(tree.seed_node)
            tree.encode_bipartitions()
            tree.resolve_polytomies(limit=2, update_bipartitions=True, rng=rng)
            # print(tree.as_string(schema='newick'))
            # print("------------")
            # print(len(tree.seed_node._child_nodes[0]._child_nodes), len(tree.seed_node._child_nodes[1]._child_nodes))
            
    with open(out_tree_path, "w") as handle:
        s = ""
        for tree in tree_list:
            s += tree.as_string(schema='newick')[5:]
            # s += "\n"
        handle.write(s)
   
def prune_gt(tree_path, out_tree_path, prune_taxa_list=['G 0']):
    tree_list = dendropy.TreeList.get(path=tree_path, schema="newick", rooting="force-rooted")
    for i, tree in enumerate(tree_list):
        tree.prune_taxa_with_labels(prune_taxa_list)
        
            
    with open(out_tree_path, "w") as handle:
        s = ""
        for tree in tree_list:
            s += tree.as_string(schema='newick')[5:]
        handle.write(s)

def run_all(dir_path, num_replicate, scale, heter=True, gt=False):
    rng = random.Random(12345)
    # for scale in ["short", "medium", "long"]:
    locus_length=2000
    for i in range(1, num_replicate+1):
        # for locus_length in [500, 2000]:
        if gt == False:
            if heter:
                iqtree_path = dir_path + scale+ "/"+str(i)+"/heter/"+str(locus_length)+"/rooted_iqtree.txt"
                resolved_tree_path = dir_path + scale+ "/"+str(i)+"/heter/"+str(locus_length)+"/rooted_iqtree_resolved.txt"
            else:
                iqtree_path = dir_path + scale+ "/"+str(i)+"/homo/"+str(locus_length)+"/rooted_iqtree.txt"
                resolved_tree_path = dir_path + scale+ "/"+str(i)+"/homo/"+str(locus_length)+"/rooted_iqtree_resolved.txt"
            collapse_short_gt_branch(iqtree_path, resolved_tree_path, rng)
        elif heter == False:
            truetree_path = dir_path + scale+ "/"+str(i)+"/homo/"+str(locus_length)+"/genetrees.txt"
            pruned_tree_path = dir_path + scale+ "/"+str(i)+"/homo/"+str(locus_length)+"/genetrees_pruned.txt"
            prune_gt(truetree_path, pruned_tree_path)

if __name__ == "__main__":
    # tree_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/long/1/heter/500/rooted_iqtree.txt"
    # out_tree_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/long/1/heter/500/rooted_iqtree_resolved.txt"
 
    # collapse_short_gt_branch(tree_path, out_tree_path)

    run_all(sys.argv[1], int(sys.argv[2]), sys.argv[3], (sys.argv[4].lower() == "true"), (sys.argv[5].lower() == "true"))
