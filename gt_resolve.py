# import collections
# import itertools

# import numpy as np
# import tskit

# def resolve_polytomy(parent_node_id, child_ids, new_nodes_by_time, rng):
#     """
#     For a polytomy and list of child node ids, return a list of (child, parent) tuples,
#     describing a bifurcating tree, rooted at parent_node_id, where the new_nodes_by_time
#     have been used to break polytomies. All possible topologies should be equiprobable.
#     """
#     assert len(child_ids) == len(new_nodes_by_time) + 2
#     edges = [[child_ids[0], None], ]  # Introduce a single edge that will be deleted later
#     edge_choice = rng.integers(0, np.arange(1, len(child_ids) * 2 - 1, 2))
#     tmp_new_node_lab = [parent_node_id] + new_nodes_by_time
#     assert len(edge_choice) == len(child_ids) - 1
#     for node_lab, child_id, target_edge_id in zip(tmp_new_node_lab, child_ids[1:], edge_choice):
#         target_edge = edges[target_edge_id]
#         # print("target", target_edge)
#         # Insert to keep edges in time order of parent
#         edges.insert(target_edge_id, [child_id, node_lab])
#         edges.insert(target_edge_id, [target_edge[0], node_lab])
#         target_edge[0] = node_lab
#     # We need to re-map the internal nodes so that they are in time order
#     real_node = iter(new_nodes_by_time)
#     edges.pop() # remove the unary node at the top
#     node_map = {c: c for c in child_ids}
#     # print("orig_edges", edges)
#     # print("parent IDs to allocate", new_nodes_by_time)
#     # last edge should have the highest node
#     node_map[edges[-1][1]] = parent_node_id
#     for e in reversed(edges):
#         # edges should be in time order - the oldest one can be give the parent_node_id
#         if e[1] not in node_map:
#             node_map[e[1]] = next(real_node)
#         if e[0] not in node_map:
#             node_map[e[0]] = next(real_node)
#         e[0] = node_map[e[0]]
#         e[1] = node_map[e[1]]
#     # print("mapped edges", edges)
#     assert len(node_map) == len(new_nodes_by_time) + len(child_ids) + 1
#     return edges


# def resolve_polytomies(ts, *, epsilon=1e-10, random_seed=None):
#     """
#     For a given parent node, an edge in or an edge out signifies a change in children
#     Each time such a change happens, we cut all existing edges with that parent,
#     and add the previous portion in to the new edge table. If, previously, there were
#     3 or more children for this node, we break the polytomy at random
#     """
#     rng = np.random.default_rng(seed=random_seed)

#     tables = ts.dump_tables()
#     edges_table = tables.edges
#     nodes_table = tables.nodes
#     # Store the left of the existing edges, as we will need to change it if the edge is split
#     existing_edges_left = edges_table.left
#     # Keep these arrays for handy reading later
#     existing_edges_right = edges_table.right
#     existing_edges_parent = edges_table.parent
#     existing_edges_child = edges_table.child
#     existing_node_time = nodes_table.time

#     edges_table.clear()

#     edges_for_node = collections.defaultdict(set)  # The edge ids dangling from each active node
#     nodes_changed = set()

#     for interval, e_out, e_in in ts.edge_diffs(include_terminal=True):
#         for edge in itertools.chain(e_out, e_in):
#             if edge.parent != tskit.NULL:
#                 nodes_changed.add(edge.parent)

#         pos = interval[0]
#         for parent_node in nodes_changed:
#             child_edge_ids = edges_for_node[parent_node]
#             if len(child_edge_ids) >= 3:
#                 # We have a previous polytomy to break
#                 parent_time = existing_node_time[parent_node]
#                 new_nodes = []
#                 child_ids = existing_edges_child[list(child_edge_ids)]
#                 remaining_edges = child_edge_ids.copy()
#                 left = None
#                 max_time = 0
#                 for edge_id, child_id in zip(child_edge_ids, child_ids):
#                     max_time = max(max_time, existing_node_time[child_id])
#                     if left is None:
#                         left = existing_edges_left[edge_id]
#                     else:
#                         assert left == existing_edges_left[edge_id]
#                     if existing_edges_right[edge_id] > interval[0]:
#                         # make sure we carry on the edge after this polytomy
#                         existing_edges_left[edge_id] = pos

#                 # ADD THE PREVIOUS EDGE SEGMENTS
#                 dt = min((parent_time - max_time)/(len(child_ids)*2), epsilon)
#                 # Each broken polytomy of degree N introduces N-2 extra nodes, each at a time
#                 # slighly less than the parent_time. Create new nodes in order of decreasing time
#                 new_nodes = [nodes_table.add_row(time=parent_time - (i * dt))
#                              for i in range(1, len(child_ids) - 1)]
#                 # print("new_nodes:", new_nodes, [tables.nodes[n].time for n in new_nodes])
#                 for new_edge in resolve_polytomy(parent_node, child_ids, new_nodes, rng):
#                     edges_table.add_row(left=left, right=pos, child=new_edge[0], parent=new_edge[1])
#                     # print("new_edge: left={}, right={}, child={}, parent={}".format(
#                     #     left, pos, new_edge[0], new_edge[1]))
#             else:
#                 # Previous node was not a polytomy - just add the edges_out, with modified left
#                 for edge_id in child_edge_ids:
#                     if existing_edges_right[edge_id] == pos:  # this edge has just gone out
#                         edges_table.add_row(
#                             left=existing_edges_left[edge_id],
#                             right=pos,
#                             parent=parent_node,
#                             child=existing_edges_child[edge_id],
#                         )

#         for edge in e_out: 
#             if edge.parent != tskit.NULL:
#                 # print("REMOVE", edge.id)
#                 edges_for_node[edge.parent].remove(edge.id)
#         for edge in e_in:
#             if edge.parent != tskit.NULL:
#                 # print("ADD", edge.id)
#                 edges_for_node[edge.parent].add(edge.id)            

#         # Chop if we have created a polytomy: the polytomy itself will be resolved
#         # at a future iteration, when any of the edges move in or out of the polytomy
#         while nodes_changed:
#             node = nodes_changed.pop()
#             edge_ids = edges_for_node[node]
#             # print("Looking at", node)

#             if len(edge_ids) == 0:
#                 del edges_for_node[node]
#             # if this node has changed *to* a polytomy, we need to cut all of the
#             # child edges that were previously present by adding the previous segment
#             # and left-truncating
#             elif len(edge_ids) >= 3:
#                 # print("Polytomy at", node, " breaking edges")
#                 for edge_id in edge_ids:
#                     if existing_edges_left[edge_id] < interval[0]:
#                         tables.edges.add_row(
#                             left=existing_edges_left[edge_id],
#                             right=interval[0],
#                             parent=existing_edges_parent[edge_id],
#                             child=existing_edges_child[edge_id],
#                         )
#                     existing_edges_left[edge_id] = interval[0]
#     assert len(edges_for_node) == 0

#     tables.edges.squash()
#     tables.sort() # Shouldn't need to do this: https://github.com/tskit-dev/tskit/issues/808

#     return tables.tree_sequence()

# import io
# import collections
# import tqdm
# import time

# import msprime
# import tsinfer

# ### Check equiprobable

# nodes_polytomy_4 = """\
# id      is_sample   population      time
# 0       1       0               0.00000000000000
# 1       1       0               0.00000000000000
# 2       1       0               0.00000000000000
# 3       1       0               0.00000000000000
# 4       0       0               1.00000000000000
# """
# edges_polytomy_4 = """\
# id      left            right           parent  child
# 0       0.00000000      1.00000000      4       0,1,2,3
# """

# poly_ts = tskit.load_text(
#     nodes=io.StringIO(nodes_polytomy_4),
#     edges=io.StringIO(edges_polytomy_4),
#     strict=False,
# )

# print(poly_ts)



# trees = collections.Counter()
# for seed in tqdm.trange(1, 100):
#     ts2 = resolve_polytomies(poly_ts, random_seed=seed)
#     print(ts2.first().as_newick())
#     trees.update([ts2.first().rank()])
# print(trees.as_newick())

# ### Time on 10000 tip inferred TS and check

# ts_old = msprime.simulate(10000, recombination_rate=100, mutation_rate=10, random_seed=123)
# sd = tsinfer.SampleData.from_tree_sequence(ts_old, use_times=False)
# ts = tsinfer.infer(sd)
# print(f"{ts.num_samples} tips ; {ts.num_trees} trees")
# start = time.time()
# ts2 = resolve_polytomies(ts, random_seed=1)
# print("Time (s):", time.time()-start)
# for tree in ts2.trees():
#     for node in tree.nodes():
#         assert tree.num_children(node) < 3

import dendropy
import random
from dendropy.calculate import treecompare
import matplotlib.pyplot as plt
import tqdm

tree = "(Z_0,(C_0:0, G_0:0, L_0:0, R_0:0, Q_0:0)"
tree = dendropy.Tree.get(data="(Z_0,((((C_0:1e-06,Q_0:1e-06):1e-06,R_0:0.0):0.0,L_0:0.0):0.0,G_0:0.0));", schema="newick", rooting="force-rooted")

non_info_gt_cnt = {}
tns = tree.taxon_namespace
rng=random.Random(0)
for i in tqdm.trange(607):
    # tree = dendropy.Tree.get(data="(Z_0,((((C_0:1e-06,Q_0:1e-06):1e-06,R_0:0.0):0.0,L_0:0.0):0.0,G_0:0.0));", schema="newick", rooting="force-rooted", taxon_namespace=tns)
    threshold=1e-5
    collapsed = False
    for e in tree.postorder_edge_iter():
        if e.length is not None and ((e.length <= threshold) and e.is_internal()):
            e.collapse()
            collapsed = True
            if len(tree.seed_node._child_nodes[0]._child_nodes) == 5 or len(tree.seed_node._child_nodes[1]._child_nodes) == 5:
                print(True)
    
    # tree.collapse_unweighted_edges(threshold=1e-5)
    if collapsed:
        tree.encode_bipartitions()
        # print(tree.as_string(schema='newick'))
        tree.resolve_polytomies(limit=2, update_bipartitions=True, rng=random)
        print(tree.as_string(schema='newick'))
        # print("------------")
    find = False
    for x in non_info_gt_cnt.keys():
        if treecompare.symmetric_difference(x, tree) == 0:
            non_info_gt_cnt[x] += 1
            find = True
            break
    if not find:
        non_info_gt_cnt[tree.clone(depth=1)] = 1


cnt_list = list(non_info_gt_cnt.values())
# print(non_info_gt_cnt)
cnt_list.sort(reverse=True)
x = len(cnt_list)
print(cnt_list[0], cnt_list[-1])

number_of_genes = 105
while x < number_of_genes:
    cnt_list.append(0)
    x += 1

total = sum(cnt_list)
cnt_list = [x/total for x in cnt_list]

fig, ax = plt.subplots(1, 1, figsize=(15, 5))
ax.bar(range(105), cnt_list)
ax.set_ylabel("Probability", fontsize=16)
ax.set_xlabel("Gene Tree", fontsize=16)
plt.savefig("gt_info_test.pdf")
plt.show()

# from ete3 import Tree
# non_info_gt_cnt = {}
# for i in tqdm.trange(10000):
#     tree = Tree("(Z_0,(C_0:1e-06,Q_0:1e-06,R_0:0.0,L_0:0.0,G_0:0.0));")
#     tree.resolve_polytomy(recursive=True)
#     print(str(tree.write()))
#     find = False
#     for x in non_info_gt_cnt.keys():
#         rf= x.robinson_foulds(tree)[0]
#         print(rf)
#         if rf < 0.0001:
#             non_info_gt_cnt[x] += 1
#             find = True
#             break
#     if not find:
#         non_info_gt_cnt[tree] = 1

# cnt_list = list(non_info_gt_cnt.values())
# # print(non_info_gt_cnt)
# cnt_list.sort(reverse=True)
# x = len(cnt_list)
# print(cnt_list[0], cnt_list[-1])

# number_of_genes = 105
# while x < number_of_genes:
#     cnt_list.append(0)
#     x += 1

# total = sum(cnt_list)
# cnt_list = [x/total for x in cnt_list]

# fig, ax = plt.subplots(1, 1, figsize=(15, 5))
# ax.bar(range(105), cnt_list)
# ax.set_ylabel("Probability", fontsize=16)
# ax.set_xlabel("Gene Tree", fontsize=16)
# plt.show()