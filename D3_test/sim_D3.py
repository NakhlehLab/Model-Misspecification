
#Author: Huw
import csv
import subprocess
import dendropy
import numpy
import scipy.stats

def calc_d3_stat(dist_matrix, outgroup_index):
	if outgroup_index == 0:
		return (dist_matrix[0,1] - dist_matrix[0,2])/(dist_matrix[0,1] + dist_matrix[0,2])
	elif outgroup_index == 1:
		return (dist_matrix[0,1] - dist_matrix[1,2])/(dist_matrix[0,1] + dist_matrix[1,2])
	else:
		return (dist_matrix[0,2] - dist_matrix[1,2])/(dist_matrix[0,2] + dist_matrix[1,2])

# use z-score of 0 to calculate p-value
def calc_p_from_bs_reps(bs_reps):
	bs_reps_mean = numpy.mean(bs_reps)
	bs_reps_std = numpy.std(bs_reps)

	abs_z_score = abs(bs_reps_mean / bs_reps_std)

	# two-tailed test
	p_value = 2 * (1 - scipy.stats.norm.cdf(abs_z_score))

	return(p_value)

pop_size = 0.05

n_genes = 100
g_length = 500
m_length = n_genes * g_length

n_blocks = 100 # also the number of bootstraps per rep
b_length = m_length // n_blocks

n_reps = 100

sg = dendropy.interop.seqgen.SeqGen()
sg.seq_len = g_length

og_matrix = numpy.zeros((3, m_length), dtype = numpy.ubyte)
bs_matrix = numpy.zeros((3, m_length), dtype = numpy.ubyte)

s_taxa = dendropy.TaxonNamespace(["a", "b", "c"])
s_tree = dendropy.Tree.get(data = "((a:0.05,b:0.05):0.05,c:0.1);", schema = "newick", taxon_namespace = s_taxa)

taxon_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(containing_taxon_namespace = s_taxa, num_contained = 1)
reverse_taxon_map = {} # only works in the simple case of one specimen per species

for g_taxon in taxon_map:
	s_taxon = taxon_map[g_taxon]
	reverse_taxon_map[s_taxon] = g_taxon

g_taxa = [reverse_taxon_map[s_taxon] for s_taxon in s_taxa]

output_file = open("d3_reps.csv", "w")
output_writer = csv.writer(output_file)
output_writer.writerow(["p.dist.p.value", "jc.dist.p.value"])

all_p_dist_p_values = numpy.zeros(n_reps)
all_jc_dist_p_values = numpy.zeros(n_reps)

for rep_i in range(n_reps):
	# construct original sequence matrix
	for g_tree_i in range(n_genes):
		gene_tree = dendropy.simulate.treesim.contained_coalescent_tree(containing_tree = s_tree, gene_to_containing_taxon_map = taxon_map, default_pop_size = pop_size)

		for root_child in gene_tree.seed_node.child_node_iter():
			if root_child.is_leaf():
				outgroup_index = g_taxa.index(root_child.taxon)
			else:
				t1 = root_child.age

		t2 = gene_tree.seed_node.age

		sg_output = sg.generate(gene_tree)

		for i in range(3):
			taxon = g_taxa[i]
			sequence = sg_output.char_matrices[0][taxon]
			sequence_bytes = str(sequence).encode("utf8")
			og_matrix[i,g_tree_i*g_length:(g_tree_i + 1)*g_length] = numpy.frombuffer(sequence_bytes, dtype = numpy.ubyte)

	jc_dist_matrix = numpy.zeros((3, 3))
	p_dist_matrix = numpy.zeros((3, 3))

	bs_p_d3 = numpy.zeros(n_blocks)
	bs_jc_d3 = numpy.zeros(n_blocks)

	for bs_rep_i in range(n_blocks):
		# construct bootstrap sequence matrix
		for bs_block_i in range(n_blocks):
			og_block_i = numpy.random.randint(n_blocks)
			bs_matrix[0:3,bs_block_i*b_length:(bs_block_i + 1)*b_length] = og_matrix[0:3,og_block_i*b_length:(og_block_i + 1)*b_length]

		for i in range(2):
			for j in range(i + 1, 3):
				p_dist = (m_length - numpy.sum(bs_matrix[i] == bs_matrix[j])) / m_length
				jc_dist = -(3/4) * numpy.log(1 - (4 / 3) * p_dist)

				p_dist_matrix[i,j] = p_dist
				jc_dist_matrix[i,j] = jc_dist

		bs_p_d3[bs_rep_i] = calc_d3_stat(p_dist_matrix, 2)
		bs_jc_d3[bs_rep_i] = calc_d3_stat(jc_dist_matrix, 2)

	p_dist_p_value = calc_p_from_bs_reps(bs_p_d3)
	jc_dist_p_value = calc_p_from_bs_reps(bs_jc_d3)

	all_p_dist_p_values[rep_i] = p_dist_p_value
	all_jc_dist_p_values[rep_i] = jc_dist_p_value

	output_writer.writerow([p_dist_p_value, jc_dist_p_value])

output_file.close()

p_dist_p_hist = numpy.histogram(all_p_dist_p_values, bins = 20, range = (0, 1))
jc_dist_p_hist = numpy.histogram(all_jc_dist_p_values, bins = 20, range = (0, 1))

print(p_dist_p_hist)
print(jc_dist_p_hist)
