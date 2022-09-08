import os
import dendropy
from dendropy.calculate import treecompare
import sys


def runIQ(marker_path, dirname, outgroup, num_species, locus_length):
    locus_to_seq = {}
    iqtree_prefix = "/Users/zhen/Desktop/Zhen/research/phylogenetics/treeAugment/proj/iqtree-1.6.10-MacOSX/bin/iqtree -s "
    iqtree_suffix = " -o " + outgroup

    with open(marker_path, "r") as markers_file:
        locus_name = ''
        for line in markers_file:
            if '[' in line:
                locus_name = line[line.find('[') + 1: line.find(',')]
                locus_to_seq[locus_name] = {}
            else:
                arr = line.split()
                name = arr[0]
                seq = arr[1]

                locus_to_seq[locus_name][name] = seq

    for locus_name in locus_to_seq.keys():
        phylip_path = os.path.join(dirname, "phylip_%s.phy" % (locus_name))
        if os.path.exists(phylip_path):
            print(phylip_path + " already exists!")
        else:
            with open(phylip_path, "w") as phylip_file:
                phylip_file.write(str(num_species) + "\t" + str(locus_length) + "\n")
                for taxon_name in locus_to_seq[locus_name].keys():
                    phylip_file.write(taxon_name)
                    phylip_file.write(" ")
                    phylip_file.write(locus_to_seq[locus_name][taxon_name])
                    phylip_file.write("\n")
        print(iqtree_prefix + phylip_path + iqtree_suffix)
        os.system(iqtree_prefix + phylip_path + iqtree_suffix)


def writeGts(treepath, dirname, num_loci):
    with open(treepath, "w") as outf:
        result = ""
        for i in range(1, num_loci + 1):
            iqpath = os.path.join(dirname, "phylip_loci" + str(i) + ".phy.treefile")
            if os.path.exists(iqpath):
                with open(iqpath, "r") as inf:
                    result += inf.read()
        outf.write(result)
        outf.close()


def rerootTree(treepath, rooted_treepath, outgroup, PRUNE_OUTGROUP):
    treelist = []
    with open(treepath, "r") as handle:
        lines = handle.read().strip().split("\n")
        for line in lines:
            tree = dendropy.Tree.get(data=line.strip(), schema='newick', rooting="force-unrooted")
            print(tree.as_ascii_plot())
            outgroup_node = tree.find_node_with_taxon_label(outgroup)
            print(outgroup_node.edge)
            tree.reroot_at_edge(outgroup_node.edge, update_bipartitions=False)
            # tree.reroot_at_node(outgroup_node, update_bipartitions=False)
            # tree.to_outgroup_position(outgroup_node, update_bipartitions=False)
            if PRUNE_OUTGROUP:
                tree.prune_taxa_with_labels([outgroup])
            print("After:")
            print(tree.as_string(schema='newick'))
            print(tree.as_ascii_plot())
            treelist.append(tree.as_string(schema='newick')[5:])

    with open(rooted_treepath, "w") as outf:
        outf.write("".join(treelist))


def run_iqtree():
    print(sys.argv)
    folder = sys.argv[1]
    
    marker_path = folder+"/markers.txt"
    dirname = folder + "/IQTREE"
    treepath = folder + "/iqtree.txt"
    rooted_treepath = folder + "/rooted_iqtree.txt"
    PRUNE_OUTGROUP = (sys.argv[2].lower() == "true")
    # HETER = int(sys.argv[3])
    RUN_IQ = (sys.argv[3].lower() == "true")
    outgroup = sys.argv[4]
    num_species = sys.argv[5]
    locus_length = sys.argv[6]
    biological = (sys.argv[7].lower() == "true")
    num_loci = sys.argv[8]
    print(biological)
    if RUN_IQ:
        if not biological:
            runIQ(marker_path, dirname, outgroup + "_0", num_species, locus_length)
        else:
            runIQ(marker_path, dirname, outgroup, num_species, locus_length)
        writeGts(treepath, dirname, num_loci)
    if not biological:
        rerootTree(treepath, rooted_treepath, outgroup + " 0", PRUNE_OUTGROUP)
    else:
        rerootTree(treepath, rooted_treepath, outgroup, PRUNE_OUTGROUP)


# def reroot_true_genetrees():
#     treepath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/data/simulation/seqgen/" + folder + "/genetrees.txt"
#     rooted_treepath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/data/simulation/seqgen/" + folder + "/rooted_genetrees.txt"
#     PRUNE_OUTGROUP = bool(sys.argv[3])
#     rerootTree(treepath, rooted_treepath, "G 0", PRUNE_OUTGROUP)


if __name__ == '__main__':
    # run_iqtree()
    # reroot_true_genetrees()
    # compare_iqtree_truetree()

    runIQ("/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/long/3/heter/500/markers6000.txt", "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/long/3/heter/500/IQTREE/", "Z_0", 6, 500)
