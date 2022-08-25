from cmath import sqrt
import numpy as np
import dendropy
import scipy
import numpy as np
import subprocess
import dendropy
from dendropy.calculate import treecompare
import sys
from pathlib import Path
import math

seqgen = "/Users/zhen/Desktop/Zhen/research/phylogenetics/treeAugment/proj/Seq-Gen/source/seq-gen"
ms = "/Users/zhen/Desktop/Zhen/research/phylogenetics/treeAugment/proj/msdir/ms"
# num_loci = 100
# SEQ_LENGTH = 200



# R/C 
# def RC_tree(ststr):
#     tree = dendropy.Tree.get(data=ststr, schema="newick", rooting="force-rooted")

# def generate_hetero(alpha, sigma=1):
#     shape = np.random.lognormal(alpha, sigma, 1)
#     print(shape)
#     x = np.random.gamma(shape, 1, 100)
#     print(x)

# # species specific heterogeneity
# def species_hetero(gtstr):
#      tree = dendropy.Tree.get(data=gtstr, schema="newick", rooting="force-rooted")


def gen_heter_gt(gtfile, true_gene_path, num_loci, NUM_SPECIES):
    treelist = []
    with open(gtfile, "r") as file:
        lines = file.read().split("\n")
        for i in range(num_loci):
            treestring = lines[i].strip()
            tree = UCLN_corrected(treestring, NUM_SPECIES)
            treestring = tree.as_string(schema="newick",)[5:]
            # print(treestring)
            treelist.append(treestring.strip())
            
    with open(true_gene_path, "w") as handle:
        handle.write("\n".join(treelist))


def run_seq_gen(gtfile, sequence_path, murate_path, LOCUS_HETEROGENEITY, num_loci, SEQ_LENGTH, base_locus_rate):
    alpha = [1 for i in range(num_loci)]
    if LOCUS_HETEROGENEITY:
        samples = np.random.dirichlet(alpha)
        # print(samples)
        samples = [x * num_loci for x in samples]
        print(samples)
        with open(murate_path, "w") as handle:
            sample_ss = [str(x) for x in samples]
            handle.write("\n".join(sample_ss))
            handle.close()
    else:
        samples = np.ones(num_loci)


    treelist = []
    with open(gtfile, "r") as file:
        lines = file.read().split("\n")
        locus2seq = {}
        for i in range(num_loci):
            treestring = lines[i].strip()
            treelist.append(treestring.strip())
            with open("./seqgen.tree", "w") as treefile:
                treefile.write(treestring)
            s = " -mGTR -s" + str(base_locus_rate * samples[i]) + " -f0.2112,0.2888,0.2896,0.2104 -r0.2173,0.9798,0.2575,0.1038,1.0,0.207 -l"+str(SEQ_LENGTH)+" -or ./seqgen.tree"
            command = seqgen + s
            result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # print(result)
            seqstring = result.stdout.decode("utf-8").strip().split("\n")[1:]
            name2seq = {}
            for ss in seqstring:
                arr = ss.split(" ")
                name2seq[arr[0]] = arr[1]
            # print(seqstring)
            locus2seq[i] = name2seq
    # print(locus2seq)
    

    with open(sequence_path, "w") as handle:
        line = ""
        for i in locus2seq.keys():
            line += "[loci"+str(i+1)+","+str(SEQ_LENGTH)+"]\n"
            for species in locus2seq[i].keys():
                line += species +" " + locus2seq[i][species] + "\n"
        handle.write(line)
    

# def run_seq_gen(gtfile, true_gene_path, sequence_path, murate_path, LINEAGE_HETEROGENEITY, LOCUS_HETEROGENEITY, num_loci, NUM_SPECIES, SEQ_LENGTH, base_locus_rate):
#     treelist = []
#     alpha = [1 for i in range(num_loci)]
#     print(alpha)
#     if LOCUS_HETEROGENEITY:
#         samples = np.random.dirichlet(alpha)
#         print(samples)
#         samples = [x * num_loci for x in samples]
#         print(samples)
#         with open(murate_path, "w") as handle:
#             sample_ss = [str(x) for x in samples]
#             handle.write("\n".join(sample_ss))
#             handle.close()
#     else:
#         samples = np.ones(num_loci)

#     with open(gtfile, "r") as file:

#         lines = file.read().split("\n")
#         locus2seq = {}
#         for i in range(num_loci):
#             treestring = lines[i].strip()
#             if LINEAGE_HETEROGENEITY:
#                 tree = UCLN_corrected(treestring, NUM_SPECIES)
#                 treestring = tree.as_string(schema="newick",)[5:]
#             print(treestring)
#             treelist.append(treestring.strip())
#             with open("./seqgen.tree", "w") as treefile:
#                 treefile.write(treestring)
#             s = " -mGTR -s" + str(base_locus_rate * samples[i]) + " -f0.2112,0.2888,0.2896,0.2104 -r0.2173,0.9798,0.2575,0.1038,1.0,0.207 -l"+str(SEQ_LENGTH)+" -or ./seqgen.tree"

#             command = seqgen + s
#             result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#             print(result)
#             seqstring = result.stdout.decode("utf-8").strip().split("\n")[1:]
#             name2seq = {}
#             for ss in seqstring:
#                 arr = ss.split(" ")
#                 name2seq[arr[0]] = arr[1]
#             print(seqstring)
#             locus2seq[i] = name2seq
#     print(locus2seq)
#     with open(true_gene_path, "w") as handle:
#         handle.write("\n".join(treelist))

#     with open(sequence_path, "w") as handle:
#         line = ""
#         for i in locus2seq.keys():
#             line += "[loci"+str(i+1)+","+str(SEQ_LENGTH)+"]\n"
#             for species in locus2seq[i].keys():
#                 line += species +" " + locus2seq[i][species] + "\n"
#         handle.write(line)

# gene lineage specific heterogeneity
def UCLN(gtstr, NUM_SPECIES):
    tree = dendropy.Tree.get(data=gtstr, schema="newick", rooting="force-rooted")
    samples = np.random.lognormal(0.1, 0.5, 2*NUM_SPECIES-1)
    # samples = np.random.lognormal(0.1, 1, 2*NUM_SPECIES-1)
    for idx, nd in enumerate(tree):
        if nd.edge.length is not None:
            nd.edge.length *= samples[idx]
    return tree

def UCLN_corrected(gtstr, NUM_SPECIES):
    tree = dendropy.Tree.get(data=gtstr, schema="newick", rooting="force-rooted")
    mu_x = 1
    sig_x = 0.5
    mu = np.log(mu_x**2/math.sqrt(mu_x**2+sig_x))
    sig = np.log(1+sig_x/(mu_x**2))
    samples = np.random.lognormal(mu, sig, 2*NUM_SPECIES-1)
    # samples = np.random.lognormal(0.1, 1, 2*NUM_SPECIES-1)
    for idx, nd in enumerate(tree):
        if nd.edge.length is not None:
            nd.edge.length *= samples[idx]
    return tree

def generate_locus():
    # sequence_path = sys.argv[3]
    # mu_rate_path = sys.argv[4]

    LINEAGE_HETEROGENEITY = (sys.argv[1].lower() == "true")
    RUN_SEQ_GEN = (sys.argv[2].lower() == "true")
    gtfile = sys.argv[3]
    num_loci = int(sys.argv[4])


    if LINEAGE_HETEROGENEITY:
        true_gene_path = sys.argv[5]
        NUM_SPECIES = int(sys.argv[6])
        gen_heter_gt(gtfile, true_gene_path, num_loci, NUM_SPECIES)
        
    elif RUN_SEQ_GEN: 
        LOCUS_HETEROGENEITY = (sys.argv[5].lower() == "true")
        sequence_path = sys.argv[6]
        mu_rate_path = sys.argv[7]
        SEQ_LENGTH = int(sys.argv[8])
        BASE_RATE = float(sys.argv[9])
        run_seq_gen(gtfile, sequence_path, mu_rate_path, LOCUS_HETEROGENEITY, num_loci, SEQ_LENGTH, BASE_RATE)

if __name__ == '__main__':
    generate_locus()