import os
from itertools import combinations
import random
import copy
import pathlib
import sys



def marker2mcmcnex(marker_path, nex_path, num_taxa, locus_length, taxon_map, sd, murate, varyps):
    s = "#NEXUS\n  Begin data;\n  Dimensions ntax="+num_taxa+" nchar="+str(int(num_taxa)*int(locus_length))+";\n  Format datatype=dna symbols=\"ACTG\" missing=? gap=-;\n  Matrix\n"
    # s_end = ";End;\nBEGIN PHYLONET;\nMCMC_SEQ -cl 50000000 -bl 5000000 -sf 5000 -sd "+ sd +" -pl 16 -mc3 (2.0) -tm "+taxon_map
    # s_end = ";End;\nBEGIN PHYLONET;\nMCMC_SEQ -cl 80000000 -bl 8000000 -sf 5000 -sd "+ sd +" -pl 16 -tm "+taxon_map
    s_end = ";End;\nBEGIN PHYLONET;\nMCMC_SEQ -cl 50000000 -bl 5000000 -sf 5000 -sps 0.01 -sd "+ sd +" -pl 24 -tm "+taxon_map
    if varyps:
        s_end += " -varyps"
    if murate:
        s_end += " -murate; \nEND;"
    else:
        s_end += "; \nEND;"
    print("what the fuck!")
    print(marker_path, nex_path)

    with open(marker_path, "r") as handle:
        content = handle.read()
        s += content

    print("done read")
    s += s_end
    print(s)
    with open(nex_path, "w") as handle:
        handle.write(s)

def genetree2MLnex(gt_path, nex_path, num_reti, taxon_map, ):
    s = "#NEXUS\n\nBEGIN TREES;\n"
    s_end = "End;\n\nBEGIN PHYLONET;\nInfernetwork_ML (all) "+str(num_reti)+" -pl 4 -po -di -x 50 -n 10 -a "+taxon_map+"; \nEND;"

    with open(gt_path, "r") as handle:
        lines = handle.read().strip().split("\n")
        i = 1
        for line in lines:
            s += "Tree genetree" + str(i) + "=" + line + "\n"
            i += 1

    s += s_end
    print(s)
    with open(nex_path, "w") as handle:
        handle.write(s)

def genetree2MLnex_fixtopo(gt_path, nex_path, num_reti, taxon_map, network_newick):
    s = "#NEXUS\n\nBEGIN TREES;\n"
    s_end = "End;\n\nBEGIN PHYLONET;\nInfernetwork_ML (all) "+str(num_reti)+" -pl 4 -po -di -x 0 -n 1 -a "+taxon_map+" -s net;\nEND;"

    with open(gt_path, "r") as handle:
        lines = handle.read().strip().split("\n")
        i = 1
        for line in lines:
            s += "Tree genetree" + str(i) + "=" + line + "\n"
            i += 1
    s += "BEGIN NETWORKS;\n\n Network net="+network_newick+"\n\nEND;\n"
    s += s_end
    print(s)
    with open(nex_path, "w") as handle:
        handle.write(s)


def prepare_beast(path, num_taxa, locus_length, taxaset=["G_0","A_0","Q_0","L_0"]):
    # taxaset = ["G_0","A_0","Q_0","L_0"]
    prefix = "begin taxa;\n dimensions ntax="+str(num_taxa)+";\ntaxlabels\n"
    for species in taxaset:
        prefix += species
        prefix += "\n"
    prefix += ";\nend;\n"
    prefix += "begin characters;\ndimensions nchar="+str(locus_length)+";\n"
    prefix2 = "format missing=? gap=- matchchar=. datatype=dna;\noptions gapmode=missing;\nmatrix\n"

    # path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/data/simulation/seqgen/tree/heter_locus/1/markers.txt"
    with open(path, 'r') as file:
        content = file.read()
        lines = content.split("\n")
        i = 0
        loci2seqs = {}
        loci2length = {}
        loci2str = {}
        begin = False
        seqs = None
        while i < len(lines):
            line = lines[i]
            if line.startswith('['):
                begin = True

                lociinfo = line.split(",")
                locus_name = lociinfo[0][1:]
                locus_length = lociinfo[1][:-1]
                loci2length[locus_name] = locus_length
                loci2str[locus_name] = "[" + locus_name + "," + str(locus_length) + ",...]\n"
                if seqs is not None and len(seqs) == len(taxaset):
                    loci2seqs[locus_name] = seqs
                seqs = {}

            elif begin:
                if len(line) == 0:
                    break
                arr = line.split(" ")
                # print(arr)
                species, seq = arr[0], arr[1]
                if species in taxaset:
                    seqs[species] = seq
                    loci2str[locus_name] += species + " " + seq[:int(locus_length)] + "\n"

            i += 1
            file.close()
    print("finished reading")
    folder = pathlib.PurePath(path).parent
    output_folder = os.path.join(folder,"nexus")
    print(output_folder)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    for locus_name in loci2str.keys():
        # print(locus_name)
        print(loci2str[locus_name])

        output_file = os.path.join(output_folder, locus_name + ".nex")
        with open(output_file, "w") as handle:
            # prefix += str(loci2length[locus_name])+";\n"
            handle.write(prefix + prefix2 + "\n" + loci2str[locus_name] + "\n;end;\n")


if __name__ == '__main__':
    print(sys.argv)
    func_type = int(sys.argv[1])
    num_taxa = sys.argv[2]
    locus_length = sys.argv[3]
    if func_type == 1:
        marker2mcmcnex(sys.argv[4], sys.argv[5], num_taxa, locus_length, sys.argv[6], sys.argv[7], sys.argv[8].lower() == "true", sys.argv[9].lower() == "true")
    elif func_type == 2:
        print(sys.argv)
        genetree2MLnex(sys.argv[4], sys.argv[5], int(sys.argv[6]), sys.argv[7])
    elif func_type == 3:
        genetree2MLnex_fixtopo(sys.argv[4], sys.argv[5], int(sys.argv[6]), sys.argv[7])
    # else:
    #     prepare_beast(sys.argv[4], num_taxa, locus_length)