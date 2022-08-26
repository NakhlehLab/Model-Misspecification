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
    
def generate_locus():
    

    RUN_SEQ_GEN = (sys.argv[1].lower() == "true")
    gtfile = sys.argv[2]
    num_loci = int(sys.argv[3])
    
        
    if RUN_SEQ_GEN: 
        LOCUS_HETEROGENEITY = (sys.argv[4].lower() == "true")
        sequence_path = sys.argv[5]
        mu_rate_path = sys.argv[6]
        SEQ_LENGTH = int(sys.argv[7])
        BASE_RATE = float(sys.argv[8])
        run_seq_gen(gtfile, sequence_path, mu_rate_path, LOCUS_HETEROGENEITY, num_loci, SEQ_LENGTH, BASE_RATE)

if __name__ == '__main__':
    generate_locus()