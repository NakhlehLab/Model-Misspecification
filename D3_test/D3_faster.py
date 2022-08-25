import dendropy
import numpy as np

import sys
from itertools import combinations
import pandas as pd
from pathlib import Path
import numpy as np
from scipy.stats import norm
import time
import numpy
import dendropy
# from dendropy.simulate import treesim
# from dendropy.interop import seqgen

# taxa = dendropy.TaxonNamespace(["a", "b", "c"])
SPECIES = ["C_0", "G_0", "L_0", "Q_0", "R_0", "Z_0"]
OUTGROUP = ["Z_0"]
INGROUP_SPECIES = list(set(SPECIES) - set(OUTGROUP))
# TRIPLET_LIST = list(combinations(INGROUP_SPECIES, 3))
print(INGROUP_SPECIES)

# n_taxa = 3
N_loci = 10000
# locus_length = 2000
N_bootstrap = 1000

# s = seqgen.SeqGen()
# s.seq_len = locus_length

# original = numpy.zeros((n_loci, n_taxa, locus_length), dtype = "uint8")
# bootstrap = numpy.zeros((n_loci, locus_length, n_taxa), dtype = "uint8")

# for locus_i in range(n_loci):
# 	random_tree = treesim.pure_kingman_tree(taxon_namespace = taxa, pop_size = 0.1)
# 	seqgen_output = s.generate(random_tree)
# 	alignment = seqgen_output.char_matrices[0]

# 	for taxon_i, taxon in enumerate(taxa):
# 		taxon_sequence = alignment[taxon].symbols_as_string()
# 		original[locus_i, taxon_i] = numpy.frombuffer(taxon_sequence.encode("utf8"), dtype = "uint8")

# start_time = time.time()

def calculate_D3(D_A_B, D_A_C, D_B_C):
    three_pairwise = [(D_A_B - D_A_C)/float(D_A_B + D_A_C + 0.0000001), (D_A_B - D_B_C)/float(D_A_B + D_B_C + 0.0000001), (D_A_C - D_B_C)/float(D_A_C + D_B_C + 0.0000001)]
    return min(three_pairwise, key=abs)

def calc_pval_from_bootstrap(D3_estimates):
    bs_reps_mean = np.mean(D3_estimates)
    bs_reps_std = np.std(D3_estimates)
    abs_z_score = abs(bs_reps_mean / bs_reps_std)
    # two-tailed test
    p_value = 2 * norm.sf(abs_z_score)  
    return bs_reps_mean, bs_reps_std, p_value


def compute_D3_concat(marker_path, locus_length):
    with open(marker_path, "r") as handle:
        D3_list = {"A":[], "B":[], "C":[], "D3": [], "D3_mean": [], "D3_stdev":[], "D3_pval": []}
        seqs_list = {}
        INGROUP_SPECIES.sort()
        for x in INGROUP_SPECIES:
            seqs_list[x] = []
        lines = handle.readlines()
        
        for line in lines:
            arr = line.split(" ")
            if arr[0] in INGROUP_SPECIES:
                seqs_list[arr[0]].append(arr[1].strip())

        N_taxa = len(seqs_list)
        seqs_array = []
        for taxon in INGROUP_SPECIES: 
            seqs_array.append(seqs_list[taxon])
        
        
        original = np.zeros((N_loci, N_taxa, locus_length), dtype = "uint8")
        bootstrap = np.zeros((N_loci, locus_length, N_taxa), dtype = "uint8")

        for locus_i in range(N_loci):
            for taxon_i, taxon in enumerate(INGROUP_SPECIES):
                taxon_sequence = seqs_array[taxon_i][locus_i]
                original[locus_i, taxon_i] = np.frombuffer(taxon_sequence.encode("utf8"), dtype = "uint8")
        
        D3_estimates = np.zeros((N_taxa, N_taxa, N_taxa, N_bootstrap), dtype="float64")
        
        for bootstrap_i in range(N_bootstrap):
            locus_indices = numpy.random.randint(N_loci, size = N_loci)

            for bootstrap_locus_i in range(N_loci):
                site_indices = numpy.random.randint(locus_length, size = locus_length)
                original_locus_i = locus_indices[bootstrap_locus_i]
                bootstrap[bootstrap_locus_i] = original[original_locus_i].T[site_indices]

            concatenated_bootstrap = bootstrap.reshape((N_loci * locus_length, N_taxa)).T
            distance_array = np.zeros((N_taxa, N_taxa), dtype = "float64")

            for x in range(N_taxa-1):
                for y in range(x+1, N_taxa):
                    distance_array[x][y] = np.mean(concatenated_bootstrap[x] != concatenated_bootstrap[y])
                    
            for x in range(N_taxa-2):
                for y in range(x+1, N_taxa-1):
                    for z in range(y+1, N_taxa):
                        D3_estimates[x][y][z][bootstrap_i] = calculate_D3(distance_array[x][y], distance_array[x][z], distance_array[y][z])
            
            
        distance_array_original = np.zeros((N_taxa, N_taxa), dtype = "float64")
        concatenated_original = original.transpose(0,2,1).reshape((N_loci * locus_length, N_taxa)).T
        for x in range(N_taxa-1):
            for y in range(x+1, N_taxa):
                distance_array_original[x][y] = np.mean(concatenated_original[x] != concatenated_original[y])

        for x in range(N_taxa-2):
            for y in range(x+1, N_taxa-1):
                for z in range(y+1, N_taxa):
                    D3_mean, D3_stdev, D3_pval = calc_pval_from_bootstrap(D3_estimates[x][y][z])
                    D3_val = calculate_D3(distance_array_original[x][y], distance_array_original[x][z], distance_array_original[y][z])
                    D3_list["A"].append(INGROUP_SPECIES[x])
                    D3_list["B"].append(INGROUP_SPECIES[y])
                    D3_list["C"].append(INGROUP_SPECIES[z])
                    D3_list["D3"].append(D3_val)
                    D3_list["D3_mean"].append(D3_mean)
                    D3_list["D3_stdev"].append(D3_stdev)
                    D3_list["D3_pval"].append(D3_pval)
            
        df = pd.DataFrame(D3_list)
        df.to_csv(str(Path(marker_path).parent.absolute())+"/D3_2stage_faster.txt")
    return df

def compute_D3_concat_faster(marker_path, locus_length):
    with open(marker_path, "r") as handle:
        D3_list = {"A":[], "B":[], "C":[], "D3": [], "D3_mean": [], "D3_stdev":[], "D3_pval": []}
        seqs_list = {}
        INGROUP_SPECIES.sort()
        for x in INGROUP_SPECIES:
            seqs_list[x] = []
        lines = handle.readlines()
        
        for line in lines:
            arr = line.split(" ")
            if arr[0] in INGROUP_SPECIES:
                seqs_list[arr[0]].append(arr[1].strip())

        N_taxa = len(seqs_list)
        seqs_array = []
        for taxon in INGROUP_SPECIES: 
            seqs_array.append(seqs_list[taxon])
        
        
        original = np.zeros((N_loci, N_taxa, locus_length), dtype = "uint8")
        bootstrap = np.zeros((N_loci, locus_length, N_taxa), dtype = "uint8")

        for locus_i in range(N_loci):
            for taxon_i, taxon in enumerate(INGROUP_SPECIES):
                taxon_sequence = seqs_array[taxon_i][locus_i]
                original[locus_i, taxon_i] = np.frombuffer(taxon_sequence.encode("utf8"), dtype = "uint8")
        
        D3_estimates = np.zeros((N_taxa, N_taxa, N_taxa, N_bootstrap), dtype="float64")
        
        for bootstrap_i in range(N_bootstrap):
            locus_indices = numpy.random.randint(N_loci, size = N_loci)

            for bootstrap_locus_i in range(N_loci):
                # site_indices = numpy.random.randint(locus_length, size = locus_length)
                original_locus_i = locus_indices[bootstrap_locus_i]
                bootstrap[bootstrap_locus_i] = original[original_locus_i].T

            concatenated_bootstrap = bootstrap.reshape((N_loci * locus_length, N_taxa)).T
            distance_array = np.zeros((N_taxa, N_taxa), dtype = "float64")

            for x in range(N_taxa-1):
                for y in range(x+1, N_taxa):
                    distance_array[x][y] = np.mean(concatenated_bootstrap[x] != concatenated_bootstrap[y])
                    
            for x in range(N_taxa-2):
                for y in range(x+1, N_taxa-1):
                    for z in range(y+1, N_taxa):
                        D3_estimates[x][y][z][bootstrap_i] = calculate_D3(distance_array[x][y], distance_array[x][z], distance_array[y][z])
            
            
        distance_array_original = np.zeros((N_taxa, N_taxa), dtype = "float64")
        concatenated_original = original.transpose(0,2,1).reshape((N_loci * locus_length, N_taxa)).T
        for x in range(N_taxa-1):
            for y in range(x+1, N_taxa):
                distance_array_original[x][y] = np.mean(concatenated_original[x] != concatenated_original[y])

        for x in range(N_taxa-2):
            for y in range(x+1, N_taxa-1):
                for z in range(y+1, N_taxa):
                    D3_mean, D3_stdev, D3_pval = calc_pval_from_bootstrap(D3_estimates[x][y][z])
                    D3_val = calculate_D3(distance_array_original[x][y], distance_array_original[x][z], distance_array_original[y][z])
                    D3_list["A"].append(INGROUP_SPECIES[x])
                    D3_list["B"].append(INGROUP_SPECIES[y])
                    D3_list["C"].append(INGROUP_SPECIES[z])
                    D3_list["D3"].append(D3_val)
                    D3_list["D3_mean"].append(D3_mean)
                    D3_list["D3_stdev"].append(D3_stdev)
                    D3_list["D3_pval"].append(D3_pval)
            
        df = pd.DataFrame(D3_list)
        df.to_csv(str(Path(marker_path).parent.absolute())+"/D3.txt")
    return df


if __name__ == "__main__":
    # marker_path= "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/long/93/markers.txt"
    # locus_length = 500

    # compute_D3_concat_faster(marker_path, locus_length)
    compute_D3_concat_faster(sys.argv[1], int(sys.argv[2]))
    # compute_D3_concat(sys.argv[1], int(sys.argv[2]))

    