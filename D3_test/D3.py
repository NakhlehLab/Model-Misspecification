import sys
from itertools import combinations
import pandas as pd
from pathlib import Path
import numpy as np
from scipy.stats import norm


SPECIES = ["C_0", "G_0", "L_0", "Q_0", "R_0", "Z_0"]
OUTGROUP = ["Z_0"]
INGROUP_SPECIES = list(set(SPECIES) - set(OUTGROUP))
BLOCK_SIZE=500
BS_NUM_REPS = 1000

# markers
def calculate_D3(A, B, C):
    D_A_B, D_A_C, D_B_C = 0, 0, 0
    for i in range(len(A)): #Count # pairwise differences
        if A[i] != B[i]:
            D_A_B += 1
        if A[i] != C[i]:
            D_A_C += 1
        if B[i] != C[i]:
            D_B_C += 1

    #Per-site divergence  
    D_A_B, D_A_C, D_B_C = float(D_A_B)/len(A), float(D_A_C)/len(A), float(D_B_C)/len(A)

    #Three pairwise divergence comparisons 
    three_pairwise = [(D_A_B - D_A_C)/float(D_A_B + D_A_C + 0.0000001), (D_A_B - D_B_C)/float(D_A_B + D_B_C + 0.0000001), (D_A_C - D_B_C)/float(D_A_C + D_B_C + 0.0000001)]

    # three_pairwise_abs = [abs(value) for value in three_pairwise] #Get absolute values of comparisons

    # index_min = min(range(len(three_pairwise_abs)), key=three_pairwise_abs.__getitem__) #get index of D3 value

    # D3 = three_pairwise[index_min] #D3 is the difference of smallest magnitude 
    D3 = min(three_pairwise, key=abs)

    return D3


def bootstrap_D3_overlapping(A, B, C, n_replicates):
    D3_estimates = [] #list for distribution of D3 statistics
    for i in range(n_replicates):

        A_bootstrap = []
        B_bootstrap = []
        C_bootstrap = []

        #sample bootstrap windows
        for j in range(100):
            bootstrap_index = np.random.randint(len(A)-10000) #index to sample alignment
            A_bootstrap.append(A[bootstrap_index:bootstrap_index+10000])
            B_bootstrap.append(B[bootstrap_index:bootstrap_index+10000])
            C_bootstrap.append(C[bootstrap_index:bootstrap_index+10000])
            
        #concatenate windows
        A_bootstrap = ''.join(A_bootstrap)
        B_bootstrap = ''.join(B_bootstrap)
        C_bootstrap = ''.join(C_bootstrap)

        #estimate D3
        D3_estimates.append(calculate_D3(A_bootstrap, B_bootstrap, C_bootstrap))

    mean_D3 = sum(D3_estimates)/float(len(D3_estimates))	

    #estimate variance
    D3_stdev = ((sum([(x - mean_D3)**2 for x in D3_estimates]))/len(D3_estimates))**(0.5)
    return D3_stdev


def calc_pval_from_bootstrap(D3_estimates):
    bs_reps_mean = np.mean(D3_estimates)
    bs_reps_std = np.std(D3_estimates)
    abs_z_score = abs(bs_reps_mean / bs_reps_std)
    # two-tailed test
    p_value = 2 * norm.sf(abs_z_score)  
    return bs_reps_mean, bs_reps_std, p_value




def bootstrap_D3(A, B, C, n_replicates):
    D3_estimates = [] #list for distribution of D3 statistics
    num_blocks = len(A) // BLOCK_SIZE

    for i in range(n_replicates):
        A_bootstrap_list = []
        B_bootstrap_list = []
        C_bootstrap_list = []

        #sample bootstrap windows
        for j in range(num_blocks):
            bootstrap_index = np.random.randint(num_blocks) #index to sample alignment
            A_bootstrap_list.append(A[bootstrap_index * BLOCK_SIZE : (bootstrap_index + 1) * BLOCK_SIZE])
            B_bootstrap_list.append(B[bootstrap_index * BLOCK_SIZE : (bootstrap_index + 1) * BLOCK_SIZE])
            C_bootstrap_list.append(C[bootstrap_index * BLOCK_SIZE : (bootstrap_index + 1) * BLOCK_SIZE])
            
        #concatenate windows
        A_bootstrap = ''.join(A_bootstrap_list)
        B_bootstrap = ''.join(B_bootstrap_list)
        C_bootstrap = ''.join(C_bootstrap_list)
        
        D3_val = calculate_D3(A_bootstrap, B_bootstrap, C_bootstrap)
        D3_estimates.append(D3_val)
    
    return calc_pval_from_bootstrap(D3_estimates)

def bootstrap_1_D3(A, B, C, n_replicates):
    D3_estimates = [] #list for distribution of D3 statistics
    num_blocks = len(A) // BLOCK_SIZE

    for i in range(n_replicates):
        A_bootstrap_list = []
        B_bootstrap_list = []
        C_bootstrap_list = []

        #sample bootstrap windows
        for j in range(num_blocks-1):
            bootstrap_index = np.random.randint(num_blocks) #index to sample alignment
            A_bootstrap_list.append(A[bootstrap_index * BLOCK_SIZE : (bootstrap_index + 1) * BLOCK_SIZE])
            B_bootstrap_list.append(B[bootstrap_index * BLOCK_SIZE : (bootstrap_index + 1) * BLOCK_SIZE])
            C_bootstrap_list.append(C[bootstrap_index * BLOCK_SIZE : (bootstrap_index + 1) * BLOCK_SIZE])
            
        #concatenate windows
        A_bootstrap = ''.join(A_bootstrap_list)
        B_bootstrap = ''.join(B_bootstrap_list)
        C_bootstrap = ''.join(C_bootstrap_list)
        
        D3_val = calculate_D3(A_bootstrap, B_bootstrap, C_bootstrap)
        D3_estimates.append(D3_val)
    
    return calc_pval_from_bootstrap(D3_estimates)

def bootknife_D3(A, B, C, n_replicates):
    D3_estimates = [] #list for distribution of D3 statistics
    num_blocks = len(A) // BLOCK_SIZE
    omit_block = np.random.randint(num_blocks)

    for i in range(n_replicates):
        A_bootstrap_list = []
        B_bootstrap_list = []
        C_bootstrap_list = []

        #sample bootstrap windows
        while len(A_bootstrap_list) < num_blocks:
            bootstrap_index = np.random.randint(num_blocks) #index to sample alignment
            if omit_block  != bootstrap_index:
                A_bootstrap_list.append(A[bootstrap_index * BLOCK_SIZE : (bootstrap_index + 1) * BLOCK_SIZE])
                B_bootstrap_list.append(B[bootstrap_index * BLOCK_SIZE : (bootstrap_index + 1) * BLOCK_SIZE])
                C_bootstrap_list.append(C[bootstrap_index * BLOCK_SIZE : (bootstrap_index + 1) * BLOCK_SIZE])
        print((len(A_bootstrap_list)))

        #concatenate windows
        A_bootstrap = ''.join(A_bootstrap_list)
        B_bootstrap = ''.join(B_bootstrap_list)
        C_bootstrap = ''.join(C_bootstrap_list)
        
        D3_val = calculate_D3(A_bootstrap, B_bootstrap, C_bootstrap)
        D3_estimates.append(D3_val)
    
    return calc_pval_from_bootstrap(D3_estimates)

def bootstrap_2_stage(A, B, C, n_replicates, locus_length):
    D3_estimates = [] #list for distribution of D3 statistics
    num_blocks = len(A) // locus_length

    for i in range(n_replicates):
        A_bootstrap_list = []
        B_bootstrap_list = []
        C_bootstrap_list = []

        #sample bootstrap windows
        while len(A_bootstrap_list) < num_blocks:
            bootstrap_index = np.random.randint(num_blocks) #index to sample alignment
            sub_A = A[bootstrap_index * locus_length : (bootstrap_index + 1) * locus_length]
            sub_B = B[bootstrap_index * locus_length : (bootstrap_index + 1) * locus_length]
            sub_C = C[bootstrap_index * locus_length : (bootstrap_index + 1) * locus_length]
            # res_A = ""
            # res_B = ""
            # res_C = ""
            # index_list = []
            # for j in range(locus_length):
            #     index = np.random.randint(locus_length)
            #     res_A += sub_A[index]
            #     res_B += sub_B[index]
            #     res_C += sub_C[index]
            indices = np.random.randint(locus_length, size=locus_length)
            res_A = np.frombuffer(sub_A.encode("utf8"), dtype = "uint8")[indices].tobytes().decode("utf8")
            res_B = np.frombuffer(sub_B.encode("utf8"), dtype = "uint8")[indices].tobytes().decode("utf8")
            res_C = np.frombuffer(sub_C.encode("utf8"), dtype = "uint8")[indices].tobytes().decode("utf8")
            
            A_bootstrap_list.append(res_A)
            B_bootstrap_list.append(res_B)
            C_bootstrap_list.append(res_C)
        print((len(A_bootstrap_list)), len(A_bootstrap_list[0]))

        #concatenate windows
        A_bootstrap = ''.join(A_bootstrap_list)
        B_bootstrap = ''.join(B_bootstrap_list)
        C_bootstrap = ''.join(C_bootstrap_list)
        
        D3_val = calculate_D3(A_bootstrap, B_bootstrap, C_bootstrap)
        D3_estimates.append(D3_val)
    
    return calc_pval_from_bootstrap(D3_estimates)


def compute_D3_concat(marker_path, locus_length, remedy = 0):
    with open(marker_path, "r") as handle:
        D3_list = {"A":[], "B":[], "C":[], "D3": [], "D3_mean": [], "D3_stdev":[], "D3_pval": []}
        seqs_list = {}
        for x in INGROUP_SPECIES:
            seqs_list[x] = ""
        lines = handle.readlines()
        triplet_list = list(combinations(INGROUP_SPECIES, 3))
        for line in lines:
            arr = line.split(" ")
            if arr[0] in INGROUP_SPECIES:
                seqs_list[arr[0]] += arr[1].strip()

        for triplet in triplet_list:
            D3_val = calculate_D3(seqs_list[triplet[0]], seqs_list[triplet[1]], seqs_list[triplet[2]])
            if remedy == 0:
                D3_mean, D3_stdev, D3_pval = bootstrap_D3(seqs_list[triplet[0]], seqs_list[triplet[1]], seqs_list[triplet[2]], BS_NUM_REPS)
            elif remedy == 1:
                D3_mean, D3_stdev, D3_pval = bootstrap_1_D3(seqs_list[triplet[0]], seqs_list[triplet[1]], seqs_list[triplet[2]], BS_NUM_REPS)
            elif remedy == 2:
                D3_mean, D3_stdev, D3_pval = bootknife_D3(seqs_list[triplet[0]], seqs_list[triplet[1]], seqs_list[triplet[2]], BS_NUM_REPS)
            elif remedy == 3:
                D3_mean, D3_stdev, D3_pval = bootstrap_2_stage(seqs_list[triplet[0]], seqs_list[triplet[1]], seqs_list[triplet[2]], BS_NUM_REPS, locus_length)
            D3_list["A"].append(triplet[0])
            D3_list["B"].append(triplet[1])
            D3_list["C"].append(triplet[2])
            D3_list["D3"].append(D3_val)
            D3_list["D3_mean"].append(D3_mean)
            D3_list["D3_stdev"].append(D3_stdev)
            D3_list["D3_pval"].append(D3_pval)
            
        df = pd.DataFrame(D3_list)
        if remedy == 0:
            df.to_csv(str(Path(marker_path).parent.absolute())+"/D3.txt")
            # df.to_csv(str(Path(marker_path).parent.absolute())+"/D3_"+str(BLOCK_SIZE)+".txt")
        elif remedy == 1:
            df.to_csv(str(Path(marker_path).parent.absolute())+"/D3_1.txt")
        elif remedy == 2:
            df.to_csv(str(Path(marker_path).parent.absolute())+"/D3_bootknife.txt")
        elif remedy == 3:
            df.to_csv(str(Path(marker_path).parent.absolute())+"/D3_2stage.txt")
    return df



def run_all(dir_path, re_end, exp_type, null_reti_num):
    # dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/"
    re_start = 1
    # re_end = 10
    # true_gt = False

    res_data = {"scale":[], "locus_length": [], "replicate":[], "true_gt":[], "pvalue_type":[], "method":[], "significant": []}
    for scale in ["short", "medium", "long"]:
        for locus_length in [500, 2000]:
            for i in range(re_start, re_end+1):
                marker_path = dir_path + scale+"/"+str(i)+"/heter/"+str(locus_length)+"/markers.txt"
                df = compute_D3_concat(marker_path)
                

if __name__ == "__main__":
    # path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/long/1/heter/500/markers.txt"
    # compute_D3_concat(path)

    compute_D3_concat(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))

