import dendropy
import numpy as np
from D3 import INGROUP_SPECIES, BLOCK_SIZE, SPECIES, OUTGROUP, calc_pval_from_bootstrap
from itertools import combinations
import pandas as pd
from pathlib import Path
import sys


NUM_LOCI = 10000
BS_NUM_REPS = 1000

def gt_calculate_D3(D_A_B, D_A_C, D_B_C):
    three_pairwise = [(D_A_B - D_A_C)/float(D_A_B + D_A_C + 0.0000001), (D_A_B - D_B_C)/float(D_A_B + D_B_C + 0.0000001), (D_A_C - D_B_C)/float(D_A_C + D_B_C + 0.0000001)]

    three_pairwise_abs = [abs(value) for value in three_pairwise] #Get absolute values of comparisons

    index_min = min(range(len(three_pairwise_abs)), key=three_pairwise_abs.__getitem__) #get index of D3 value

    D3 = three_pairwise[index_min] #D3 is the difference of smallest magnitude 

    return D3

def gt_bootstrap_D3(dist_list_ab, dist_list_ac, dist_list_bc, num_blocks, gt_window):
    D3_estimates = [] #list for distribution of D3 statistics

    for i in range(BS_NUM_REPS):
        ab_bootstrap_list = []
        ac_bootstrap_list = []
        bc_bootstrap_list = []

        #sample bootstrap windows
        for j in range(num_blocks):
            bootstrap_index = np.random.randint(num_blocks) #index to sample alignment
            ab_bootstrap_list.append(np.sum(dist_list_ab[bootstrap_index * gt_window : (bootstrap_index+1) * gt_window]))
            ac_bootstrap_list.append(np.sum(dist_list_ac[bootstrap_index * gt_window : (bootstrap_index+1) * gt_window]))
            bc_bootstrap_list.append(np.sum(dist_list_bc[bootstrap_index * gt_window : (bootstrap_index+1) * gt_window]))
            
        ab_bootstrap = np.sum(ab_bootstrap_list)
        ac_bootstrap = np.sum(ac_bootstrap_list)
        bc_bootstrap = np.sum(bc_bootstrap_list)
        
        D3_val = gt_calculate_D3(ab_bootstrap, ac_bootstrap, bc_bootstrap)
        D3_estimates.append(D3_val)
    D3_val_concat = gt_calculate_D3(np.sum(dist_list_ab), np.sum(dist_list_ac), np.sum(dist_list_bc))
    bs_reps_mean, bs_reps_std, p_value = calc_pval_from_bootstrap(D3_estimates)    

    return D3_val_concat,  bs_reps_mean, bs_reps_std, p_value
        


def compute_D3_gt(gt_path, locus_length,  theta_divide_2):
    tree_list = dendropy.TreeList.get(path=gt_path, schema="newick", rooting="force-rooted")
    distance_matrix_list = []
    num_blocks = NUM_LOCI * locus_length // BLOCK_SIZE
    gt_window = BLOCK_SIZE // locus_length

    # construct pairwise distance list 
    for tree in tree_list:
        pdc = tree.phylogenetic_distance_matrix()
        distance_matrix = {}
        for i, t1 in enumerate(tree.taxon_namespace[:-1]):
            for t2 in tree.taxon_namespace[i+1:]:
                # print(t1.label, t2.label, pdc(t1, t2))
                if t1.label < t2.label:
                    distance_matrix[t1.label.replace(" ", "_")+","+t2.label.replace(" ", "_")] = pdc(t1, t2) * theta_divide_2
                else:
                    distance_matrix[t2.label.replace(" ", "_")+","+t1.label.replace(" ", "_")] = pdc(t1, t2) * theta_divide_2
        distance_matrix_list.append(distance_matrix) 
        
    D3_list = {"A":[], "B":[], "C":[], "D3": [], "D3_mean": [], "D3_stdev":[], "D3_pval": []}
    triplet_list = list(combinations(INGROUP_SPECIES, 3))
    df_distance = pd.DataFrame(distance_matrix_list)

    for triplet in triplet_list:
        triplet_list = list(triplet)
        triplet_list.sort()
        # D3_val = gt_calculate_D3()
        D3_val, D3_mean, D3_stdev, D3_pval = gt_bootstrap_D3(df_distance.loc[:, triplet_list[0]+","+triplet_list[1]], 
        df_distance.loc[:, triplet_list[0]+","+triplet_list[2]],  df_distance.loc[:, triplet_list[1]+","+triplet_list[2]], num_blocks, gt_window)
        
        D3_list["A"].append(triplet_list[0])
        D3_list["B"].append(triplet_list[1])
        D3_list["C"].append(triplet_list[2])
        D3_list["D3"].append(D3_val)
        D3_list["D3_mean"].append(D3_mean)
        D3_list["D3_stdev"].append(D3_stdev)
        D3_list["D3_pval"].append(D3_pval)
        
    df = pd.DataFrame(D3_list)
    df.to_csv(str(Path(gt_path).parent.absolute())+"/D3_gt.txt")
    print(df)
    return df


if __name__ == "__main__":
    # path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/long/1/homo/genetrees.txt"
    # x = pd.read_csv(path)
    # compute_D3_gt(path, 500, 0.0005)
    compute_D3_gt(sys.argv[1], int(sys.argv[2]), float(sys.argv[3]))
    