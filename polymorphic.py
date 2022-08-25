import matplotlib.pyplot as plt
import seaborn as sns
import sys
import pandas as pd
import dendropy
from dendropy.calculate import treecompare

# SEQ_LEN = int(sys.argv[4])

SEQ_LEN = 0

def read_markers(marker_path):
    with open(marker_path, "r") as handle:
        markers = []
        locus_index = 0
        for line in handle.readlines():
            if line.startswith("["):
                locus_index += 1
                markers.append({})
               
            else:
                arr = line.split(" ")
                markers[locus_index-1][arr[0]] = arr[1]
    
        return markers


def diff_letters(a, b):
    return sum ( a[i] != b[i] for i in range(len(a)) )

def diff_letters_array(seqs):
    cnt = 0
    for x in range(min(SEQ_LEN, len(seqs[0]))):
        diff = False
        for i in range(len(seqs) - 1):
            j = i + 1
            while j < len(seqs):
                if seqs[i][x] != seqs[j][x]:
                    cnt += 1
                    diff = True
                    break
                j += 1
            if diff:
                break
    return cnt 

def compute_pdistance_pairs(markers, species):
    diffs = {"species_i_j":[], "p-distance":[]}
    for marker in markers:
        for i in range(0, len(species)-1):
            for j in range(i+1, len(species)):
                diff_cnt = diff_letters(marker[species[i]], marker[species[j]])
                diffs["species_i_j"].append(species[i]+"-"+species[j])
                diffs["p-distance"].append(diff_cnt/SEQ_LEN)
    return diffs


def compute_polymorphism_array(markers, species):
    # print(species)
    diffs = []
    for marker in markers:
        seqs = [marker[x] for x in species]
        cnt = diff_letters_array(seqs)
        diffs.append(cnt/SEQ_LEN)
    return diffs

def polymorphism_all_replicate(dir_markers, ntaxa, num_replicate, sub_dir=None, heter=1):
    # dir_markers = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa3/locus_genebranch-0.1_0.5/"
    # dir_markers = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa3_tall/01/"
    data = {"replicate":[], "diff": []}
    if ntaxa == 5:
        species = ['R_0', 'Q_0', 'L_0', 'C_0', 'G_0']
    elif ntaxa == 4: 
        # species = ['R_0', 'Q_0', 'L_0', 'C_0']
        species = ['A_0', 'Q_0', 'L_0', 'G_0']
    else:
        species = ['L_0', 'A_0', 'Q_0']
    diffs = []
    # heter = 1
    for i in range(1, num_replicate+1):
        print(i)
        if heter == 1:
            if sub_dir is not None:
                marker_path = dir_markers + str(i) +"/heter/" + sub_dir+"/markers.txt"
            else:
                marker_path = dir_markers + str(i) + "/heter/markers.txt"
        elif heter == 0:
            if sub_dir is not None:
                marker_path = dir_markers + str(i) +"/homo/" + sub_dir+"/markers.txt"
            else:
                marker_path = dir_markers + str(i) + "/homo/markers.txt"
        elif heter == 2:
            marker_path = dir_markers + str(i) + "/heter_locus/markers.txt"
        else:
            marker_path = dir_markers + str(i) + "/markers.txt"

        markers = read_markers(marker_path)
        diffs_i = compute_polymorphism_array(markers, species)
        for x in range(len(diffs_i)):
            data["replicate"].append(i)
            data["diff"].append(diffs_i[x])
        
        diffs.extend(diffs_i)
    res = pd.DataFrame(data)
    plt.clf()
    sns.histplot(x=diffs, bins=30, stat="probability")
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel("Proportion of polymorphic sites", fontsize=16)
    plt.ylabel("Proportion", fontsize=16)
    plt.tight_layout()
    if heter == 1:
        if sub_dir is not None:
            res.to_csv(dir_markers+"/polymorphic_heter_"+sub_dir+".csv")
            plt.savefig(dir_markers+"polymorphism_heter_"+sub_dir+".pdf")
        else:
            plt.savefig(dir_markers+"polymorphism_heter.pdf")
            res.to_csv(dir_markers+"/polymorphic_heter.csv")
    elif heter == 0:
        plt.savefig(dir_markers+"polymorphism_homo.pdf")
        res.to_csv(dir_markers+"/polymorphic_homo.csv")
    elif heter == 2:
        plt.savefig(dir_markers+"polymorphism_heter_locus.pdf")
        res.to_csv(dir_markers+"/polymorphic_heter_locus.csv")
    else:
        plt.savefig(dir_markers+"polymorphism.pdf")
        
def p_dist_all_replicate(dir_markers, ntaxa, num_replicate, sub_dir="", heter = 1):
    data = {"replicate":[], "species_i_j": [], "p-distance": []}
    if ntaxa == 5:
        species = ['R_0', 'Q_0', 'L_0', 'C_0', 'G_0']
    elif ntaxa == 4:
        # species = ['R_0', 'Q_0', 'L_0', 'C_0']
        species = ['A_0', 'Q_0', 'L_0', 'G_0']
    else:
        species = ['L_0', 'A_0', 'Q_0']
    # heter = 1
    for i in range(1, num_replicate+1):
        print(i)
        if heter == 1:
            if sub_dir is not None:
                marker_path = dir_markers + str(i) +"/heter/" + sub_dir+"/markers.txt"
            else:
                marker_path = dir_markers + str(i) + "/heter/markers.txt"
        elif heter == 0:
            marker_path = dir_markers + str(i) + "/homo/"+ sub_dir+"/markers.txt"
        elif heter == 2:
            marker_path = dir_markers + str(i) + "/heter_locus/markers.txt"
        else:
            marker_path = dir_markers + str(i) + "/markers.txt"

        markers = read_markers(marker_path)
        diffs_i = compute_pdistance_pairs(markers, species)
        data["replicate"].extend([i for x in range(len(diffs_i["p-distance"]))])
        data["species_i_j"].extend(diffs_i["species_i_j"])
        data["p-distance"].extend(diffs_i["p-distance"])
    
    df = pd.DataFrame.from_dict(data)  
    plt.clf()
    ax = sns.violinplot(x="species_i_j", y="p-distance", data=df)
    plt.ylim(-0.05, 0.2)
    # sns.histplot(x=diffs, bins=30, stat="probability")
    # plt.xticks(fontsize=14)
    # plt.yticks(fontsize=14)
    # plt.xlabel("Proportion of polymorphic sites", fontsize=16)
    # plt.ylabel("Proportion", fontsize=16)
    plt.tight_layout()
    if heter == 1:
        if sub_dir is not None:
            plt.savefig(dir_markers+"pdist_heter_"+sub_dir+".pdf")
            df.to_csv(dir_markers+"pdist_heter_"+sub_dir+".csv")
        else:
            plt.savefig(dir_markers+"pdist_heter.pdf")
    elif heter == 0:
        plt.savefig(dir_markers+"pdist_homo.pdf")
        df.to_csv(dir_markers+"pdist_homo_.csv")
    elif heter == 2:
        plt.savefig(dir_markers+"pdist_heter_locus.pdf")
        df.to_csv(dir_markers+"pdist_heter_.csv")
    else:
        plt.savefig(dir_markers+"pdist.pdf")

def plot_hybrid_coalescence():
    # csv_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/net1_2/pdist_heter_500.csv"
    csv_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/net1/medium/pdist_heter_500.csv"
    df = pd.read_csv(csv_path)
    fig, axes = plt.subplots(1, 10, figsize=(30, 3))
    for rep in range(10):
        df_cur = df[df["replicate"] == rep+1]
        sns.histplot(y="p-distance", data=df_cur[(df_cur["species_i_j"] =="L_0-Q_0") | (df_cur["species_i_j"] =="Q_0-L_0")], ax=axes[rep])
        axes[rep].set_ylim(0, 0.4)
    plt.tight_layout()
    # plt.savefig("/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/net1_2/pdist_LQ.pdf")
    plt.savefig("/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/net1/medium/pdist_LQ.pdf")
    plt.show()

if __name__ == "__main__":
    # marker_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_gt/1000/1/heter/018/markers.txt"
    # marker_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_gt/1000/1/homo/markers.txt"

    # markers = read_markers(marker_path)
    # species = ['R_0', 'Q_0', 'L_0', 'C_0', 'G_0']
    # diff = compute_polymorphism_array(markers, species)
    # print(diff)
    # sns.histplot(x=diff, bins=30, stat="probability")
    # plt.show()
    

    # polymorphism_all_replicate(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])

    for scale in ["medium", "long"]:
        marker_dir = sys.argv[1]+scale+"/"
        n_taxa = int(sys.argv[5])
        num_replicate = int(sys.argv[2])
        locus_length = sys.argv[3]
        SEQ_LEN = int(locus_length)
        heter = int(sys.argv[4])
        polymorphism_all_replicate(marker_dir, n_taxa, num_replicate, locus_length, heter)
        p_dist_all_replicate(marker_dir, n_taxa, num_replicate, locus_length, heter)
    
    # marker_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/data/simulation/seqgen/replica/net/"
    # num_replicate = 10
    # n_taxa=4
    # SEQ_LEN=500
    # heter=0
    # polymorphism_all_replicate(marker_dir, n_taxa, num_replicate, heter=heter)
    # p_dist_all_replicate(marker_dir, n_taxa, num_replicate, heter=heter)
    # plot_hybrid_coalescence()