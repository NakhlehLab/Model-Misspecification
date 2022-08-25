import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import dendropy
from dendropy.calculate import treecompare
from scipy.special import rel_entr
import scipy.stats as stats
import sys
import os
from scipy.stats import power_divergence

from calculator import *

N_TAXA = 5
DEBUG=False
if DEBUG:
    NUM_REPLICATE = 3
else:
    NUM_REPLICATE = 100
NUM_GT = 10000
NUM_GT_TOPO = 105
gt_name = {True: "truegt", False: "iqgt"}
DOF = NUM_GT_TOPO - N_TAXA + 1
DDOF = N_TAXA - 2

def get_obs(tree_path, outgroup=None):
    tree_list = dendropy.TreeList.get(path=tree_path, schema="newick", rooting="force-rooted")
    pruned_tree_list = []
    for tree in tree_list:
        tree.prune_taxa_with_labels([outgroup])
        pruned_tree_list.append(tree)
    pruned_tree_list = dendropy.TreeList(pruned_tree_list)
    unique_topologies = pruned_tree_list.as_tree_array().topologies(sort_descending = True, frequency_attr_name='frequency')
    return unique_topologies

def get_exp(csv_path, taxon_map=None):
    df = pd.read_csv(csv_path)
    exp_dict = df.set_index('gt').to_dict()['prob_exp']
    return exp_dict

def map_obs_exp(obs, exp):
    data = {"gt":[], "obs": [], "exp": []}
    treelist_string = "".join(exp.keys())
    treelist_exp = dendropy.TreeList.get(data=treelist_string, rooting="force-rooted", schema="newick")
    i = 0
    for tree1 in treelist_exp:
        i += 1
        tns = obs[0].taxon_namespace
        tree11 = dendropy.Tree.get(data=tree1.as_string(schema="newick"), rooting="force-rooted", schema="newick", taxon_namespace=tns)
        tree1str = tree1.as_string(schema="newick")[5:].strip()
        find = False
        for tree2 in obs:
            # print(tree2.as_string(schema="newick"))
            if treecompare.symmetric_difference(tree11, tree2) == 0:
                data["gt"].append(tree1str)
                data["exp"].append(exp.get(tree1str))
                data["obs"].append(tree2.frequency) 
                find = True  
                break
        if not find:
            data["gt"].append(tree1str)
            data["exp"].append(exp.get(tree1str))
            data["obs"].append(0)
             
    # print(data)
    return data

def process_data(obs_path, exp_path, outgroup):
    unique_topologies = get_obs(obs_path, outgroup)
    exp_dict = get_exp(exp_path)
    data = map_obs_exp(unique_topologies, exp_dict)
    print(sum(data["exp"]))
    return data
    # return compute_Pvalue(data["obs"], data["exp"], N_TAXA-2)

def compute_KL(obs, exp):
    KL = rel_entr(obs, exp)
    # print(f"KL={sum(KL)}")
    return sum(KL)

def get_rf(iqtree_error_path):
    data = pd.read_csv(iqtree_error_path)
    return data["RF"]

# def get_ML_species_tree(path):
#     with open(path, "r") as handle:
#         begin = False
#         for line in handle.readlines():
#             if line.startswith("Inferred Network #1"):
#                 begin = True
#             elif begin:
#                 print(line)
#                 inferred_tree = dendropy.Tree.get(data=line, schema="newick", rooting="force-rooted")
#                 break


def compute_heter_replica(obs_dir, exp_path, outgroup, csv_path, iqtree_error_path=None, dg="KL", plot="pvalue"):
    gt_errors = pd.read_csv(csv_path)
    print(gt_errors)
    obs_list = []
    exp_list = []
    res_dict = {"KL":[], "RF": [], "pvalue":[], "chi2":[]}

    for i in range(1, NUM_REPLICATE+1):
        # x_list.append(gt_errors[gt_errors["index"] == i]["RF"].values[0])
        obs_path = obs_dir + "/" + str(i) + "/heter/rooted_iqtree.txt"
        data = process_data(obs_path, exp_path, outgroup)
        print(data)
        print("*****max*******")
        print(max(data["obs"]))
        print(max(data["exp"]))

        # if dg == "KL":
        #     dg_list.append(compute_KL(data["obs"], data["exp"]))
        res_dict["KL"].append(compute_KL(data["obs"], data["exp"]))

        obs = [NUM_GT*i for i in data["obs"]]
        exp = [NUM_GT*i for i in data["exp"]]
        print("--------------")
        print(obs)
        print(exp)
        
        
        chi2, pvalue = compute_Pvalue(obs, exp, DDOF)
        res_dict["chi2"].append(chi2)
        res_dict["pvalue"].append(pvalue)
        obs_list.extend(data["obs"])
        exp_list.extend(data["exp"])
        # y_list.append(pvalue)
    
    # if dg == "RF" and iqtree_error_path is not None:
    #     dg_list = get_rf(iqtree_error_path)
    res_dict["RF"] = get_rf(iqtree_error_path)

    df = pd.DataFrame(res_dict)
    print(df)

    if plot == "pvalue":
        fig = sns.scatterplot(data=df, x=dg, y="pvalue", ci=None)
        fig.axhline(0.05)
        plt.ylabel("p-value")
        plt.xlabel(dg)
        plt.savefig(obs_dir+"/pvalue_truest_"+dg+".pdf")
        plt.show()
    elif plot == "gtprob":
        index_list = np.argsort(obs_list)
        obs_list = [obs_list[i] for i in index_list]
        exp_list = [exp_list[i] for i in index_list]
        fig = sns.scatterplot(x=obs_list, y=exp_list, s=5)
        sns.lineplot( x=obs_list, y=obs_list,color="blue")

        # fig = sns.scatterplot(x=range(10*105), y=exp_list, s=2)
        plt.ylabel("Expected gene tree probability")
        plt.xlabel("Observed gene tree probability")

        plt.savefig(obs_dir+"/gtprob_scatter_truest"+".pdf")
        plt.show()
    

def true_st_iqgt_x2():
    # obs_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5/1/heter/rooted_iqtree.txt"
    # exp_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/experiment/5taxa/gtprob.csv"
    # taxon_map = {"Q":"4", "R":"5", "C":"3", "G":"2", "L":"1"}
    # taxon_map = {"4":"Q_0", "5":"R_0", "3":"C_0", "2":"G_0", "1":"L_0"}
    # compute_sigificance(obs_path, exp_path, outgroup="Z 0", taxon_map=taxon_map)

    obs_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_tall10/036 2/"
    exp_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/experiment/5taxa/gtprob.csv"
    csv_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5/iqtree_error.csv"
    iq_err_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5/iqtree_error.csv"
    compute_heter_replica(obs_dir, exp_path, "Z 0", csv_path, iq_err_path, "KL", "gtprob")


def compute_heter_inferred_replica(obs_dir, true_gt, outgroup, sub_dir="", iqtree_error_path=None):
    if iqtree_error_path is not None:
        res_dict = {"KL":[], "RF": [], "chi2":[], "pvalue-chi2":[], "g-test":[], "pvalue-g":[]}
        res_dict["RF"] = get_rf(iqtree_error_path)
    else:
        res_dict = {"KL":[], "chi2":[], "pvalue-chi2":[], "g-test":[], "pvalue-g":[]}
    exp_list = []
    obs_list = []
    for i in range(1, NUM_REPLICATE + 1):
        if true_gt:
            obs_path = obs_dir + str(i)+"/heter/"+sub_dir+"/genetrees.txt"
            exp_path = obs_dir + str(i)+"/heter/"+sub_dir+"/gtprob_true_0.csv"
        else:
            obs_path = obs_dir + str(i)+"/heter/"+sub_dir+"/rooted_iqtree.txt"
            exp_path = obs_dir + str(i)+"/heter/"+sub_dir+"/gtprob_iq_0.csv"
        data = process_data(obs_path, exp_path, outgroup)
        

        obs = [NUM_GT*i for i in data["obs"]]
        exp = [NUM_GT*i for i in data["exp"]]
        print("-----------DDOF------------")
        print(DDOF)
        chi2, pvalue = compute_Pvalue(obs, exp, DDOF)
        g_test, p_g = power_divergence(obs, exp, DDOF, lambda_='log-likelihood')
        obs_list.extend(data["obs"])
        exp_list.extend(data["exp"])
        # print(data)
        res_dict["KL"].append(compute_KL(data["obs"], data["exp"]))
        res_dict["chi2"].append(chi2)
        res_dict["pvalue-chi2"].append(pvalue)  
        res_dict["g-test"].append(g_test)
        res_dict["pvalue-g"].append(p_g)
    print(res_dict)
    df = pd.DataFrame(res_dict)
    print(df)
    df.to_csv(obs_dir+"/res_"+gt_name[true_gt]+"_"+sub_dir+".csv")
    return obs_list, exp_list
  
def plot_figures(obs_dir, true_gt, plot, dg, obs_list, exp_list, rate_dir=""):
    res_csv_path = obs_dir+"/res_"+gt_name[true_gt]+"_"+rate_dir+".csv"
    df = pd.read_csv(res_csv_path)

    plt.clf()
    if plot == "p-dg":
        fig = sns.scatterplot(data=df, x=dg, y="pvalue-chi2", ci=None)
        fig.axhline(0.05)
        plt.ylabel("p-value", fontsize=16)
        plt.xlabel(dg, fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.tight_layout()
        plt.savefig(obs_dir+"/pvalue_ml_"+gt_name[true_gt]+"_"+dg+"_"+rate_dir+".pdf")
        
    elif plot == "data":
        index_list = np.argsort(obs_list)
        obs_list = [obs_list[i] for i in index_list]
        exp_list = [exp_list[i] for i in index_list]
        fig = sns.scatterplot(x=range(NUM_REPLICATE*NUM_GT_TOPO), y=obs_list, s=2)
        fig = sns.scatterplot(x=range(NUM_REPLICATE*NUM_GT_TOPO), y=exp_list, s=2)
        plt.ylabel("Gene tree probability", fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.tight_layout()
        plt.savefig(obs_dir+"/data_ml_"+gt_name[true_gt]+"_"+rate_dir+".pdf")

    elif plot == "data2":
        index_list = np.argsort(obs_list)
        obs_list = [obs_list[i] for i in index_list]
        exp_list = [exp_list[i] for i in index_list]
        sns.lineplot(x=obs_list, y=obs_list,color="pink")
        fig = sns.scatterplot(x=obs_list, y=exp_list, s=5)
        plt.ylabel("Expected gene tree probability", fontsize=16)
        plt.xlabel("Observed gene tree probability", fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.tight_layout()
        plt.savefig(obs_dir+"/data_ml_scatter_"+gt_name[true_gt]+"_"+rate_dir+".pdf")

    elif plot == "chi2":
        fig = sns.histplot(data=df, x="chi2", bins=30, stat="probability", log_scale=True, palette="Set2")
        x = np.arange(0, 1000, .0001)
        plt.plot(x, stats.chi2.pdf(x, df=DOF), color='r', lw=2)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel("$\chi^2$", fontsize=16)
        plt.ylabel("Proportion", fontsize=16)
        plt.tight_layout()
        plt.savefig(obs_dir+"/chi2_"+gt_name[true_gt]+"_"+rate_dir+".pdf")

    elif plot == "p-chi2":
        fig = sns.histplot(data=df, x="pvalue-chi2", bins=20, stat="probability")
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel("P-value", fontsize=16)
        plt.ylabel("Proportion", fontsize=16)
        # plt.xlim(0,1)
        plt.tight_layout()
        plt.savefig(obs_dir+"/pvalue-chi2_"+gt_name[true_gt]+"_"+rate_dir+".pdf")

    elif plot == "g-test":
        fig = sns.histplot(data=df, x="g-test", bins=30, stat="probability", palette="Set2")
        x = np.arange(0, 1000, .0001)
        plt.plot(x, stats.chi2.pdf(x, df=DOF), color='r', lw=2)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel("G-test", fontsize=16)
        plt.ylabel("Proportion", fontsize=16)
        plt.tight_layout()
        plt.savefig(obs_dir+"/Gtest_"+gt_name[true_gt]+"_"+rate_dir+".pdf")

    elif plot == "pvalue-g":
        fig = sns.histplot(data=df, x="pvalue-g", bins=20, stat="probability")
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel("P-value", fontsize=16)
        plt.ylabel("Proportion", fontsize=16)
        # plt.xlim(0,1)
        plt.tight_layout()
        plt.savefig(obs_dir+"/pvalue-gtest_"+gt_name[true_gt]+"_"+rate_dir+".pdf")


def ML_st_iqgt_x2(obs_dir, sub_dir=""):
    # obs_list, exp_list = compute_heter_inferred_replica(obs_dir, True, "Z 0","", None)
    # plot_figures(obs_dir, True, "g-test", "KL", obs_list, exp_list)
    # plot_figures(obs_dir, True, "pvalue-g", "KL", obs_list, exp_list)
    # plot_figures(obs_dir, True, "p-chi2", "KL", obs_list, exp_list)
    # plot_figures(obs_dir, True, "chi2", "KL", obs_list, exp_list)
    # plot_figures(obs_dir, True, "data", "KL", obs_list, exp_list)
    # plot_figures(obs_dir, True, "data2", "KL", obs_list, exp_list)

    obs_list, exp_list = compute_heter_inferred_replica(obs_dir, False, "Z 0",sub_dir, None)
    plot_figures(obs_dir, False, "g-test", "KL", obs_list, exp_list, sub_dir)
    plot_figures(obs_dir, False, "pvalue-g", "KL", obs_list, exp_list, sub_dir)
    plot_figures(obs_dir, False, "p-chi2", "KL", obs_list, exp_list, sub_dir)
    plot_figures(obs_dir, False, "chi2", "KL", obs_list, exp_list, sub_dir)
    plot_figures(obs_dir, False, "data", "KL", obs_list, exp_list, sub_dir)
    plot_figures(obs_dir, False, "data2", "KL", obs_list, exp_list, sub_dir)



def true_st_truegt_x2():
    # obs_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5/1/heter/rooted_iqtree.txt"
    # exp_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/experiment/5taxa/gtprob.csv"
    # taxon_map = {"4":"Q_0", "5":"R_0", "3":"C_0", "2":"G_0", "1":"L_0"}
    # compute_sigificance(obs_path, exp_path, outgroup="Z 0", taxon_map=taxon_map)

    obs_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5/"
    exp_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/experiment/5taxa/gtprob.csv"
    csv_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5/iqtree_error.csv"
    iq_err_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5/iqtree_error.csv"
    compute_heter_replica(obs_dir, exp_path, "Z 0", csv_path, iq_err_path, "KL")

if __name__ == "__main__":
    # true_st_iqgt_x2()
    if len(sys.argv) == 3:
        ML_st_iqgt_x2(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 2:
        ML_st_iqgt_x2(sys.argv[1])

    