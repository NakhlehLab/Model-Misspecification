import dendropy
from calculator import *
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from calculator import compute_Pvalue
import scipy.stats as stats
from scipy.special import rel_entr

gt_name = {True: "truegt", False:"infergt"}
N_TAXA = 3
NUM_GT_TOP = 3
NUM_REPLICATE = 100
def compute_KL(obs, exp):
    KL = rel_entr(obs, exp)
    print(f"KL={sum(KL)}")
    return sum(KL)

def compute_freqs(tree_path, outgroup=None):
    tree_list = dendropy.TreeList.get(path=tree_path, schema="newick", rooting="force-rooted")
    pruned_tree_list = []
    for tree in tree_list:
        tree.prune_taxa_with_labels([outgroup])
        pruned_tree_list.append(tree)

    pruned_tree_list = dendropy.TreeList(pruned_tree_list)
    unique_topologies = pruned_tree_list.as_tree_array().topologies(sort_descending = True, frequency_attr_name='frequency')
    freqs = []
    for tree in unique_topologies:
        freqs.append(tree.frequency)
    return freqs

def get_obs_exp(tree_path, outgroup):
    obs = compute_freqs(tree_path, outgroup)
    obs = [x*100 for x in obs]
    t = internal_time(obs[0]/sum(obs))
    expct = compute_exp(t)
    exps = [x*sum(obs) for x in expct]
    # chi2, pvalue = chisquare(obs, exps, ddof=1)
    # print(f"pvalue={pvalue}")
    print(f"obs={obs}, exps={exps}")
    return obs, exps, t

    # obs = compute_freqs(tree_path, outgroup)
    # t = internal_time(obs[0])
    # exps = compute_exp(t)
    # print(f"obs={obs}, exps={exps}")
    # return obs, exps

def compute_all_replica(dirct, outgroup, csv_path, dg="KL"):
    scenario="net"
    csv_data = pd.read_csv(csv_path)
    gt_errors = csv_data[(csv_data["type"]==scenario) & (csv_data["homo"]==False)]
    print(gt_errors)
    
    
    res_dict = {"KL":[], "RF": [], "pvalue":[], "chi2":[]}
    for i in range(1, 11):
        path = dirct + scenario + "/" + str(i) + "/heter_locus/rooted_iqtree.txt"
        obs, exp = get_obs_exp(path, outgroup=outgroup)
        res_dict["KL"].append(compute_KL(obs, exp))
        chi2, pvalue = compute_Pvalue(obs, exp, 1)
        res_dict["chi2"].append(chi2)
        res_dict["pvalue"].append(pvalue)
    res_dict["RF"] = list(gt_errors["RF"])
    df = pd.DataFrame(res_dict)
    print(df)

    # fig = sns.scatterplot(x=x_list, y=y_list, ci=None)
    # fig = plt.plot(x_list, y_list)
    fig = sns.scatterplot(data=df, x=dg, y="pvalue", ci=None)
    fig.axhline(0.05)
    # plt.xlabel(dg)
    # plt.ylabel("p-value")
    plt.savefig("/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/experiment/heter/"+scenario+"_"+dg+".pdf")
    plt.show()

def compute_heter_replica(obs_dir, outgroup, csv_path, dg="KL", true_gt=True, plot = "time"):
    gt_errors = pd.read_csv(csv_path)
    print(gt_errors)
    # x_list = []
    # y_list = []
    res_dict = {"KL":[], "RF": [], "time":[], "pvalue":[], "chi2":[]}
    obs_list = []
    exp_list = []
    for i in range(1, NUM_REPLICATE + 1):
        if true_gt:
            path = obs_dir + "/" + str(i) + "/heter/genetrees.txt"
        else:
            path = obs_dir + "/" + str(i) + "/heter/rooted_iqtree.txt"
        obs, exp, t = get_obs_exp(path, outgroup=outgroup)
        res_dict["KL"].append(compute_KL(obs, exp))
        res_dict["time"].append(1.0-t)
        chi2, pvalue = compute_Pvalue(obs, exp, 1)
        res_dict["chi2"].append(chi2)
        res_dict["pvalue"].append(pvalue)
        obs_list.extend(obs)
        exp_list.extend(exp)
        
    res_dict["RF"] = list(gt_errors["RF"])
    df = pd.DataFrame(res_dict)
    print(df)
    df.to_csv(obs_dir+"/res_"+gt_name[true_gt]+".csv")

    plt.clf()
    if plot == "p-dg":
        fig = sns.scatterplot(data=df, x=dg, y="pvalue", ci=None)
        fig.axhline(0.05)
        plt.ylabel("p-value", fontsize=16)
        plt.xlabel(dg, fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.tight_layout()
        plt.savefig(obs_dir+"/pvalue_"+gt_name[true_gt]+"_"+dg+".pdf")
    #     fig = sns.scatterplot(data=df, x=dg, y="pvalue", ci=None)
    #     fig.axhline(0.05)
        
    #     plt.savefig(obs_dir+"/pvalue_"+dg+"_"+truegt_name[truegt]+".pdf")
    #     plt.show()
    elif plot == "time":
        fig = sns.scatterplot(data=df, x=dg, y="time", ci=None)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.savefig(obs_dir+"/brl_"+gt_name[true_gt]+".pdf")
        # plt.show()
    elif plot == "time2":
        # fig = sns.scatterplot(data=df, x=dg, y="time", ci=None)
        df.boxplot(column=["time"])
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.savefig(obs_dir+"/brl_"+gt_name[true_gt]+"_box.pdf")
        # plt.show()
    elif plot == "data":
        # index_list = np.argsort(obs_list)
        # obs_list = [obs_list[i] for i in index_list]
        # exp_list = [exp_list[i] for i in index_list]
        # fig = sns.scatterplot(x=range(NUM_REPLICATE*NUM_GT_TOP), y=obs_list, s=10)
        # fig = sns.scatterplot(x=range(NUM_REPLICATE*NUM_GT_TOP), y=exp_list, s=10)
        # fig.set(ylim=(0, 100))

        # # plt.scatter(obs, exp)
        # plt.savefig(dirct+"/obs_exp_"+truegt_name[truegt]+".pdf")
        # plt.show()

        index_list = np.argsort(obs_list)
        obs_list = [obs_list[i] for i in index_list]
        exp_list = [exp_list[i] for i in index_list]
        fig = sns.scatterplot(x=range(NUM_REPLICATE*NUM_GT_TOP), y=obs_list, s=2)
        fig = sns.scatterplot(x=range(NUM_REPLICATE*NUM_GT_TOP), y=exp_list, s=2)
        plt.ylabel("Gene tree probability", fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.tight_layout()
        plt.savefig(obs_dir+"/data_"+gt_name[true_gt]+".pdf")
    elif plot == "data2":
        # index_list = np.argsort(obs_list)
        # obs_list = [obs_list[i] for i in index_list]
        # exp_list = [exp_list[i] for i in index_list]
        # fig = sns.scatterplot(x=obs_list, y=exp_list, s=5)
        # sns.lineplot(x=obs_list, y=obs_list,color="blue")

        # # fig = sns.scatterplot(x=range(10*105), y=exp_list, s=2)
        # plt.ylabel("Expected gene tree probability")
        # plt.xlabel("Observed gene tree probability")

        
        # plt.savefig(dirct+"/data_ml_scatter_"+truegt_name[truegt]+".pdf")

        index_list = np.argsort(obs_list)
        obs_list = [obs_list[i] for i in index_list]
        exp_list = [exp_list[i] for i in index_list]
        sns.lineplot(x=obs_list, y=obs_list,color="pink", alpha=0.5)
        fig = sns.scatterplot(x=obs_list, y=exp_list, s=5)
        # fig = sns.scatterplot(x=range(10*105), y=exp_list, s=2)
        plt.ylabel("Expected gene tree probability", fontsize=16)
        plt.xlabel("Observed gene tree probability", fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.tight_layout()

        plt.savefig(obs_dir+"/data_scatter_"+gt_name[true_gt]+".pdf")

    elif plot == "X2":
        # fig = sns.histplot(data=df, x="chi2", bins=30)
        # plt.savefig(dirct+"/chi2_"+truegt_name[truegt]+".pdf")

        fig = sns.histplot(data=df, x="chi2", bins=30, stat="probability", palette="Set2")
        x = np.arange(0, 10, .01)
        plt.plot(x, stats.chi2.pdf(x, df=N_TAXA-2), color='r', lw=2)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel("$\chi^2$", fontsize=16)
        plt.ylabel("Proportion", fontsize=16)
        plt.ylim([0, 1])
        plt.xlim(0, )

        plt.tight_layout()
        plt.savefig(obs_dir+"/chi2_"+gt_name[true_gt]+".pdf")

    elif plot=="pvalue":
        fig = sns.histplot(data=df, x="pvalue", bins=20, stat="probability")
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel("P-value", fontsize=16)
        plt.ylabel("Proportion", fontsize=16)
        plt.tight_layout()
        plt.savefig(obs_dir+"/pvalue_"+gt_name[true_gt]+".pdf")


def run(obs_dir):
    csv_path = obs_dir+"iqtree_error.csv"

    compute_heter_replica(obs_dir, outgroup="Z 0", csv_path=csv_path, dg="KL", true_gt=True, plot="time2")
    compute_heter_replica(obs_dir, outgroup="Z 0", csv_path=csv_path, dg="KL", true_gt=True, plot="data2")
    compute_heter_replica(obs_dir, outgroup="Z 0", csv_path=csv_path, dg="KL", true_gt=True, plot="data")
    compute_heter_replica(obs_dir, outgroup="Z 0", csv_path=csv_path, dg="KL", true_gt=True, plot="pvalue")
    compute_heter_replica(obs_dir, outgroup="Z 0", csv_path=csv_path, dg="KL", true_gt=True, plot="p-dg")
    compute_heter_replica(obs_dir, outgroup="Z 0", csv_path=csv_path, dg="RF", true_gt=True, plot="p-dg")
    compute_heter_replica(obs_dir, outgroup="Z 0", csv_path=csv_path, dg="KL", true_gt=True, plot="X2")

    compute_heter_replica(obs_dir, outgroup="Z 0", csv_path=csv_path, dg="KL", true_gt=False, plot="time2")
    compute_heter_replica(obs_dir, outgroup="Z 0", csv_path=csv_path, dg="KL", true_gt=False, plot="data2")
    compute_heter_replica(obs_dir, outgroup="Z 0", csv_path=csv_path, dg="KL", true_gt=False, plot="data")
    compute_heter_replica(obs_dir, outgroup="Z 0", csv_path=csv_path, dg="KL", true_gt=False, plot="pvalue")
    compute_heter_replica(obs_dir, outgroup="Z 0", csv_path=csv_path, dg="KL", true_gt=False, plot="p-dg")
    compute_heter_replica(obs_dir, outgroup="Z 0", csv_path=csv_path, dg="RF", true_gt=False, plot="p-dg")
    compute_heter_replica(obs_dir, outgroup="Z 0", csv_path=csv_path, dg="KL", true_gt=False, plot="X2")


    

if __name__ == "__main__":

    # csv_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/data/simulation/seqgen/replica/iqtree_err.csv"
    # dirct = "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/data/simulation/seqgen/replica/"
    # compute_all_replica(dirct, outgroup="G 0", csv_path=csv_path, dg="KL")


    # dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa3/locus_genebranch-0.1_0.5/"
    # dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa3_tall/036/"

    # csv_path = dir+"iqtree_error.csv"
    # compute_heter_replica(dir, outgroup="Z 0", csv_path=csv_path, dg="KL", truegt=True, plot="time2")

    # path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/data/simulation/seqgen/tree/heter/rooted_iqtree.txt"
    # print(compute_sigificance(path, outgroup="G 0"))

    run(sys.argv[1])