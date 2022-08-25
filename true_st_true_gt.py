from eval_5taxa import process_data, DDOF, DOF
from calculator import compute_Pvalue
from scipy.stats import power_divergence
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
from pathlib import Path

def one_replicate(obs_path, exp_path, outgroup):
    data = process_data(obs_path, exp_path, outgroup)
    obs = [NUM_GT*i for i in data["obs"]]
    exp = [NUM_GT*i for i in data["exp"]]
    print(data)
    print(sum(obs), sum(exp))
    df = pd.DataFrame(data)
    df.to_csv(str(Path(obs_path).parent.absolute())+"/gt_prob_iqgt_MLst_obs_exp.csv")
    chi2, pvalue = compute_Pvalue(obs, exp, DDOF)
    g_test, p_g = power_divergence(obs, exp, DDOF, lambda_='log-likelihood')
    print(chi2, pvalue, g_test, p_g)
    return chi2, pvalue, g_test, p_g

def run_all(obs_dir, exp_path, outgroup, rate_dir):
    # obs_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_gt/10000/"
    # exp_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_gt/10000/gtprob_true_st.csv"
    # outgroup = "Z 0"
    # obs_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_abcde/10000/"
    # exp_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_abcde/10000/gtprob_true_st2.csv"
    # outgroup = "Z"

    chi2_list = []
    p_chi2_list = []
    g_list = []
    p_g_list = []
    for i in range(1, 101):
        # obs_path = obs_dir + str(i)+"/heter/genetrees.txt"
        # obs_path = obs_dir + str(i)+"/heter/genetrees.txt"
        # exp_path = obs_dir + str(i)+"/heter/gtprob_true_0.csv"

        obs_path = obs_dir + str(i)+"/heter/"+rate_dir+"/rooted_iqtree.txt"
        exp_path = obs_dir + str(i)+"/heter/"+rate_dir+"/gtprob_iq_0.csv"
        chi2, pvalue, g_test, p_g = one_replicate(obs_path, exp_path, outgroup)
    #     chi2_list.append(chi2)
    #     p_chi2_list.append(pvalue)
    #     g_list.append(g_test)
    #     p_g_list.append(p_g)
    # print(sum(chi2_list)/len(chi2_list))
    # print(sum(p_chi2_list)/len(p_chi2_list))
    # print(sum(g_list)/len(g_list))
    # print(sum(p_g_list)/len(p_g_list))
    # res_dict = {"chi2": chi2_list, "p-chi2": p_chi2_list, "g": g_list, "p-g": p_g_list}
    # df = pd.DataFrame(res_dict)
    # df.to_csv(obs_dir+"/res_true_st_true_gt.csv")

def plot():
    obs_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_gt/10000/"
    # obs_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_abcde/10000/"

    df = pd.read_csv(obs_dir+"res_true_st_true_gt.csv")
    fig = sns.histplot(data=df, x="chi2", bins=30, stat="probability", palette="Set2")
    x = np.arange(0, 1000, .0001)
    plt.plot(x, scipy.stats.chi2.pdf(x, df=DOF), color='r', lw=2)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel("$\chi^2$", fontsize=16)
    plt.ylabel("Proportion", fontsize=16)
    plt.tight_layout()
    plt.savefig(obs_dir+"/chi2_true_st_true_gt.pdf")

    plt.clf()
    fig = sns.histplot(data=df, x="p-chi2", bins=30, stat="probability", palette="Set2")
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel("p-value chi2", fontsize=16)
    plt.ylabel("Proportion", fontsize=16)
    plt.tight_layout()
    plt.savefig(obs_dir+"/pvalue-chi2_true_st_true_gt.pdf")


    plt.clf()
    fig = sns.histplot(data=df, x="g", bins=30, stat="probability", palette="Set2")
    x = np.arange(0, 1000, .0001)
    plt.plot(x, scipy.stats.chi2.pdf(x, df=DOF), color='r', lw=2)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel("g-test", fontsize=16)
    plt.ylabel("Proportion", fontsize=16)
    plt.tight_layout()
    plt.savefig(obs_dir+"/g-test_true_st_true_gt.pdf")


    plt.clf()
    fig = sns.histplot(data=df, x="p-g", bins=30, stat="probability", palette="Set2")
    x = np.arange(0, 1000, .0001)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel("p-value g", fontsize=16)
    plt.ylabel("Proportion", fontsize=16)
    plt.tight_layout()
    plt.savefig(obs_dir+"/pvalue-g-test_true_st_true_gt.pdf")

def plot_subplot(obs_dir):
    

    df = pd.read_csv(obs_dir+"res_true_st_true_gt.csv")
    print(df)
    fig, axes = plt.subplots(2, 2, figsize=(10,8))
    sns.histplot(data=df, x="chi2", stat="probability", bins=200, ax=axes[0][0])
    x = np.arange(0, 200, .005)
    axes[0][0].plot(x, scipy.stats.chi2.pdf(x, df=104), color='r', lw=2)
    sns.histplot(data=df, x="p-chi2", stat="probability", bins=20, ax=axes[0][1])
    sns.histplot(data=df, x="g", stat="probability",  bins=200, ax=axes[1][0])
    sns.histplot(data=df, x="p-g", stat="probability", bins=20, ax=axes[1][1])
    axes[0][0].set_xlim(0,200)
    axes[1][0].set_xlim(0,200)
    # axes[0][0].set_xlim(0,500)
    # axes[1][0].set_xlim(0,500)
    axes[0][1].set_xlim(0,1)
    axes[1][1].set_xlim(0,1)
    axes[1][0].plot(x, scipy.stats.chi2.pdf(x, df=104), color='r', lw=2)
    plt.savefig(obs_dir+"/test_null2.pdf")
    plt.show()

obs_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_gt/10000/"
NUM_GT=10000
# obs_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_abcde/10000/"
run_all(obs_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_gt/10000/",
    exp_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_gt/10000/gtprob_true_st.csv",
    outgroup = "Z 0",
    rate_dir="005"
    )
# plot()
# plot_subplot(obs_dir)