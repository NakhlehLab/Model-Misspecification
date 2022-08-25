import os
import matplotlib.pyplot as plt
from scipy.stats import bayes_mvs
import sys
import pandas as pd
# from dendropy.calculate import compare
# import dendropy


CHAIN_LEN = 20000000
BURN_IN = 2000000
SAMPLE_FREQ = 5000
NUM_SEED = 2
NUM_RUN = 8
NUM_REPLICATE = 10
# TRUE_NET = dendropy.Tree.get(data="((((Q:0.4)#H1:0.2::0.7,R:0.6)I3:1.0,(L:0.8,#H1:0.4::0.3):0.8)I1:0.4,C:2);", schema="newick", rooting="force-rooted")

def mcmc_posterior(mcmc_path, burn_in_proportion=0):
    with open(mcmc_path, "r") as handle:
        begin = False
        posterior_list = []
        max_posterior = float("-inf")
        max_posterior_net = None
        lastESS = 0
        lines = handle.read().split("\n")
        for index, line in enumerate(lines):
            if line.strip().startswith("Iteration;"):
                begin = True
            elif begin:
                if line.startswith("-"):
                    break
                elif not line.startswith("["):
                    arr = line.strip().split()
                    posterior_prob = float(arr[1][:-1])
                    lastESS = float(arr[2][:-1])
                    if max_posterior < posterior_prob:
                        max_posterior = posterior_prob
                        max_posterior_net = lines[index+1]
                    posterior_list.append(posterior_prob)
        
        return posterior_list, max_posterior, max_posterior_net, lastESS


def mcmc_trace_chains(dir_path, seed_index):
    posterior_list = []
    net_dict = {"seed":[], "num_run":[], "posterior_prob":[], "net":[], "lastESS":[]}
    for i in range(NUM_RUN):
        sub_dir = dir_path + "/mcmcseq/0_"+str(seed_index) + "/"+str(CHAIN_LEN)+"_"+str(i)+"/"
        if not os.path.exists(sub_dir):
            break
        
        file_names = os.listdir(sub_dir)
        slurm_name = None
        for f in file_names:
            if f.startswith("slurm"):
                if slurm_name is None:
                    slurm_name = f
                elif slurm_name < f:
                    slurm_name = f
        if slurm_name != None:
            mcmc_path = sub_dir + slurm_name
            print(mcmc_path)
            posteriors, max_posterior, max_posterior_net, lastESS = mcmc_posterior(mcmc_path)
            # with open(sub_dir+"check.txt", "r") as handle:
            #     lastESS = float(handle.read().split("\n")[3].strip())
            if i == 0:
                posterior_list.extend(posteriors[BURN_IN//SAMPLE_FREQ + 1 : ])
            # print(i)
            net_dict["seed"].append(seed_index)
            net_dict["num_run"].append(i)
            net_dict["posterior_prob"].append(max_posterior)
            net_dict["net"].append(max_posterior_net)
            # dendropy.Tree.get
            net_dict["lastESS"].append(lastESS)
            # print(net_dict)
        else:
            break
    return posterior_list, pd.DataFrame(net_dict)
    # f, axes = plt.subplots(figsize=(10, 8))
    # plt.plot(range(len(posterior_list)), posterior_list, c="black")
    # plt.ylim(-94008, )
    # plt.savefig()
    # plt.show()
    
def mcmc_check_trace_replicate(dir_path, seed_index, locus_length, murate, heter):
    print(dir_path, seed_index, murate, heter)
    f, axes = plt.subplots(2, 5, figsize=(20, 8))
    df_list = []
    for rep in range(NUM_REPLICATE):
        if heter:
            mcmc_dir = dir_path + str(rep+1)+"/heter/"+locus_length+"/"
            fig_path = dir_path+"/posterior_prob_trace_heter"
        else:
            mcmc_dir = dir_path + str(rep+1)+"/homo/"+locus_length+"/"
            fig_path = dir_path+"/posterior_prob_trace_homo"

        if murate:
            mcmc_dir += "murate/"
            fig_path += "_murate"
        # print(mcmc_dir, fig_path, rep)
        posterior_list, net_df = mcmc_trace_chains(mcmc_dir, seed_index)
        net_df["replicate"] = [rep for x in range(len(net_df["seed"]))]
        df_list.append(net_df)
        if rep%2 == 0:
            axes[0][rep//2].plot(range(len(posterior_list)), posterior_list, c="black")
            # axes[0][rep//2].set_ylim(-99900, )
            # axes[0][rep//2].set_ylim(-95000, )

        else:
            axes[1][rep//2].plot(range(len(posterior_list)), posterior_list, c="black")
            # axes[1][rep//2].set_ylim(-99900, )
            # axes[0][rep//2].set_ylim(-95000, )

    # plt.ylim(-100000, )
    
    res_df = pd.concat(df_list)
    res_df.to_csv(fig_path + ".csv")
    plt.savefig(fig_path+".pdf")



if __name__ == '__main__':
    
    # dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/long/1/heter/murate/"
    # seed_index = 0
    # mcmc_trace_chains(dir_path, seed_index)
   
    mcmc_check_trace_replicate(sys.argv[1], int(sys.argv[2]), sys.argv[3], (sys.argv[4].lower() == "true"), (sys.argv[5].lower() == "true"))