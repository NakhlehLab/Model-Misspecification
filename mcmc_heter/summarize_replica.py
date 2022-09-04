import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from ML_result_summarize import get_best_ML_net
import os.path
import MCMC_result_summarize
import sys


types = ["tree", "net", "AZ"]
types_map = {"tree":"Scenario A", "net": "Scenario B", "AZ": "Scenario C"}


def plot_replica_mcmc_net_hist_result(path1, path2, figpath):
    data1 = pd.read_csv(path1)
    data2 = pd.read_csv(path2)
    data = pd.concat([data1, data2])
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    print(data)
    print(data[data["type"] == "tree"])
    tree_data = data[(data["type"] == "tree") & (data["murate"] == True) & (data["homo"] == True)]

    print(tree_data)
    types = ["tree", "net"]
   

    diffs = []
    for index, row in data.iterrows():
        # print(index, row)
        if bool(row["murate"]) == True:       
            if bool(row["homo"]) == True:
                diffs.append(0)
            else:
                diffs.append(1)
        else:
            if bool(row["homo"]) == True:
                diffs.append(2)
            else:
                diffs.append(3)

    data["diff"] = diffs
    # distance
    for i in range(len(types)):
        for homo in [True, False]:
            for murate in [True, False]:
                print(data[(data["type"] == types[i]) &  (data["homo"] == homo) & (data["murate"] == murate)]["distance"].mean())
        
        sns.histplot(x='distance', hue='diff',  hue_order=[0,1,2,3],
                    multiple="dodge", discrete=True, binwidth=0.1,
                    data=data[data["type"] == types[i]], 
                    ax=axes[0][i], 
                    
                    )
        if i == 0:
            labels = axes[0][0].get_legend()
            handles = labels.legendHandles
        axes[0][i].set_xlabel("Distance", fontsize=16)
        axes[0][i].set_ylabel("", fontsize=20)
        print(axes[0][i].get_yticklabels())
        # axes[0][i].set_xticklabels([homo_map[t.get_text()] for t in axes[0][i].get_xticklabels()], size=20)
        axes[0][i].set_ylim(0, 11)
        axes[0][i].set_xlim(-1, 9)
        axes[0][i].set_title(types_map[types[i]], size=16, fontweight="bold")
        axes[0][i].get_legend().remove()
        axes[0][i].tick_params(axis='y', labelsize=14)
        axes[0][i].tick_params(axis='x', labelsize=14)
        print(types[i])
        
        # print(data[data["type"] == types[i] & data["murate"] == ""] )
    
    # reticulation
    for i in range(len(types)):
        sns.histplot(x='inferred_reti', hue='diff', hue_order=[0,1,2,3],
                     multiple="dodge", discrete=True, binwidth=0.1,
                    data=data[data["type"] == types[i]],
                    ax=axes[1][i], 
                    )
       

        axes[1][i].set_xlabel("Inferred # of reticulations", fontsize=16)
        axes[1][i].set_ylabel("", fontsize=20)
        # axes[1][i].set_xticklabels([homo_map[t.get_text()] for t in axes[1][i].get_xticklabels()], size=20)
        axes[1][i].set_ylim(0, 11)
        axes[1][i].set_xlim(-0.6, 5)
        # axes[1][i].set_title(types[i], size=20, fontweight="bold")
        axes[1][i].get_legend().remove()
        axes[1][i].tick_params(axis='y', labelsize=14)
        axes[1][i].tick_params(axis='x', labelsize=14)
    
    axes[0][0].set_ylabel("Count", fontsize=16)
    axes[1][0].set_ylabel("Count", fontsize=16)
    
    fig.legend(handles, ["DR inference & SR data","DR inference & DR data","SR inference & SR data","SR inference & DR data"], loc='lower center', fontsize=16, ncol=2)
    # plt.tight_layout()
    fig.subplots_adjust(bottom=0.15, hspace=0.18)
    plt.savefig(figpath, dpi=300)

def plot_replica_mcmc_net_hist_result1(path, figpath):
    data = pd.read_csv(path)
    fig, axes = plt.subplots(2, 1, figsize=(7, 10))
    
    print(data)
    print(data[data["type"] == "tree"])
    tree_data = data[(data["type"] == "tree") & (data["murate"] == True) & (data["homo"] == True)]

    print(tree_data)
    types = ["net"]
    # homo_map = {"False": "SR data", "True": "DR data"}
    # murate_map = {"False": "SR inference", "True": "DR inference"}

    diffs = []
    for index, row in data.iterrows():
        # print(index, row)
        if bool(row["murate"]) == True:       
            if bool(row["homo"]) == True:
                diffs.append(0)
            else:
                diffs.append(1)
        else:
            if bool(row["homo"]) == True:
                diffs.append(2)
            else:
                diffs.append(3)

    data["diff"] = diffs
    # distance
    
    i = 0
    for homo in [True, False]:
        for murate in [True, False]:
            print(data[(data["type"] == types[i]) &  (data["homo"] == homo) & (data["murate"] == murate)]["distance"].mean())
    
    sns.histplot(x='distance', hue='diff',  hue_order=[0,1,2,3],
                multiple="dodge", discrete=True, binwidth=0.1,
                data=data[data["type"] == types[i]], 
                ax=axes[0], 
                
                )
    if i == 0:
        labels = axes[0].get_legend()
        handles = labels.legendHandles
    axes[0].set_xlabel("Distance", fontsize=16)
    axes[0].set_ylabel("", fontsize=20)
    # print(axes[0][i].get_yticklabels())
    # axes[0][i].set_xticklabels([homo_map[t.get_text()] for t in axes[0][i].get_xticklabels()], size=20)
    axes[0].set_ylim(0, 11)
    axes[0].set_xlim(-1, 9)
    axes[0].set_title(types_map[types[i]], size=16, fontweight="bold")
    axes[0].get_legend().remove()
    axes[0].tick_params(axis='y', labelsize=14)
    axes[0].tick_params(axis='x', labelsize=14)
    print(types[i])
        
        # print(data[data["type"] == types[i] & data["murate"] == ""] )
    
    # reticulation
    for i in range(len(types)):
        sns.histplot(x='inferred_reti', hue='diff', hue_order=[0,1,2,3],
                     multiple="dodge", discrete=True, binwidth=0.1,
                    data=data[data["type"] == types[i]],
                    ax=axes[1], 
                    )
       

        axes[1].set_xlabel("Inferred # of reticulations", fontsize=16)
        axes[1].set_ylabel("", fontsize=20)
        # axes[1][i].set_xticklabels([homo_map[t.get_text()] for t in axes[1][i].get_xticklabels()], size=20)
        axes[1].set_ylim(0, 11)
        axes[1].set_xlim(-0.6, 5)
        # axes[1][i].set_title(types[i], size=20, fontweight="bold")
        axes[1].get_legend().remove()
        axes[1].tick_params(axis='y', labelsize=14)
        axes[1].tick_params(axis='x', labelsize=14)
    
    axes[0].set_ylabel("Count", fontsize=16)
    axes[1].set_ylabel("Count", fontsize=16)
    
    fig.legend(handles, ["DR inference & SR data","DR inference & DR data","SR inference & SR data","SR inference & DR data"], loc='lower center', fontsize=16, ncol=2)
    # plt.tight_layout()
    fig.subplots_adjust(bottom=0.15, hspace=0.18)
    plt.savefig(figpath, dpi=300)


def summarize_iqtree_error_to_csv(mcmc_dir, locus_length, type_net, csv_path):
    # types = ["tree", "net"]
    homos = ["heter", "homo"]
    data = {"type": [], "homo": [], "index": [], "RF": [], "nrBS": [], "right_gt":[]}
    # for type in types:
    for i in range(1, 11):
        for h in homos:
            path = mcmc_dir + "/" + str(i) + "/"+h+"/"+str(locus_length) + "/iqtree_err.csv"
            df = pd.read_csv(path)
            data["nrBS"].append(np.mean(df["nrBS"]))
            data["RF"].append(np.mean(df["RF-out"]))
            data["type"].append(type_net)
            data["homo"].append(h == "homo")
            data["index"].append(i)
            data["right_gt"].append((df['RF-in'] == 0).sum())

    print(data)
    df = pd.DataFrame.from_dict(data)
    df.to_csv(csv_path)

def summarize_replica_ML_result_to_csv(dir, type, locus_length, csv_path):
    homos = ["heter", "homo"]
    gt_types = ["true", "iq"]
    # gt_types = ["iqtree"]
    data = {"type": [], "index": [], "homo": [], "gt": [], "reti_set": [], "net": [], "likelihood": []}
    for i in range(1, 11):
        for h in homos:
            for reti in range(4):
                for gtt in gt_types:
                    if gtt == "true" and h == "heter":
                        continue
                    path = dir + "/" + str(i) + "/" + h + "/"+ str(locus_length) + "/ML_"+gtt+"_"+str(reti)+".out"
                    if not os.path.exists(path):
                        print(path)
                        continue
                    net, score = get_best_ML_net(path)
                    if net is None or score is None:
                        print(path)

                    data["type"].append(type)
                    data["index"].append(i)
                    data["homo"].append((h == "homo"))
                    data["gt"].append((gtt=="true"))
                    data["reti_set"].append(reti)
                    data["net"].append(net)
                    data["likelihood"].append(score)

    df = pd.DataFrame.from_dict(data)
    df.to_csv(csv_path)


def plot_ML_ll_trend_difference(csv_path1, csv_path2, figpath):
    pd.set_option("display.max_rows", None, "display.max_columns", None)
    data = pd.concat([pd.read_csv(csv_path1), pd.read_csv(csv_path2)])
    types = ["tree", "net"]
    diff_type = 1
    homos_map = {True: "SR data", False: "DR data"}
    gtt_map = {True: "True", False: "IQTREE"}

    fig, axes = plt.subplots(3, 2, figsize=(14, 10))
    
    
    diffs = []
    for index, row in data.iterrows():
        # print(index, row)
        if int(row["reti_set"]) == 0:       
            if diff_type == 1:
                loc0 = row["likelihood"]
                diffs.append(0)
            else:
                diffs.append(None)
        elif int(row["reti_set"]) in [1,2,3]:
            if diff_type == 1:
                diffs.append(row["likelihood"] - loc0)
            else:
                diffs.append(row["likelihood"] - data.loc[index-1]["likelihood"])
    data["diff"] = diffs

    print(data[["likelihood", "diff"]])
    for i in range(len(types)):
        df = data[data["type"] == types[i]]
        j = 0
        for homo in [True, False]:
            for gtt in [True, False]:
                if homo == False and gtt == True:
                    continue
                subplot_data = df[(df["homo"] == homo) & (df["gtt"] == gtt)]
                sns.lineplot(data=subplot_data, x="reti_set", y="diff", hue="index", ax=axes[j][i], alpha=0.5, linewidth=1.0)
                sns.lineplot(data=subplot_data, x="reti_set", y="diff", ax=axes[j][i], ci=None, linewidth=2)

                axes[j][i].set_xlabel("", fontsize=20)
                axes[j][i].set_ylabel("", fontsize=16)
                # print(subplot_data["likelihood"].min()-1)
                if i == 0:
                    if homo == True and gtt == True:
                        axes[j][i].set_ylabel("True GTs", size=16)
                    elif homo == True and gtt == False:    
                        axes[j][i].set_ylabel("Inferred SR GTs", size=16)
                    elif homo == False and gtt == False:    
                        axes[j][i].set_ylabel("Inferred DR  GTs", size=16)
                # handles, labels = axes[i][j].get_legend_handles_labels()
                axes[j][i].get_legend().remove()
                axes[j][i].set_xticks(range(4))
                lower_value = subplot_data["likelihood"].min() - 1
                # print(y_range[i])
                axes[j][i].set_ylim(-1, 25)
                axes[j][i].set_xlim(0, 3.1)
                axes[j][i].tick_params(axis='y', labelsize=14)
                axes[j][i].tick_params(axis='x', labelsize=14)
                j += 1

            # new_labels = [murate_map[x] for x in labels[0:2]]
        axes[0][i].set_title(types_map[types[i]], size=16,  fontweight="bold")

    
    fig.text(0.5, 0.04, 'Set # of reticulations', ha='center', fontsize=18, fontweight="bold")
    fig.text(0.04, 0.5, 'Relative log-likelihood', va='center', rotation='vertical', fontsize=18, fontweight="bold")
    # plt.tight_layout(pad=1.5)
    handles, labels = axes[2][1].get_legend_handles_labels()
    # new_labels = [murate_map[x] for x in labels[0:2]]

    print(handles, labels)

    # fig.legend(handles, labels, loc='lower center', fontsize=20, ncol=5)
    plt.tight_layout(rect=[0, 5, 1, 1])
    plt.savefig(figpath, dpi=300)
    plt.show()


def plot_ML_distance_hist(csv_path1, csv_path2, figpath):
    data = pd.concat([pd.read_csv(csv_path1, dtype={"distance":np.float64}),pd.read_csv(csv_path2, dtype={"distance":np.float64})])
    types = ["tree", "net"]
    pd.set_option("display.max_rows", None, "display.max_columns", None)
    homos = [True, False]
    list_labels = [(True, True), (False, True), (False, False)]

    diffs = []
    for index, row in data.iterrows():
        # print(index, row)
        if bool(row["gtt"]) == True:       
            if bool(row["homo"]) == True:
                diffs.append(0)
            else:
                diffs.append(1)
        else:
            if bool(row["homo"]) == True:
                diffs.append(2)
            else:
                diffs.append(3)

    data["diff"] = diffs

    data = data[data['diff'] != 1]

    f, axes = plt.subplots(1, 2, figsize=(14, 5))
    reti_map = {"tree": 0, "net": 1}
    for i in range(len(types)):
        print(types[i])
        sns.histplot(x='distance', hue="diff", hue_order=[0,2,3], multiple="dodge", discrete=True, binwidth=0.1,
                    data=data[(data["type"] == types[i]) & (data["reti_set"] == reti_map[types[i]])],
                    ax=axes[i], bins=5,
                    )
        if i==0:
            labels = axes[0].get_legend()
            handles = labels.legendHandles
        axes[i].set_xlabel("Distance", fontsize=16)
        axes[i].set_ylabel("", fontsize=16)
        # axes[i].set_xticklabels([homo_map[t.get_text()] for t in axes[i].get_xticklabels()], size=20)
        axes[i].set_xlim(-1, 7)
        axes[i].set_ylim(0, 11)
        axes[i].get_legend().remove()
        # if j == 1:
        #     axes[i].set_xlabel("Distance", fontsize=20)
        # else:
        axes[i].set_title(types_map[types[i]], size=16, fontweight="bold")
        axes[i].tick_params(axis='y', labelsize=14)
        axes[i].tick_params(axis='x', labelsize=14)

    
    labels.remove()
    axes[0].set_ylabel("Count", fontsize=16)
    # axes[1][0].set_ylabel("Hetero", fontsize=20)

    f.legend(handles, ["True GTs", "Inferred SR GTs", "Inferred DR GTs"], loc='lower center', fontsize=16, ncol=3, )
    f.subplots_adjust(bottom=0.28, wspace = 0.15)
    plt.tight_layout(rect=[0, 5, 1, 1])
    plt.savefig(figpath, dpi=300)


def compute_mcmc_gt_error(mcmc_dir, locus_length, bl, sf, burnin=0):
    # types = ["tree", "net"]
    homos = ["heter", "homo"]
    murates = [True, False]
    
    for j in range(1, 11):
        for h in homos:
            for mu in murates:
                path = mcmc_dir + "/" + str(j) + "/"+h+"/"+str(locus_length) + "/"
                true_tree_path = mcmc_dir + "/" + str(j) + "/homo/"+str(locus_length)+"/genetrees.txt"
                if mu:
                    path += "murate/"
                print(path)
                print(true_tree_path)
                MCMC_result_summarize.trace_mcmc_gene_tree_errors(path, true_tree_path, bl, sf, burnin=burnin)

def summarize_mcmc_gt_error_to_csv(mcmc_dir, locus_length, type_net, csv_path):
    # types = ["tree", "net"]
    homos = ["heter", "homo"]
    murates_map = {"murate": True, "": False}
    data = {"type": [], "homo": [], "murate": [], "index": [], "RF": [], "nrBS": []}
    
    for h in homos:
        for mu in murates_map.keys():
            for j in range(1, 11):
                path = mcmc_dir + "/" + str(j) + "/"+h+"/"+str(locus_length) + "/"
                if mu:
                    path += "/murate"
                path += "/mcmc_gt_err.txt"
                print(path)
                if not os.path.exists(path):
                    continue
                print(path)
                with open(path, "r") as handle:
                    lines = handle.read().split("\n")
                    data["type"].append(type_net)
                    data["homo"].append(h == "homo")
                    data["index"].append(j)
                    data["murate"].append(murates_map[mu])

                    for line in lines:
                        if "nrBS" in line:
                            data["nrBS"].append(line.strip().split(":")[1].strip())
                        elif "RF" in line:
                            data["RF"].append(line.strip().split(":")[1].strip())

                    handle.close()
    print(data)
    for k in data.keys():
        print(len(data[k]))
    df = pd.DataFrame.from_dict(data)
    df.to_csv(csv_path)


def plot_gt_error_mcmc_iqtree_subfig(csv_mcmc1, csv_iq1, csv_mcmc2, csv_iq2, figpath, metric):
    data = pd.concat([pd.read_csv(csv_mcmc1), pd.read_csv(csv_mcmc2)])
    data_iq = pd.concat([pd.read_csv(csv_iq1), pd.read_csv(csv_iq2)])
    data = data.append(data_iq, ignore_index = True)
    data['murate'] = data['murate'].replace(True, 1).replace(False, 2).fillna(3)
    data = data.astype({'murate': 'int32'})
    # print(data)

    f, axes = plt.subplots(1, 2, figsize=(16, 6))

    types = ["tree", "net", "AZ"]
    # metric = "nrBS"
    # metric = "RF"
    homo = [False, True]
    homo_map = {"False": "DR data", "True": "SR data"}
    murate_map = {'2': "SR inference", '1': "DR inference", '3': "IQ-TREE"}
    print(data)
    
    for j in range(len(homo)):
        for murate in [1, 2, 3]:
            for k in range(len(types)):
                print(types[k], homo[j], murate_map[str(murate)])
                # print(data[(data["type"] == types[k]) & (data["homo"] == homo[j]) &  (data["murate"] == murate)])
                print(data[(data["type"] == types[k]) & (data["homo"] == homo[j]) &  (data["murate"] == murate)][metric].mean())
            print("------------------------------")
            print(homo[j], murate_map[str(murate)],metric)
            print(data[(data["homo"] == homo[j]) &  (data["murate"] == murate)][metric].mean())

        sns.boxplot(x='type', y=metric, hue='murate', 
                    data=data[data["homo"] == homo[j]], 
                    width=0.2, orient='v',
                    ax=axes[j], palette="Set2",
                    )

        axes[j].set_xlabel("", fontsize=16)
        if metric == "RF":
            axes[j].set_ylabel("nrRF", fontsize=16)
        else:
            axes[j].set_ylabel("nrBS", fontsize=16)
        axes[j].tick_params(axis='y', labelsize=14)
        axes[j].tick_params(axis='x', labelsize=14)
        axes[j].set_xticklabels([types_map[t.get_text()] for t in axes[j].get_xticklabels()], size=16)
        
        axes[j].set_title(homo_map[str(homo[j])], size=16, fontweight="bold")
        if metric == "RF":
            axes[j].set_ylim(0, 0.15)
        else:
            axes[j].set_ylim(0.3, 0.5)
        
        axes[j].get_legend().remove()
    handles, labels = axes[0].get_legend_handles_labels()

    new_labels = [murate_map[x] for x in labels[0:3]]
    f.legend(handles, new_labels, loc='lower center', fontsize=16, ncol=3, )
    # f.set_tight_layout(True)
    f.subplots_adjust(bottom=0.2, wspace = 0.2)
    plt.savefig(figpath, dpi=300)


def summarize_mcmc_murate_error_to_csv(mcmc_dir, locus_length, net_type, csv_path, burnin=0):
    # types = ["tree", "net"]
    data = {"type": [], "index": [], "loci":[], "true_mu":[], "inferred_mu":[]}

    for j in range(1, 11):
        path = mcmc_dir + "/" + str(j) + "/heter/"+str(locus_length)+"/murate/"
        true_mu_path = mcmc_dir + "/" + str(j) + "/heter/"+str(locus_length)+"/murate.txt"
        true_mu_rates = MCMC_result_summarize.read_true_mu_path(true_mu_path)
        locus2rates = MCMC_result_summarize.read_mcmc_mu_rates(path)
        start = int(burnin / 100 * len(locus2rates["loci1"]))
        print(path)
        for k in range(len(true_mu_rates)):
            if len(locus2rates["loci1"]) == 0:
                continue
            data["type"].append(net_type)
            data["index"].append(j)
            data["loci"].append(k+1)
            data["true_mu"].append(true_mu_rates[k])
            data["inferred_mu"].append(
                sum(locus2rates["loci" + str(k+1)][start:]) / len(locus2rates["loci" + str(k+1)][start:]))

    df = pd.DataFrame.from_dict(data)
    df.to_csv(csv_path)


def plot_mcmc_murates(csv_path1, csv_path2, fig_path):
    data = pd.concat([pd.read_csv(csv_path1), pd.read_csv(csv_path2)])
    label_map = {"tree":"Scenario A", "net":"Scenario B"}
    fig, ax = plt.subplots()
    sns.scatterplot(data=data, x="true_mu", y="inferred_mu", hue="type", s=10, ax=ax, palette="Set2")
    sns.lineplot(data=data, x="true_mu", y="true_mu", ax=ax, color="blue")
    ax.tick_params(axis='y', labelsize=14)
    ax.tick_params(axis='x', labelsize=14)
    ax.set_xlabel("True relative substitution rate", fontsize=16)
    ax.set_ylabel("Inferred relative substitution rate", fontsize=16)
    handles, labels = ax.get_legend_handles_labels()
    new_labels = [label_map[x] for x in labels[0:3]]
    ax.legend(handles, new_labels, loc='upper left', fontsize=16)
    # ax.set(xscale="log", yscale="log")
    errors = abs(data["true_mu"] - data["inferred_mu"])/data["true_mu"]

    print(errors.mean())

    plt.savefig(fig_path, dpi=300)

def plot_mcmc_ess(csv_path, fig_path):
    data = pd.read_csv(csv_path)
    fig, ax = plt.subplots()
    sns.scatterplot(data=data, x="ESS", y="distance", hue="type", s=15, ax=ax, palette="Set2")
    ax.tick_params(axis='y', labelsize=12)
    ax.tick_params(axis='x', labelsize=12)
    ax.set_xlabel("ESS", fontsize=16)
    ax.set_ylabel("Distance", fontsize=16)
    handles, labels = ax.get_legend_handles_labels()
    ax.set(xscale="log")
    ax.legend(handles, labels, loc='upper left', fontsize=14)

    plt.savefig(fig_path, dpi=300)
    for index, row in data.iterrows():
        if row["ESS"] <100 and row["distance"] > 0 and row["murate"] == True:
            print(row["type"]+","+str(row["index"])+","+str(row["homo"])+","+str(row["murate"])+","+str(row["num"])+","+str(row["ESS"]))


if __name__ == '__main__':
    bl = 2000000
    sf = 5000
    burnin=50

    fig_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/"
    # path1 = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/net0/long/mcmc_result.csv"
    # path2 = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/net1_2/long/mcmc_result.csv"
    ifmap="MAP"
    scale = "long"
    path1 = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/net0/"+scale+"/mcmc_result"+ifmap+".csv"
    path2 = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/net1_2/"+scale+"/mcmc_result"+ifmap+".csv"
    figpath_mcmc = fig_dir+"replica_mcmc_"+scale+"_"+ifmap+".pdf"
    plot_replica_mcmc_net_hist_result(path1, path2, figpath_mcmc)
    

    experiment_dir1 = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/net0/long/"
    experiment_dir2 = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/net1_2/long/"
    
    # experiment_dir1 = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/net0/medium/"
    # experiment_dir2 = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/net1_2/medium/"
    mcmc_murate_csv1 = experiment_dir1+"MCMC_murates.csv"
    mcmc_murate_csv2 = experiment_dir2+"MCMC_murates.csv"
    summarize_mcmc_murate_error_to_csv(experiment_dir1, 2000, "tree", mcmc_murate_csv1, burnin=burnin)
    summarize_mcmc_murate_error_to_csv(experiment_dir2, 2000, "net", mcmc_murate_csv2, burnin=burnin)
    mcmc_murate_fig = fig_dir+"/MCMC_murates_medium.pdf"
    plot_mcmc_murates(mcmc_murate_csv1, mcmc_murate_csv2, mcmc_murate_fig)  

    
    ML_result_csv1 = experiment_dir1+"/ML_result.csv"
    ML_result_csv2 = experiment_dir2+"/ML_result.csv"
   

    ML_fig_trend_diff_path = fig_dir + "replica_ML_trend_diff.pdf"
    ML_dist_hist_fig_path = fig_dir + "ML_distance_hist.pdf"
    plot_ML_ll_trend_difference(ML_result_csv1, ML_result_csv2, ML_fig_trend_diff_path)
    plot_ML_distance_hist(ML_result_csv1, ML_result_csv2, ML_dist_hist_fig_path)

    gt_error_csv1 = experiment_dir1 + "/MCMC_gt_err.csv"
   
    gt_err_fig_rf = fig_dir + "/replica_gt_err_rf.pdf"
    gt_err_fig_nrbs = fig_dir + "/replica_gt_err_nrbs.pdf"

    summarize_mcmc_gt_error_to_csv(sys.argv[1], 2000, "tree", sys.argv[1] + "/MCMC_gt_err.csv")
    summarize_mcmc_gt_error_to_csv(sys.argv[2], 2000, "net", sys.argv[2] + "/MCMC_gt_err.csv")
    summarize_iqtree_error_to_csv(sys.argv[1], 2000, "tree", sys.argv[1] + "/iq_gt_err.csv")
    summarize_iqtree_error_to_csv(sys.argv[2], 2000, "net", sys.argv[2] + "/iq_gt_err.csv")
    
    mcmc_gt_csv1 = experiment_dir1+"/MCMC_gt_err.csv"
    mcmc_gt_csv2 = experiment_dir2+"/MCMC_gt_err.csv"
    iqtree_err_csv1 = experiment_dir1+"/iq_gt_err.csv"
    iqtree_err_csv2 = experiment_dir2+"/iq_gt_err.csv"

    plot_gt_error_mcmc_iqtree_subfig(mcmc_gt_csv1, iqtree_err_csv1, mcmc_gt_csv2, iqtree_err_csv2, gt_err_fig_rf, "RF")
    plot_gt_error_mcmc_iqtree_subfig(mcmc_gt_csv1, iqtree_err_csv1, mcmc_gt_csv2, iqtree_err_csv2, gt_err_fig_nrbs, "nrBS")
