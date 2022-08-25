from tkinter import FALSE
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import dendropy
from dendropy.calculate import treecompare
from pathlib import Path
import sys


gt_name = {True: "truegt", False: "infergt"}

def plot_branches(brl_path, truegt):
    path = Path(brl_path)
    folder = path.parent.absolute()
    df = pd.read_csv(brl_path, index_col=False)
    # print(df)
    data = df.to_dict()
    print(data)
    for k in data.keys():
        k_arr = k.split("-")
        # print(k_arr)
        if k_arr == None or len(k_arr) < 2:
            continue
        # print(k_arr)
        if k_arr[1].startswith("I"):
            branch_list = list(data[k].values())
            sns.scatterplot(range(len(branch_list)),branch_list)
            plt.savefig(str(folder)+"/brl_"+gt_name[truegt]+"_"+k+".pdf")
            plt.show()
            print(branch_list)

def plot_branches_box(brl_path, truegt):
    path = Path(brl_path)
    folder = path.parent.absolute()
    df = pd.read_csv(brl_path, index_col=False, usecols=["I0-I1","I0-I2", "I1-I3"])
    # print(df)
    df.boxplot()
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.savefig(str(folder)+"/brl_"+gt_name[truegt]+"_box"+".pdf")
    # plt.show()
    # data = df.to_dict()
    # print(data)
    # for k in data.keys():
    #     k_arr = k.split("-")
    #     # print(k_arr)
    #     if k_arr == None or len(k_arr) < 2:
    #         continue
    #     # print(k_arr)
    #     if k_arr[1].startswith("I"):
    #         branch_list = list(data[k].values())
    #         sns.scatterplot(range(len(branch_list)),branch_list)
    #         plt.savefig(str(folder)+"/brl_"+gt_name[truegt]+"_"+k+".pdf")
    #         plt.show()
    #         print(branch_list)

if __name__ == "__main__":
    directory = sys.argv[1]
    truegt = (sys.argv[2].lower() == "true")
    print(truegt)
    brl_path = directory+"/st_brl_"+gt_name[truegt]+".csv"
    plot_branches_box(brl_path, truegt)

    # truegt = True
    # brl_path = directory+"/st_brl_"+gt_name[truegt]+".csv"
    # plot_branches_box(brl_path, truegt)
