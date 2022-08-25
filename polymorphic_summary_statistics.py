
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statistics
import sys


def plot_box_polymorphic(dir_path):
    # dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/"
    # locus_length = 500
    plt.clf()
    df_list = []
    for scale in ["short", "medium", "long"]:
        for locus_length in [200, 1000]:
            csv_path = dir_path+scale+"/polymorphic_heter_"+str(locus_length)+".csv"
            df = pd.read_csv(csv_path)
            df["scale"] = [scale for x in range(len(df["replicate"]))]
            df["locus_length"] = [locus_length for x in range(len(df["replicate"]))]
            df_list.append(df)
            mean = sum(df["diff"])/len(df["diff"])
            std = statistics.pstdev(df["diff"])
            print(scale, mean, std)
    res_df = pd.concat(df_list, ignore_index=True)
    print(res_df)
    res_df.to_csv(dir_path+"polymorphic_"+".csv")
    ax = sns.boxplot(x="scale", y="diff", data=res_df, hue="locus_length", showfliers = False)
    ax.set_xticklabels(ax.get_xticklabels(),rotation=30)
    plt.savefig(dir_path+"polymorphic_boxplot.pdf")

def plot_box_pdistance(dir_path):
    # dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/"
    # locus_length = 500
    
    for locus_length in [200, 1000]:
        plt.clf()
        df_list = []
        for scale in ["short", "medium", "long"]:
            csv_path = dir_path+scale+"/pdist_heter_"+str(locus_length)+".csv"
            df = pd.read_csv(csv_path)
            df["scale"] = [scale for x in range(len(df["replicate"]))]
            df_list.append(df)
        res_df = pd.concat(df_list, ignore_index=True)
        print(res_df)
        res_df.to_csv(dir_path+"pdist_"+str(locus_length)+".csv")
        ax = sns.boxplot(x="species_i_j", y="p-distance", data=res_df, hue="scale", showfliers = False)
        ax.set_xticklabels(ax.get_xticklabels(),rotation=30)
        plt.savefig(dir_path+"pdist_"+str(locus_length)+"_boxplot.pdf", dpi=300)

def plot_box_gt_error(dir_path, ingroup=False):
    # dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/"
    df_list = []
    plt.clf()
    for scale in ["short", "medium", "long"]:
        for locus_length in [200, 1000]:
            csv_path = dir_path+scale+"/iqtree_error"+str(locus_length)+".csv"

            df = pd.read_csv(csv_path)
            df["scale"] = [scale for x in range(len(df["index"]))]
            df["locus_length"] = [locus_length for x in range(len(df["index"]))]
            df_list.append(df)
        # mean = sum(df["diff"])/len(df["diff"])
        # std = statistics.pstdev(df["diff"])
        # print(scale, mean, std)
    res_df = pd.concat(df_list, ignore_index=True)
    print(res_df)
    if ingroup:
        # res_df.to_csv(dir_path+"iqtree_error_ingroup.csv")
        res_df.to_csv(dir_path+"iqtree_error.csv")
        ax = sns.boxplot(x="scale", y="RF-in", data=res_df, hue="locus_length", showfliers = False)
        plt.savefig(dir_path+"iqtree_error_ingroup"+"_boxplot.pdf")
    else:
        res_df.to_csv(dir_path+"iqtree_error.csv")        
        ax = sns.boxplot(x="scale", y="RF-out", data=res_df, hue="locus_length", showfliers = False)
        plt.savefig(dir_path+"iqtree_error_boxplot.pdf")
    
def plot():
    dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/test_outgroup/iqtree_error.csv"
    res_df = pd.read_csv(dir_path)
    ax = sns.boxplot(x="scale", y="RF-in", data=res_df, hue="locus_length", showfliers = False)
    plt.savefig(dir_path+"iqtree_error_ingroup_boxplot.pdf")

def plot_gt_all(dir_path, re_end, unroot=True, ingroup=False, longer=True):
    df_list = []
    if longer:
        locus_length_list = [500, 2000]
    else:
        locus_length_list = [200, 1000]
    
    for scale in ["short", "medium", "long"]:
        for locus_length in locus_length_list:
        # for locus_length in [200, 1000]:
            for i in range(1, re_end+1):
                csv_path = dir_path+scale+"/"+str(i)+"/heter/"+str(locus_length)+"/gt_stat.csv"
                df = pd.read_csv(csv_path)
                df["index"] = [i for x in range(len(df["RF-out"]))]
                df["scale"] = [scale for x in range(len(df["RF-out"]))]
                df["locus_length"] = [locus_length for x in range(len(df["RF-out"]))]
                df_list.append(df)
        # mean = sum(df["diff"])/len(df["diff"])
        # std = statistics.pstdev(df["diff"])
        # print(scale, mean, std)
    res_df = pd.concat(df_list, ignore_index=True)
    res_df.to_csv(dir_path+"gt_stat.csv")
    
    plt.clf()
    if unroot:
        fig, axes = plt.subplots(1, 2, figsize=(20,5))
        ax = sns.boxplot(x="scale", y="RF-in_unroot", data=res_df, hue="locus_length", showmeans=True, showfliers = False)
        for idx, locus_length in enumerate(locus_length_list):
            df_locus = res_df[res_df["locus_length"] == locus_length]
            counted_data = df_locus.groupby(["scale"])["RF-in_unroot"].value_counts(normalize=True, sort=False)
            dict_res=[{"RF category":scale, 'scale':count, 'percent (%)':percent*100} for (count, scale), percent in dict(counted_data).items()]
            df_bar = pd.DataFrame(dict_res)
            sns.barplot(x="RF category", y="percent (%)", hue="scale", data=df_bar, ax = axes[idx])
            axes[idx].set_title(label="locus_length="+str(locus_length))
            axes[idx].set_ylim(0,100)
        
        # plt.grid(axis="y")
        plt.savefig(dir_path+"gt_stat_ingroup_unroot_count.pdf")
        # ax = sns.boxplot(x="scale", y = "RF-in_unroot", data=res_df, hue="locus_length", showmeans=True, showfliers = False)
        # plt.savefig(dir_path+"gt_stat_ingroup_unroot_count_box.pdf")

    elif ingroup:
        fig, axes = plt.subplots(1, 2, figsize=(20,5))
        # ax = sns.boxplot(x="scale", y="RF-in", data=res_df, hue="locus_length", showmeans=True, showfliers = False)
        for idx, locus_length in enumerate(locus_length_list):
            df_locus = res_df[res_df["locus_length"] == locus_length]
            counted_data = df_locus.groupby(["scale"])["RF-in"].value_counts(normalize=True, sort=False)
            dict_res=[{"RF category":scale, 'scale':count, 'percent (%)':percent*100} for (count, scale), percent in dict(counted_data).items()]
            df_bar = pd.DataFrame(dict_res)
            sns.barplot(x="RF category", y="percent (%)", hue="scale", data=df_bar, ax = axes[idx])
            axes[idx].set_title(label="locus_length="+str(locus_length))
            axes[idx].set_ylim(0,100)
            axes[idx].grid(axis="y")
        plt.savefig(dir_path+"gt_stat_ingroup_count.pdf")
        
        # ax = sns.boxplot(x="scale", y="RF-in", data=res_df, hue="locus_length", showmeans=True, showfliers = False)
        # plt.savefig(dir_path+"gt_stat_ingroup_box.pdf")
    
    else:
        
        fig, axes = plt.subplots(1, 2, figsize=(20,5))
        # ax = sns.boxplot(x="scale", y="RF-out", data=res_df, hue="locus_length", showmeans=True, showfliers = False)
        for idx, locus_length in enumerate(locus_length_list):
            df_locus = res_df[res_df["locus_length"] == locus_length]
            counted_data = df_locus.groupby(["scale"])["RF-out"].value_counts(normalize=True, sort=False)
            dict_res=[{"RF category":scale, 'scale':count, 'percent (%)':percent*100} for (count, scale), percent in dict(counted_data).items()]
            df_bar = pd.DataFrame(dict_res)
            sns.barplot(x="RF category", y="percent (%)", hue="scale", data=df_bar, ax = axes[idx])
            axes[idx].set_title(label="locus_length="+str(locus_length))
            axes[idx].set_ylim(0,100)
            axes[idx].grid(axis="y")
        plt.savefig(dir_path+"gt_stat_outgroup_count.pdf")

        # ax = sns.boxplot(x="scale", y="RF-out", data=res_df, hue="locus_length", showmeans=True, showfliers = False)
        # plt.savefig(dir_path+"gt_stat_box.pdf")
    

if __name__ == "__main__":
    dir_path = sys.argv[1]
    re_end = int(sys.argv[2])
    unroot = (sys.argv[3].lower() == "true")
    ingroup = (sys.argv[4].lower() == "true")
    longer = (sys.argv[5].lower() == "true")
    # plot_box_polymorphic(dir_path)
    # plot_box_pdistance(dir_path)
    # plot_box_gt_error(dir_path)
    # print(unroot, ingroup)
    plot_gt_all(dir_path, re_end, unroot, ingroup, longer)
