import matplotlib.pyplot as plt
import seaborn as sns
import sys
import pandas as pd
import dendropy
from dendropy.calculate import treecompare
import os

def plot_polymorphism():
    dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/"
    outputpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/polymorphism_heter_500.pdf"
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # plt.clf()
    for i, scale in enumerate(["short", "medium", "long"]):
        inputpath = dir_path+scale+"/polymorphic_heter_500.csv"
        diffs = pd.read_csv(inputpath)
        print(diffs)
        sns.histplot(x=diffs["diff"], bins=30, stat="probability", ax=axes[i])
        axes[i].set_ylabel("Proportion", fontsize=16)
        axes[i].set_xlabel("Proportion of polymorphic sites", fontsize=16)
        axes[i].set_title(scale, fontsize=16, fontweight="bold")
        axes[i].set_ylim(0, 0.25)
        axes[i].set_xlim(0, 1.0)
        axes[i].tick_params(axis='y', labelsize=14)
        axes[i].tick_params(axis='x', labelsize=14)
    plt.tight_layout()
    plt.savefig(outputpath)
    plt.show()


def plot_non_info_sites_cnt(plot=True):
    dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/"
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    if not os.path.exists(dir_path + "non_info_sites_cnt.csv"):
        res_dict = {"scale": [], "replicate": [], "non_info_sites_cnt": []}
        for i, scale in enumerate(["short", "medium", "long"]):
            input_path = dir_path + scale + "/polymorphic_heter_500.csv"
            diffs = pd.read_csv(input_path)
            cnt = 0
            prev_rep = 0
            subcnt = 0
            for index, row in diffs.iterrows():
                if prev_rep != row["replicate"]:
                    res_dict["scale"].append(scale)
                    res_dict["replicate"].append(int(row["replicate"]))
                    res_dict["non_info_sites_cnt"].append(subcnt/10000)
                    # print(row["replicate"], subcnt/10000)
                    subcnt = 0
                    prev_rep += 1
                    
                if row["diff"] < 0.0001:
                    subcnt += 1
                    # print(row["replicate"])
                    cnt += 1
            print(f"number of loci with no informative sites: {cnt/1000000}")
            # if plot:
            #     sns.barplot(x="replicate", y="non_info_sites_cnt" , data=res_dict[res_dict["scale"] == scale], ax=axes[i])
            #     axes[i].set_ylabel("Probability", fontsize=16)
            #     axes[i].set_xlabel("Gene Tree", fontsize=16)
            #     axes[i].set_title(scale, fontsize=16, fontweight="bold")
            #     axes[i].set_ylim(0, 0.75)

        res_df = pd.DataFrame(res_dict)
        res_df.to_csv(dir_path + "non_info_sites_cnt.csv")
        
    else :
        res_df = pd.read_csv(dir_path + "non_info_sites_cnt.csv")
            
    if plot:
        for i, scale in enumerate(["short", "medium", "long"]):
            sns.barplot(x="replicate", y="non_info_sites_cnt" , data=res_df[res_df["scale"] == scale], ax=axes[i])
            axes[i].set_ylabel("Non informative sites proportion", fontsize=16)
            axes[i].set_xlabel("Replicate", fontsize=16)
            axes[i].set_title(scale, fontsize=16, fontweight="bold")
            axes[i].set_ylim(0, 0.05)
            axes[i].tick_params(axis='y', labelsize=12)
            axes[i].tick_params(axis='x', labelsize=12)
            axes[i].xaxis.set_ticks([x-1 for x in [20,40,60,80,100]])

        plt.tight_layout()
        plt.savefig(dir_path+"non_info_sites_cnt.pdf")
        plt.show()


def gt_non_info_distribution(polymorphic_path, gt_dir, num_replicate, locus_length, gt_type, non_info=True):
    tree_list = []
    tns = None
    print(gt_type)
    for i in range(1, num_replicate+1):
        if gt_type.lower() == "true":
            gt_path = gt_dir + str(i) + "/homo/"+str(locus_length)+"/genetrees.txt"
        elif gt_type.lower() == "false":
            gt_path = gt_dir + str(i) + "/heter/"+str(locus_length)+"/rooted_iqtree.txt"
        elif gt_type.lower() == "resolve":
            gt_path = gt_dir + str(i) + "/heter/"+str(locus_length)+"/rooted_iqtree_resolved.txt"
        if tns is None:
            replicate_tree_list = dendropy.TreeList.get(path=gt_path, schema="newick", rooting="force-rooted")
            tree_list.extend(replicate_tree_list)
            tns = tree_list[0].taxon_namespace
            if len(replicate_tree_list) != 10000:
                print(f"gt inference error:{i}")
        else:
            replicate_tree_list = dendropy.TreeList.get(path=gt_path, schema="newick", rooting="force-rooted", taxon_namespace=tns)
            tree_list.extend(replicate_tree_list)
            if len(replicate_tree_list) != 10000:
                print(f"gt inference error:{i}")
    
    if non_info:
        non_info_gt_cnt = {}
        diffs = pd.read_csv(polymorphic_path)
        for index, row in diffs.iterrows():
            if row["diff"] < 0.0001:
                find = False
                for x in non_info_gt_cnt.keys():
                    if treecompare.symmetric_difference(x, tree_list[index]) <= 0.00001:
                        non_info_gt_cnt[x] += 1
                        find = True
                        break
                if not find:
                    non_info_gt_cnt[tree_list[index]] = 1
        df = pd.DataFrame(list(non_info_gt_cnt.items()), columns = ['gt','count'])
        df.to_csv(gt_dir+"/non_info_gt_"+str(locus_length)+"_"+gt_type+".csv")
    else:
        gt_cnt = {}
        for index in range(len(tree_list)):
            find = False
            for x in gt_cnt.keys():
                if treecompare.symmetric_difference(x, tree_list[index]) <= 0.0001:
                    gt_cnt[x] += 1
                    find = True
                    break
            if not find:
                gt_cnt[tree_list[index]] = 1
        df = pd.DataFrame(list(gt_cnt.items()), columns = ['gt','count'])
        df.to_csv(gt_dir+"/gt_info_"+str(locus_length)+"_"+gt_type+".csv")


        
def run_all_gt_info(dir_path, num_replicate, locus_length, gt_type, non_info):
    for scale in ["short", "medium", "long"]:
        path = dir_path + scale + "/"
        # locus_length = 2000
        # truegt = True
        gt_non_info_distribution(path + "polymorphic_heter_"+str(locus_length)+".csv", path, num_replicate, locus_length, gt_type, non_info)
    
def plot_gt_non_info_distribution():
    dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/"
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    number_of_genes = 105
    locus_length = 500
    info_type = False
    truegt = "resolve"
    for i, scale in enumerate(["short", "medium", "long"]):
        if info_type == True:
            input_path = dir_path + scale + "/gt_info_"+str(locus_length)+"_"+str(truegt)+".csv"
            fig_path = dir_path+"gt_info_distribution_"+str(locus_length)+"_"+str(truegt)+".pdf"
        else:
            input_path = dir_path + scale + "/non_info_gt_"+str(locus_length)+"_"+str(truegt)+".csv"
            fig_path = dir_path+"non_info_gt_distribution_"+str(locus_length)+"_"+str(truegt)+".pdf"
        # elif truegt == "false" or truegt == "resolve":
        #     input_path = dir_path + scale + "/gt_info_"+str(locus_length)+"_"+str(truegt)+".csv"
        #     fig_path = dir_path+"gt_info_distribution_"+str(locus_length)+"_"+str(truegt)+".pdf"
        print(input_path)
        df = pd.read_csv(input_path)
        # print(df["count"])
        cnt_list = list(df["count"])
        cnt_list.sort(reverse=True)
        x = len(cnt_list)
        while x < number_of_genes:
            cnt_list.append(0)
            x += 1
        
        total = sum(cnt_list)
        cnt_list = [x/total for x in cnt_list]
        print(cnt_list)
        print(total)
        # sns.barplot(x=range(len(cnt_list)), y = cnt_list)
        axes[i].bar(range(number_of_genes), cnt_list)
        axes[i].set_ylabel("Probability", fontsize=16)
        axes[i].set_xlabel("Gene Tree", fontsize=16)
        axes[i].set_title(scale, fontsize=16, fontweight="bold")
        axes[i].set_ylim(0, 0.2)
        # plt.savefig(dir_path+scale+"/non_info_gt_distribution.pdf")
    plt.tight_layout()
    plt.savefig(fig_path)
    # plt.show()


if __name__ == "__main__":
    # plot_polymorphism()
    # run_all_gt_info(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5].lower() == "true")
    plot_gt_non_info_distribution()
    # plot_non_info_sites_cnt()
