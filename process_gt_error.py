import pandas as pd
import sys
import seaborn as sns
import matplotlib.pyplot as plt

def summarize_iqtree_error_to_csv(dir, csv_path, rate_dir):
    data = {"index": [], "RF": [], "nrBS": [], "right_gt":[], "CoV":[]}
    
    for i in range(1, 11):
        with open(dir+"/"+"/"+str(i)+"/heter/"+rate_dir+"/gt_error.txt", "r") as handle:
            lines = handle.read().split("\n")
            data["index"].append(i)

            for line in lines:
                if "nrBS" in line:
                    data["nrBS"].append(line.strip().split(":")[1].strip())
                elif "RF" in line:
                    data["RF"].append(line.strip().split(":")[1].strip())
                elif "CoV of mutation" in line:
                    data["CoV"].append(line.strip().split(":")[1].strip())
                elif "correctly inferred" in line:
                    data["right_gt"].append(line.strip().split(":")[1].strip())
            
            handle.close()
    print(data)
    df = pd.DataFrame.from_dict(data)
    df.to_csv(csv_path)

def summarize_iqtree_error_to_csv(dir, csv_path, sub_dir, num_replicate, **kwargs):
    data = {"index": [], "RF-in": [], "RF-out":[], "nrBS": [], "right_gt":[], "CoV":[]}
    print(kwargs)
    for i in range(1, num_replicate+1):
        with open(dir+"/"+"/"+str(i)+"/heter/"+sub_dir+"/gt_error.txt", "r") as handle:
        # with open(dir+"/"+"/"+str(i)+"/heter/"+sub_dir+"/gt_error_ingroup.txt", "r") as handle:
            lines = handle.read().split("\n")
            data["index"].append(i)

            for line in lines:
                if "nrBS" in line:
                    data["nrBS"].append(line.strip().split(":")[1].strip())
                elif "ingroup RF" in line:
                    data["RF-in"].append(line.strip().split(":")[1].strip())
                elif "outgroup RF" in line:
                    data["RF-out"].append(line.strip().split(":")[1].strip()) 
                elif "CoV of mutation" in line:
                    data["CoV"].append(line.strip().split(":")[1].strip())
                elif "correctly inferred" in line:
                    data["right_gt"].append(line.strip().split(":")[1].strip())
            
            handle.close()
    print(data)

    for key, value in kwargs.items() :
        data[key] = [value for x in range(num_replicate)]
    df = pd.DataFrame.from_dict(data)
    df.to_csv(csv_path)

def plot_box_gt_error():
    dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/"
    locus_length = 1000
    df_list = []
    for scale in ["short", "medium", "long"]:
        csv_path = dir_path+scale+"/iqtree_error_ingroup"+str(locus_length)+".csv"
        df = pd.read_csv(csv_path)
        df["scale"] = [scale for x in range(len(df["index"]))]
        df_list.append(df)
        # mean = sum(df["diff"])/len(df["diff"])
        # std = statistics.pstdev(df["diff"])
        # print(scale, mean, std)
    res_df = pd.concat(df_list, ignore_index=True)
    print(res_df)
    res_df.to_csv(dir_path+"iqtree_error_ingroup"+str(locus_length)+".csv")
    ax = sns.boxplot(x="scale", y="RF", data=res_df)
    # ax.set_xticklabels(ax.get_xticklabels(),rotation=30)
    plt.savefig(dir_path+"iqtree_error_ingroup"+str(locus_length)+"_boxplot.pdf")
    plt.show()

if __name__ == "__main__":
# dir="/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_tall10/036/"
    # dir="/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa3_tall/036/"
    # dir = sys.argv[1]
    # if len(sys.argv) == 3:
    #     sub_dir = sys.argv[2]
    # else:
    #     sub_dir = ""
    for scale in ["short", "medium", "long"]:
        for locus_length in [500, 2000]:
            kwargs = {"scale":scale, "locus_length": locus_length}
            dir = sys.argv[1]+"/"+kwargs["scale"]+"/"
            csv_path = dir+"/iqtree_error"+str(kwargs["locus_length"])+".csv"
            # csv_path = dir+"/iqtree_error_ingroup"+str(kwargs["locus_length"])+".csv"
            num_replicate=int(sys.argv[2])
            summarize_iqtree_error_to_csv(dir, csv_path, str(kwargs["locus_length"]), num_replicate, **kwargs)
    # plot_box_gt_error()