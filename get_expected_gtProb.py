import pandas as pd


def read_gt_exp_csv(path):
    df = pd.read_csv(path, usecols=["gt"])
    print(df)
    df.to_csv("/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_gt/gt_topo.csv", index=False)

def compare_gt_prob(path1, path2):
    df1 = pd.read_csv(path1)
    df2 = pd.read_csv(path2, delimiter="\t", header=None)
    print(df1)
    print(df2)
    for topo in df1["gt"]:
        print(topo)

    dict2 = df2.to_dict()
    # print(dict1, dict2)

if __name__ == "__main__":
    # path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_gt/10000/gtprob_true_st.csv"
    path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_gt/gtprob_true_st2.csv"
    read_gt_exp_csv(path)
    
    # path1 = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_gt/10000/gtprob_true_st.csv"
    # path2 = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/tools/PRANC/BIN/outUnrGT.txt"
    # compare_gt_prob(path1, path2)