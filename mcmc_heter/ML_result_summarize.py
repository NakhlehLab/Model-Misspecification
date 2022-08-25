import matplotlib.pyplot as plt
import seaborn as sns
import sys

def get_best_ML_net(path):
    with open(path, "r") as handle:
        net = None
        score = None
        lines = handle.read().strip().split("\n")
        begin = False
        for line in lines:
            if "Inferred Network #1:" in line:
                begin = True
                continue
            if begin:
                if line.startswith("("):
                    net = line
                elif line.startswith("Total log probability:"):
                    score = line.split(": ")[1].strip()
                    break
        if score is not None:
            return net, float(score)
        else:
            return net, score

def get_all_ML_nets(rootpath):
    score_list = []
    for i in range(4):
        path = rootpath+str(i)+".out"
        net, score = get_best_ML_net(path)
        score_list.append(score)

    plt.plot(range(4), score_list, 'bo-')
    plt.xlabel("reticulations")
    plt.ylabel("log likelihood")
    plt.xticks([0,1,2,3])
    plt.ylim([ -60, -55])

    plt.savefig(rootpath+"_ll.pdf")
    plt.show()



if __name__ == '__main__':
    rootpaths = []
    # rootpaths.append("/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/ML/heter/truegt/ML_truegt")
    # rootpaths.append("/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/ML/heter/iqtreegt/ML_iqtree")
    # rootpaths.append("/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/ML/homo/truegt/ML_truegt")
    # rootpaths.append("/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/ML/homo/iqtreegt/ML_iqtree")
    # rootpaths.append("/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/ML/heter_locus/truegt/ML_truegt")
    # rootpaths.append("/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/ML/heter_locus/iqtreegt/ML_iqtree")
    # rootpaths.append(
    #     "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/net/ML/homo/truegt/ML_truegt")
    # rootpaths.append(
    #     "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/net/ML/homo/iqtreegt/ML_iqtree")
    # rootpaths.append(
    #     "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/net/ML/heter_locus/truegt/ML_truegt")
    # rootpaths.append(
    #     "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/net/ML/heter_locus/iqtreegt/ML_iqtree")
    # rootpaths.append(
    #     "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/net/ML/heter/truegt/ML_truegt")
    # rootpaths.append(
    #     "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/net/ML/heter/iqtreegt/ML_iqtree")
    #
    # for rootpath in rootpaths:
    #     get_all_ML_nets(rootpath)
    get_all_ML_nets(sys.argv[1])

