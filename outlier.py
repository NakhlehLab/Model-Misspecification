from eval_5taxa import process_data, TAXON_MAP
import pandas as pd
from scipy.stats import chi2_contingency, power_divergence
import numpy as np
import dendropy
import scipy.special

NUM_LOCI = 10000
def read_obs(obs_dir):
    # i = 84
    i=1
    obs_path = obs_dir + str(i)+"/heter/genetrees.txt"
    exp_path = obs_dir + str(i)+"/heter/gtprob_true_0.csv"
    data = process_data(obs_path, exp_path, "Z 0", taxon_map=TAXON_MAP)
    df = pd.DataFrame(data)
    print(df)
    sumx = 0
    for i in range(105):
        x = (data["obs"][i]*NUM_LOCI - data["exp"][i]*NUM_LOCI)**2/(data["exp"][i]*NUM_LOCI)
        # print(np.array([[data["obs"][i]*100], [data["exp"][i]*100]]))
        # g, p, dof, expctd = chi2_contingency(np.array([[data["obs"][i]*100], [data["exp"][i]*100]]), lambda_="log-likelihood")

        print(i, data["obs"][i]*NUM_LOCI, data["exp"][i]*NUM_LOCI, x)
        sumx += x
    print(sumx)
    
    print(sum(data["obs"]*NUM_LOCI), sum(data["exp"]*NUM_LOCI))
    # g, p = power_divergence(data["obs"]*NUM_LOCI, data["exp"]*NUM_LOCI, lambda_='log-likelihood')
    # g2 = G(data["obs"]*NUM_LOCI, data["exp"]*NUM_LOCI)

    g, p = power_divergence([1,2], [2,1], lambda_='log-likelihood')
    g2 = G([1,2], [2,1])
    print(g, p, g2)
   


def G(obs, exp):
    s = 0
    for i in range(len(obs)):
        # print(obs[i], exp[i])
        # if exp[i] != 0:
        s += obs[i]*np.log(obs[i]/exp[i])
    return s*2

    # obs = np.asanyarray(obs)
    # exp = np.asanyarray(exp)
    # return 2.0 * scipy.special.xlogy(obs, obs / exp)

if __name__ == "__main__":
    # obs_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_tall10/036/"
    obs_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa5_gt/10000/"
    read_obs(obs_dir)
    # treestring = "(((wOo-HOST-Onchocerca_ochengi:0.001475838,wOv-HOST-Onchocerca_volvulus_str._Cameroon:0.0011135014)100:0.1163358521,wCfeJ-HOST-Ctenocephalides_felis:0.0769825256)39:0.010252609,wCle-HOST-Cimex_lectularius_JESC:0.0799948073,(wWb-HOST-Wuchereria_bancrofti:0.0087574512,(wBm-HOST-Brugia_malayi_strain_TRS:0.0000009234,wBpFR3-HOST-Brugia_pahangi:0.0012867426)73:0.0072423008)100:0.0815715978);"
    # t = dendropy.Tree.get(data=treestring, schema="newick", rooting="force-unrooted")
    # print(t)