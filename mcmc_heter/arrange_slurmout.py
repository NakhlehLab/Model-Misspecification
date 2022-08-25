import os
import sys


CHAIN_LEN = 20000000

def slurm_folder(sub_dir_path, file):
    slurm_path = sub_dir_path + file
    with open(slurm_path, 'r') as handle:
        last_line = handle.readlines()[-1].strip()
        # print(last_line)
        if (last_line.startswith("/scratch/zc36") or last_line.startswith("/home/zc36")) and ((str(CHAIN_LEN)+"_2/" in last_line) or (str(CHAIN_LEN)+"_1/" in last_line) or (str(CHAIN_LEN)+"_0" in last_line) or (str(CHAIN_LEN)+"_3" in last_line)):
            #do sth
            # print(slurm_path)
            last_line = last_line.replace("/scratch/", "/home/")
            if ("net0" not in last_line) and ("net1" not in last_line):
                netx = slurm_path[22]
                last_line = last_line.replace("/mcmc/", "/mcmc/net"+netx)
            if file.startswith("slurm"):
                print("cp "+ slurm_path + " " + last_line)
                os.system("cp "+ slurm_path + " " + last_line)
            else:
                print("cp "+ slurm_path + " " + last_line+"/slurm-123.out")
                os.system("cp "+ slurm_path + " " + last_line+"/slurm-123.out")


def arrange_files(dir_path):
    # dir_path = "/home/zc36/X2/mcmc/net0/long/"
    for i in range(1, 11):
        rep_dir_path = dir_path + str(i)
        for heter in ["heter", "homo"]:
            for murate in ["/murate/", "/"]:
                sub_dir_path = rep_dir_path + "/" + heter + "/2000" + murate
                for file in os.listdir(sub_dir_path):
                    if file.startswith("slurm") or file=="mcmc.out":
                        # path = sub_dir_path + file
                        # outpath = sub_dir_path
                        slurm_folder(sub_dir_path, file)
def arrange_all(netx):
    
    for scale in [ "long"]:
        dir_path =  "/home/zc36/X2/mcmc/"+netx+"/"+scale+"/"
        arrange_files(dir_path)

if __name__ == "__main__":
    # path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/long/1/heter/mcmcseq/0_0/10000000_0/slurm-4912525.out"
    # slurm_folder(path)
    # arrange_files()
    arrange_all(sys.argv[1])
    