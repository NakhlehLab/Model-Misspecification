
def read_errors(path):
    with open(path, "r") as handle:
        for line in handle.readlines():
            if "Inferred Network #1:" in line:
                return True
        return False

def check_all(directory, num_reti, re_start, re_end):
    iq_rerun_list = []
    true_rerun_list = []
    for scale in ["short", "medium", "long"]:
        for i in range(re_start, re_end+1):
            for locus_length in [500, 1000]:
                path1 = directory +"/"+scale+"/"+ str(i) + "/heter/"+str(locus_length)+"/ML_true_"+num_reti+".out"
                if not read_errors(path1):
                    true_rerun_list.append(scale+"/"+str(i)+"/"+str(locus_length)+"_"+num_reti)
                
                path2 = directory+"/"+scale+"/"+ str(i) + "/heter/"+str(locus_length)+"/ML_iq_"+num_reti+".out"
                if not read_errors(path2):
                    iq_rerun_list.append(scale+"/"+str(i)+"/"+str(locus_length)+"_"+num_reti)
    print("rerun ML using IQTREE:")
    print(" ".join(iq_rerun_list))
    print("return ML using true gt:")
    print(" ".join(true_rerun_list))


if __name__ == "__main__":
    directory = "/home/zc36/X2/taxa5/net1/"
    num_reti = "0"
    re_start = 21
    re_end = 100
    check_all(directory, num_reti, re_start, re_end)