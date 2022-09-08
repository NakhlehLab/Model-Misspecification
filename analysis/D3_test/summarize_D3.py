import pandas as pd
import sys


def summarize(directory, re_start, re_end,  D3_type = 1):
    df_list = []
    for scale in ["short", "medium", "long"]:
        for i in range(re_start, re_end+1):
            for locus_length in [500, 2000]:
                if D3_type == 1:
                    path = directory +"/"+scale+"/"+ str(i) + "/heter/"+str(locus_length)+"/D3.txt"
                else:
                    path = directory +"/"+scale+"/"+ str(i) + "/heter/"+str(locus_length)+"/D3_2stage_faster.txt"
                try:
                    df_sub =pd.read_csv(path)
                    df_sub["scale"] = [scale for x in df_sub["D3"]]
                    df_sub["locus_length"] = [locus_length for x in df_sub["D3"]]
                    df_sub["id"] = [i for x in df_sub["D3"]]
                    df_list.append(df_sub)
                except:
                    print("Exception:"+path)

    df_res = pd.concat(df_list, ignore_index=True)
    if D3_type == 1:
        df_res.to_csv(directory+"/D3_res.csv")
    else:
        df_res.to_csv(directory+"/D3_res_2stage_faster.csv")

def summarize_subset():
    scale="long"
    locus_length = 500
    df_list = []
    remedy=3
    directory = "/home/zc36/X2/test_outgroup4/net0/"
    for i in [93, 95, 96]:
        if remedy == 0:
            path = directory +"/"+scale+"/"+ str(i) + "/heter/"+str(locus_length)+"/D3_500.txt"
        elif remedy == 1:
            path = directory +"/"+scale+"/"+ str(i) + "/heter/"+str(locus_length)+"/D3_1.txt"
        elif remedy == 2:
            path = directory +"/"+scale+"/"+ str(i) + "/heter/"+str(locus_length)+"/D3_bootknife.txt"
        elif remedy == 3:
            path = directory +"/"+scale+"/"+ str(i) + "/heter/"+str(locus_length)+"/D3_2stage.txt"
        df_sub =pd.read_csv(path)
        df_sub["scale"] = [scale for x in df_sub["D3"]]
        df_sub["locus_length"] = [locus_length for x in df_sub["D3"]]
        df_sub["id"] = [i for x in df_sub["D3"]]
        df_list.append(df_sub)

    df_res = pd.concat(df_list, ignore_index=True)
    if remedy == 0:
        df_res.to_csv(directory+"/D3_res_sub_500.csv")
    elif remedy == 1:
        df_res.to_csv(directory+"/D3_res_sub_ab1.csv")
    elif remedy == 2:
        df_res.to_csv(directory+"/D3_res_sub_ab_bootknife.csv")
    elif remedy == 3:
        df_res.to_csv(directory+"/D3_res_sub_ab_2stage.csv")


if __name__ == "__main__":
    # dir_path = "/home/zc36/X2/test_outgroup4/net0"
    dir_path = sys.argv[1]
    re_start = 1
    re_end = 100
    summarize(dir_path, re_start, re_end, int(sys.argv[2]))
    # summarize_subset()