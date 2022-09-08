import pandas as pd
import sys

def summarize(directory, re_start, re_end):
    df_list = []
    for scale in ["short", "medium", "long"]:
        for i in range(re_start, re_end+1):
            for locus_length in [500, 2000]:
                path = directory +"/"+scale+"/"+ str(i) + "/homo/"+str(locus_length)+"/D3_gt.txt"
                df_sub =pd.read_csv(path)
                df_sub["scale"] = [scale for x in df_sub["D3"]]
                df_sub["locus_length"] = [locus_length for x in df_sub["D3"]]
                df_sub["id"] = [i for x in df_sub["D3"]]
                df_list.append(df_sub)

    df_res = pd.concat(df_list, ignore_index=True)
    df_res.to_csv(directory+"/D3_gt_res.csv")

if __name__ == "__main__":
    # dir_path = "/home/zc36/X2/test_outgroup4/net0"
    dir_path = sys.argv[1]
    re_start = 1
    re_end = 100
    summarize( dir_path, re_start, re_end)