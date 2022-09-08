import pandas as pd
import sys

def summarize_quarnetGof_pval(directory, re_start, re_end, truegt, ingroup):
    df_list = []
    df_summary = {"scale":[], "locus_length": [], "id": [], "p_value":[], "uncorrected_z": [], "sigma": []}
    for scale in ["short", "medium", "long"]:
        for i in range(re_start, re_end+1):
            for locus_length in [500, 2000]:
                if truegt:
                    if ingroup:
                        path_pval = directory +"/"+scale+"/"+ str(i) + "/heter/"+str(locus_length)+"/quarnetGoF_true_0_res_ingroup.txt"
                        path_csv = directory +"/"+scale+"/"+ str(i) + "/heter/"+str(locus_length)+"/quarnetGoF_true_0_stat_ingroup.csv"
                    else:
                        path_pval = directory +"/"+scale+"/"+ str(i) + "/heter/"+str(locus_length)+"/quarnetGoF_true_0_res.txt"
                        path_csv = directory +"/"+scale+"/"+ str(i) + "/heter/"+str(locus_length)+"/quarnetGoF_true_0_stat.csv"
                else:
                    if ingroup:
                        path_pval = directory +"/"+scale+"/"+ str(i) + "/heter/"+str(locus_length)+"/quarnetGoF_iq_0_res_ingroup.txt"
                        path_csv = directory +"/"+scale+"/"+ str(i) + "/heter/"+str(locus_length)+"/quarnetGoF_iq_0_stat_ingroup.csv"
                    else:
                        path_pval = directory +"/"+scale+"/"+ str(i) + "/heter/"+str(locus_length)+"/quarnetGoF_iq_0_res.txt"
                        path_csv = directory +"/"+scale+"/"+ str(i) + "/heter/"+str(locus_length)+"/quarnetGoF_iq_0_stat.csv"
                df_sub =pd.read_csv(path_csv)
                df_sub["scale"] = [scale for x in df_sub["p_value"]]
                df_sub["locus_length"] = [locus_length for x in df_sub["p_value"]]
                df_sub["id"] = [i for x in df_sub["p_value"]]
                df_list.append(df_sub)


                df_summary["scale"].append(scale)
                df_summary["locus_length"].append(locus_length)
                df_summary["id"].append(i)
                with open(path_pval, "r") as handle:
                    lines = handle.readlines()
                    df_summary["p_value"].append(lines[0].strip())
                    df_summary["uncorrected_z"].append(lines[1].strip())
                    df_summary["sigma"].append(lines[2].strip())

    df_res = pd.concat(df_list, ignore_index=True)
    df_summary = pd.DataFrame(df_summary)
    if truegt:
        if ingroup:
            df_res.to_csv(directory+"/quarnetGoF_true_outlierp_ingroup.csv")
            df_summary.to_csv(directory+"/quarnetGoF_true_summary_ingroup.csv")
        else:
            df_res.to_csv(directory+"/quarnetGoF_true_outlierp.csv")
            df_summary.to_csv(directory+"/quarnetGoF_true_summary.csv")
    else:
        if ingroup:
            df_res.to_csv(directory+"/quarnetGoF_iq_outlierp_ingroup.csv")
            df_summary.to_csv(directory+"/quarnetGoF_iq_summary_ingroup.csv")
            
        else:
            df_res.to_csv(directory+"/quarnetGoF_iq_outlierp.csv")
            df_summary.to_csv(directory+"/quarnetGoF_iq_summary.csv")

if __name__ == "__main__":
    # dir_path = "/home/zc36/X2/test_outgroup4/net0"
    re_start = 1
    re_end = 100
    truegt = True
    ingroup = True
    summarize_quarnetGof_pval(sys.argv[1], re_start, re_end, truegt, ingroup)

