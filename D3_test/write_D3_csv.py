import pandas as pd

path="/home/zc36/X2/D3/main/gamma_0_sims/"
dict_res = {"D":[], "D3":[], "D_stdev":[], "D3_stdev":[]}
for i in range(1,101):
    sub_path = path + "gamma_"+str(i)+".out"
    with open(sub_path, "r") as handle:
        content = handle.read().split(" ")
        dict_res["D"].append(content[0])
        dict_res["D3"].append(content[1])
        dict_res["D_stdev"].append(content[2][1:len(content[2])-2])
        dict_res["D3_stdev"].append(content[3][1:len(content[3])-2])

df = pd.DataFrame(dict_res)
df.to_csv("D3_res_gamma_0.csv")