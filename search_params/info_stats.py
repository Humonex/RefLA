
import os
import sys

import numpy as np
from collections import defaultdict
def stats_dic(dic, out_file):
    fo = open(out_file, "w")
    step = 0.001
    quantile_ls = [i * step for i in range(int(1 / step))]
    final_out = [[] for i in range(len(quantile_ls))]
    key_ls= [] 
    for key, ls in dic.items():
        key_ls.append(key)
        # new_ls = [a for a in ls if a != None]
        new_ls = []
        for a in ls:
            if a == None: 
                continue
            else:
                new_ls.append(a)
        # print("# Reg num: {}, Error_portion{:.4f}".format(len(new_ls), 1 - len(new_ls)))
        # print("# mean: {}, median: {}".format(np.mean(new_ls), np.median(new_ls)))
        fo.write("# --------{}--------\n".format(key))
        fo.write("# Raw reg num: {}, Normal reg num: {}, Error reg portion: {:.4f}\n".format(len(ls), len(new_ls), 1 - len(new_ls)/len(ls)))
        fo.write("# mean: {}, median: {}, min: {}, max: {}\n".format(np.mean(new_ls), np.median(new_ls), min(new_ls), max(new_ls)))
        ### 
        
        quantile_num_ls = np.quantile(new_ls, quantile_ls)
        # if key == "correct_portion":
        #     # print(new_ls)
        #     print(quantile_num_ls)
        # print(key, ":", new_ls)
        # if reverse_flag_dic[key]:
        #     quantile_num_ls = np.flip(quantile_num_ls)
        # ressults = []
        for i in range(len(quantile_ls)):
            # print("{:.3f}:{:.4f}".format(quantile_ls[i], quantile_num_ls[i]))
            # ressults.append(quantile_num_ls[i])
            final_out[i].append(quantile_num_ls[i])
    # print(final_out)
    ## 
    fo.write("#threashold\t" + "\t".join(key_ls) + "\n")
    for idx, ls in enumerate(final_out):
        fo.write("{:.3f}\t".format(quantile_ls[idx]))
        for a in ls:
            fo.write("{:.4f}\t".format(a))
        fo.write("\n")
    fo.close()

def add_ls(ls, a):
    if a == "nan":
        ls.append(None)
    else:
        # ls.append(float(a))
        try:
            ls.append(float(a))
        except:
            ls.append(None)
            pass
            # print(a)
def stats_info(file_in, file_out, key_ls):
    print("--------------------Start stats for {}--------------------".format(os.path.abspath(file_in)))
    with open(file_in, "r") as f:
        dic = defaultdict(list)
        remove_ls = ["contig", "start_pos", "end_pos"]
        # for key in key_ls:
        #     dic[key] = []
        print(key_ls)
        key1 = key_ls[0]
        for line in f:
            if line.startswith("#"): 
                # key_ls = line.strip("\n").strip("#").split("\t")
                continue
            if line.startswith(key1):
                continue
            if not line:
                continue
            # print(line)
            fields = line.strip("\n").split("\t")
            for i in range(len(fields)):
                if key_ls[i] in remove_ls:
                    continue
                add_ls(dic[key_ls[i]], fields[i])
        stats_dic(dic, file_out)

file_in = sys.argv[1]
key_ls = sys.argv[2].split(",")
file_out = os.path.abspath(file_in) + ".stats"

stats_info(file_in, file_out, key_ls)
# print("1".isalpha())
# dp,clip_num,correct_portion,differ_portion,disagree_portion
# contig,start,end,dp,clip_num,correct_portion,differ_portion,disagree_portion
# contig,start_pos,end_pos,dp,clip_num,correct_portion,differ_portion,disagree_portion