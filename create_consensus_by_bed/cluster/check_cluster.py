import numpy as np
import os
from collections import namedtuple
bed_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi/cluster/cluster0.bed"
all_ls = []
ctg_asm_ls = []
pre_ctg = None
with open(bed_in, "r") as f:
    for line in f:
        fields = line.strip().split("\t")
        ctg, start, end = fields[0], int(fields[1]), int(fields[2])
        if ctg == pre_ctg:
            ctg_asm_ls.append([start, end])
        else:
            all_ls.append([pre_ctg, ctg_asm_ls])
            ctg_asm_ls = [[start, end]]
            pre_ctg = ctg
    all_ls.append([pre_ctg, ctg_asm_ls])

fo = open("/public/home/hpc214712170/Test/tests/chm13_hifi/test_clu.bed", "w")
for ctg, ctg_asm_ls in all_ls:
    import cluster_for_reg
    if not ctg: continue
    print(ctg)
    ctg_asm_ls = [[ctg, reg[0], reg[1]] for reg in ctg_asm_ls]
    # print(ctg_asm_ls, "->")
    # print(ctg, cal_density(ctg_asm_ls, 2000000, 100000, 10000))
    print("clu_by_density: ", cluster_for_reg.cluster_by_density(ctg_asm_ls, 5000000, 300000))
    ls = cluster_for_reg.cluster_by_reg_size(ctg_asm_ls, 10000, 500000)
    # print(ls, "->")
    ls2 = cluster_for_reg.cluster_by_reg_size(ls, 50000, 1000000)
    print("clu_by_dis: ", ls2)
fo.close()


print(type(4//2))