# import numpy as np
# import os
from collections import namedtuple

def remove_overlaps(ls):
    ls.sort(key=lambda reg:reg[1])
    ls2 = []
    start = -1
    end = -1
    for reg in ls:
        if start >= 0:
            if reg[1] - end <= 0:
                end = max(end, reg[2])
            else:
                ls2.append([chr_id, start, end])
                start, end =  reg[1], reg[2]
        else:
            chr_id, start, end =  reg[0], reg[1], reg[2]
    if start >= 0:
        ls2.append([chr_id, start, end])
    return ls2

def cal_density(reg_ls, radius, t1, t2):    # t1表示超大区间，t2表示中等区间
    density_ls = []
    for idx, reg in enumerate(reg_ls):
        density = 0
        ## left
        for i in range(idx - 1, -1, -1):
            now_reg = reg_ls[i]
            if reg[1] - now_reg[2] > radius:break
            if now_reg[2] - now_reg[1] > t1:density += 3
            elif now_reg[2] - now_reg[1] > t2:density += 2
            else: density += 1
        ## right
        for i in range(idx, len(reg_ls)):   # 还包括自身
            now_reg = reg_ls[i]
            if now_reg[1] - reg[2] > radius:break
            if now_reg[2] - now_reg[1] > t1:density += 3
            elif now_reg[2] - now_reg[1] > t2:density += 2
        density_ls.append(density)
    return density_ls

def cluster_by_density(reg_ls, clu1, clu2, t1, t2, radius):
    ls = []
    # radius = 2000000
    # t1, t2 = 100000, 10000
    density_ls = cal_density(reg_ls, radius, t1, t2)
    pre_reg = None
    max_density = -1
    for idx, reg in enumerate(reg_ls):
        if pre_reg:
            max_density = max(density_ls[idx], density_ls[idx - 1])
            if max_density >= 3:clu = clu1
            else:clu = clu2
            ## 
            if reg[1] - pre_reg[2] < clu:
                pre_reg = [pre_reg[0], pre_reg[1], max(pre_reg[2], reg[2])]
            else:
                ls.append(pre_reg)
                pre_reg = reg
        else:pre_reg = reg
    if pre_reg:
        ls.append(pre_reg)
    return ls

def cluster_by_reg_size(reg_ls_in, reg_size, cluster_dis):
    ls = []
    ls2 = []
    for reg in reg_ls_in:
        if reg[2] - reg[2] > reg_size: ls.append([chr_id, reg[2], reg[2]])
    reg = None
    start = -1
    end = -1
    for reg in ls:
        if start >= 0:
            if reg[1] - end < cluster_dis:
                end = max(end, reg[2])
            else:
                ls2.append([chr_id, start, end])
                start, end =  reg[1], reg[2]
        else:
            chr_id, start, end =  reg[0], reg[1], reg[2]
    if start >= 0:
        ls2.append([chr_id, start, end])
    ls2.extend(reg_ls_in)
    return remove_overlaps(ls2)



if __name__ == "__main__":
    ## Debug
    pass
