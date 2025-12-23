'''
统计待局部组装的区间大小、read nums

'''
import numpy as np
# import psutil
import os

def reg_to_id(reg):
    return reg[0] + ":" + str(reg[1]) + "-" + str(reg[2])
def id_to_reg(reg_str:str):
    ctg = reg_str.split(":")[0]
    start, end = reg_str.split(":")[1].split("-")
    start, end = int(start), int(end)
    return [ctg, start, end]
def stats(work_dir):
    dic = {}
    ls = []
    reg_size_ls =[]
    # if reg_ls is None:
    reg_ls = os.listdir(work_dir)
    # print(reg_ls)
    # return
    for reg_id in reg_ls:
        # reg_id = reg_to_id(reg)
        reg = id_to_reg(reg_id)
        ids_file = work_dir + "/" + reg_id + "/reg_denovo_ids.bed"
        read_num = 0
        with open(ids_file, "r") as f:
            for lines in f:
                read_num += 1
            ls.append(read_num)
            reg_size_ls.append(reg[2] - reg[1])
            dic[reg_id] = read_num
    print(np.sort(reg_size_ls))
    print(np.sort(ls))
    print(dic)
    pass
work_dir ="/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/thaliana_ont/my_pipe2/step3_SV_consensus/denovo_asm1/asm"
# bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/thaliana_ont/my_pipe2/step1_mapping/aln.sorted.bam"
stats(work_dir)