'''
1、判断区间是否能够被读数跨过
能：使用读数修补
不能：进行局部组装
能够跨过条件: 1)区间为clip_reg; 2)区间长度<=MAX_SKIP_LEN; 3) 满足有n条高质量读数横跨区间, 区间两侧匹配较好(n>=min_support_reads, 质量>=min_MQ);
否则:不能跨过

2、局部组装修补   抽取剩余区间读数与unmapping reads合并, 进行local assembly, 然后对区间进行修补
1)可以修补的区间, 直接进行修补    可以修补的条件:略
2)无法修补: n填充


3、局部组装新思路：
抽取读数，构建string graph，进行组装

'''
'''
bed format2:    consensus.bed(记录对于参考的修改操作)
# chr   start   end     operation   info    
chr1    0       10      DEL         .         
chr1    12      13      INS         AGGGGGGG  
chr1    15      20      N_fill      NNNNNNNNNN
chr1    33      44      INV         .         
chr1    46      47      DUP         GGGAAATT  

//new final consensus bed | candidate bed 将
chr     start   end     operation       info                patch_id
chr1    100     120     reads_patch     AGCTGCT             read_id
chr1    100     120     asm_patch       AGCTGCT             asm_id
chr1    100     120     N_fill          NNNNNNNNNN(.)       .
'''
import argparse
import math
import os
import sys
import pysam
from multiprocessing import Pool
import random
import time
import re
import subprocess
from collections import namedtuple,defaultdict
from Bio import Seq
import shutil
import logging
import gzip
import copy
from create_consensus_by_bed import fasta_parser, solve_region_by_reads, merge_clip_regions, cluster_for_reg, apply_rec, solve_by_denovo, Scaffolding # linux pipe执行
from create_consensus_by_bed.Utils import Record, _enable_logging, Run_for_denovo, get_unmapped_reads, run_cmd_ls, make_dir, DepthRec, Connect_info, run_minimap2, run_samtools, reg_to_id, id_to_reg
from find_candidate_regions.find_reg_by_depth import Depth_info
from create_consensus_by_bed import extract_read_fastq_utils
# from create_consensus_by_bed.extract_read_fastq_utils import get_region_read_ids, select_reads_from_names
# import fasta_parser
# import solve_region_by_reads
# import merge_clip_regions
# import cluster_for_reg
# import apply_rec
# from Utils import Record, _enable_logging, Run_for_denovo
def call_back(res):
    print(res)

def error_call_back(error_code):
    print(error_code)

Region = namedtuple('Region', ["chr_id", "start", "end"])
Chr_info = namedtuple('Chr_info', ["chr_id", "chr_len"])

logger = logging.getLogger()

#################################################

###
def check_flank(bedinfo_ls, i):
    '''intervals'''
    target_reg = [bedinfo_ls[i][0], bedinfo_ls[i][1], bedinfo_ls[i][2]]  # target low_dep reg
    merge_ls = []
    dis_bias = 100
    merge_ls.append(target_reg)
    ## 向前找
    for j in range(i-1, -1, -1):
        if target_reg[1] - bedinfo_ls[j][2] <= dis_bias:
            merge_ls.append([bedinfo_ls[j][0], bedinfo_ls[j][1], bedinfo_ls[j][2]])
        else:
            break
    ## 向后找
    for j in range(i+1, len(bedinfo_ls)):
        if bedinfo_ls[j][1] - target_reg[2] <= dis_bias:
            merge_ls.append([bedinfo_ls[j][0], bedinfo_ls[j][1], bedinfo_ls[j][2]])
        else:
            break
    return merge_ls

def clip_cov_reg_merge(bedinfo_ls_in):  # ctg bedibnfo_ls from candidate bed
    '''
    specify ctg, merge low_dep reg and flank clip_reg
    bed format: ctg, start, end, info(clip_reg、low_dep)
    遍历区间, 找到low_dep就往左右搜索, 看是否能够合并
    '''
    if len(bedinfo_ls_in) <= 1:
        return bedinfo_ls_in
    merge_ls = []   # 
    bedinfo_ls_out = []
    for i in range(len(bedinfo_ls_in)):
        if bedinfo_ls_in[i][3] == 'low_dep':
            res = check_flank(bedinfo_ls_in, i)
            if len(res) > 1:
                res.sort(key=lambda i:i[1])
                r_bound = max([i[2] for i in res])  # res[-1][2]
                merge_ls.extend(res)
                bedinfo_ls_out.append([res[0][0], res[0][1], r_bound, "low_dep_clip"])
                # print("merge regions from:", res, " To: {}:{}-{}".format(res[0][0], res[0][1], r_bound))
    ## 
    for bedinfo in bedinfo_ls_in:
        if [bedinfo[0],bedinfo[1],bedinfo[2]] in merge_ls:
            continue
        else:
            bedinfo_ls_out.append(bedinfo)
    # for 
    return bedinfo_ls_out

###     
def reg_expand(ls_in, chr_len, config):
    ls = []
    reg_expand = config["reg_expand_dis"]   # 500
    for reg in ls_in:
        chr_id, start, end, info = reg[:4]  # 
        new_start = start - reg_expand if start - reg_expand > 0 else 0
        new_end = end + reg_expand if end + reg_expand < chr_len else chr_len
        ls.append([chr_id, new_start, new_end])
    # print("Reg expand:\n{}\n->{}".format(ls_in, ls))
    return ls


def reg_cluster_by_dis(reg_ls_in, dis): # 根据距离进行聚类，注意前面区间的end可能比后面区间要大
    t0 = time.time()
    if len(reg_ls_in) <= 1:
        return reg_ls_in
    reg_start = -1
    reg_end = -1
    chr_id = reg_ls_in[0][0]
    reg_ls_out = []
    reg_ls_in = sorted(reg_ls_in, key=lambda x: x[1])
    for reg in reg_ls_in:
        if reg_start > -1:  # 非首次
            if reg[1] - reg_end <= dis:  # 合并距离
                reg_end = max(reg[2], reg_end)
                need_to_cluster.append(reg)
            else:   # new_reg
                # reg_ls_out.append([chr_id, reg_start, reg_end, "clip_reg"])
                reg_ls_out.append([chr_id, reg_start, reg_end])
                reg_start = reg[1]
                reg_end = max(reg[2], reg_end)
                need_to_cluster = [reg]
        else:
            reg_start = reg[1]
            reg_end = reg[2]
            need_to_cluster = [reg]
    if reg_start > -1:
        # reg_ls_out.append([chr_id, reg_start, reg_end, "clip_reg"])   # 
        reg_ls_out.append([chr_id, reg_start, reg_end])
    print("reg_cluster_by_dis done, cost:{}s".format(time.time() - t0))
    return reg_ls_out

def rec_cluster_by_dis(rec_ls_in:list, dis):
    ls_out = []
    rec_ls_in.sort(key=lambda rec: rec.start)   # 凡聚类，先排序
    start, end = -1, -1
    for rec in rec_ls_in:
        if start > -1:
            if rec.start - end < dis:
                end = max(end, rec.end)   # 注意是选最大end; 原因在于clip 和lowdep聚类分别找到的区间并没有真正的一前一后，后面区间的end可能比前面区间的end 要小
            else:
                new_rec = Record(rec.chr_id, start, end)
                new_rec.add_operation("reg_asm")
                new_rec.add_info(".")
                new_rec.add_patch_id(".")
                ls_out.append(new_rec)
                # 
                start = rec.start
                end = rec.end
        else:
            start = rec.start
            end = rec.end
    if start > -1:
        new_rec = Record(rec.chr_id, start, end)
        new_rec.add_operation("reg_asm")
        new_rec.add_info(".")
        new_rec.add_patch_id(".")
        ls_out.append(new_rec)
    # print("cluster:\n{}->\n{}".format(rec_ls_in, ls_out))
    print("Cluster:from {}->{}".format(len(rec_ls_in), len(ls_out)))
    return ls_out

def write_rec_to_bed(bed_file, rec_ls):
    with open(bed_file, "w") as f:
        for rec in rec_ls:
            f.write("{}\t{}\t{}\t{}\n".format(rec.chr_id, rec.start, rec.end, rec.operation))

def write_reg_to_bed(bed_file, reg_ls):
    with open(bed_file, "w") as f:
        for reg in reg_ls:
            f.write("{}\t{}\t{}\n".format(reg[0], reg[1], reg[2]))

def cluster_by_asm_candidate(ls_in, chr_id, chr_len, bed_out_dir, config):
    ''' 根据asm 区域的大小等信息，使用asm区域合并周围小的asm区域，并将asm区域进行往两侧进行一定程度的延伸，优化asm region'''
    logger.info("{} cluster by asm candidate!!!".format(chr_id))
    ## 一些调试用的文件夹
    dir0 = bed_out_dir + "/merge0"
    dir1 = bed_out_dir + "/merge1"
    dir2 = bed_out_dir + "/merge2"
    dir3 = bed_out_dir + "/merge3"
    dir4 = bed_out_dir + "/merge4"
    # dir5 = bed_out_dir + "/no_merge"
    bed0 = dir0 + "/" + chr_id + ".bed"
    bed1 = dir1 + "/" + chr_id + ".bed"
    bed2 = dir2 + "/" + chr_id + ".bed"
    bed3 = dir3 + "/" + chr_id + ".bed"
    bed4 = dir4 + "/" + chr_id + ".bed"
    # bed5 = dir5 + "/" + chr_id + ".bed"
    for dir_x in [dir0, dir1, dir2, dir3, dir4]:
        if not os.path.isdir(dir_x):os.makedirs(dir_x)
    
    ## 聚类尺度参数 # for hifi data
    params = config["cluster_by_asm_candidate"]
    radius_1 = params["radius_1"]     # 30000      # 搜索周围20000bp的区域，合并reads_span区域
    radius_2 = params["radius_2"]     # 100000      # 搜索50000bp内的asm for cluster      建议比2*(radius1 + radius2)大   100000|150000
    radius_3 = params["radius_3"]         # 1000000
    # expand区间阈值
    reg_len_2 = params["reg_len_2"]    # 50000      # 
    reg_len_3 = params["reg_len_3"]        # 1000000     # 1MB以上的区间
    # 不同大小的区间作不同的扩充延展
    reg_expand_1 = params["reg_expand_1"]     # 20000   # 对于最后asm的区间处理，往左右两边进行一定的延展，最小的区间
    reg_expand_2 = params["reg_expand_2"]       # 40000    # 50000
    reg_expand_3 = params["reg_expand_3"]     # 100000    # 
    #
    ##
    # write_rec_to_bed(bed5, ls_in)
    ls = []
    target_ls = []    # 找到目标区域（asm 区域，合并侧翼的）    [[rec_idx, rec]....]
    ls_in.sort(key=lambda rec:rec.start)    # 聚类前先排序
    for idx, rec in enumerate(ls_in):
        op = rec.operation
        if op.endswith("asm"):
            target_ls.append([idx, rec])
            ls.append(rec)
    # 
    write_rec_to_bed(bed0, ls)   # 未聚类或者经过初始简单聚类
    # 第一轮简单聚类，引入周边reads_span
    merge_ls = []
    new_ls = []
    for idx, target_rec in target_ls:   # Merge 1
        new_start, new_end = target_rec.start, target_rec.end
        merge_ls = [target_rec]
        left = target_rec.start - radius_1 if target_rec.start - radius_1 > 0 else 0
        right = target_rec.end + radius_1 if target_rec.end + radius_1 < chr_len else chr_len
        new_radius_1 = radius_1 // 2
        for i in range(idx - 1, -1, -1):    # for j in range(i-1, -1, -1):
            if ls_in[i].end >= left:   # 合并
                merge_ls.append(ls_in[i])
                new_start = min(ls_in[i].start, new_start)  # 
                new_end = max(new_end, ls_in[i].end)    # 
                left = ls_in[i].start - new_radius_1 if ls_in[i].start - new_radius_1 > 0 else 0
                new_radius_1 = new_radius_1 // 2    # 防止存在过多小变异导致不断的合并
            else:
                break
        new_radius_1 = radius_1 // 2
        for i in range(idx + 1, len(ls_in)):   # for j in range(i+1, len(bedinfo_ls)):
            if ls_in[i].start <= right:
                merge_ls.append(ls_in[i])
                new_start = min(ls_in[i].start, new_start)  # 
                new_end = max(new_end, ls_in[i].end)
                right = ls_in[i].end + new_radius_1 if ls_in[i].end + new_radius_1 < chr_len else chr_len
                new_radius_1 = new_radius_1 // 2
            else:break

        ##  对merge_ls
        # merge_ls.sort(key=lambda rec:rec.start) 
        new_rec = Record(chr_id, new_start, new_end)    # 新的rec
        new_rec.add_info(".")
        new_rec.add_operation("asm")
        new_rec.add_patch_id(".")
        new_ls.append(new_rec)
        # if len(merge_ls) > 1:
        #     print("Merge from {} -> {}".format(merge_ls, ))   
    write_rec_to_bed(bed1, new_ls)  # for debug

    ## 两轮聚类: 大的聚类 big_cluster
    big_cluster_params = config["cluster_by_asm_candidate"]["big_cluster"]
    if big_cluster_params["apply_big_cluster"]:
        print("Apply big_cluster")
        reg_size = big_cluster_params["reg_size"]
        cluster_dis = big_cluster_params["cluster_dis"]
        radius = big_cluster_params["radius"]
        clu1 = big_cluster_params["clu1"]  # 一级聚类
        clu2 = big_cluster_params["clu2"]
        t1 = big_cluster_params["t1"]
        t2 = big_cluster_params["t2"]
        asm_reg_ls = [[rec.chr_id, rec.start, rec.end] for rec in new_ls]   # 转换成reg_ls
        # 聚类1
        ls = cluster_for_reg.cluster_by_reg_size(asm_reg_ls, reg_size, cluster_dis)   # 对所有asm都聚一遍
        write_reg_to_bed(bed2, ls)
        # 聚类2
        ls2 = cluster_for_reg.cluster_by_density(ls, clu1, clu2, t1, t2, radius)
        write_reg_to_bed(bed3, ls2)
        merge_rec_ls = []   # 转换成rec_ls
        for reg in ls2:
            rec = Record(reg[0], reg[1], reg[2])
            rec.add_info(".")
            rec.add_operation("asm")
            rec.add_patch_id(".")
            merge_rec_ls.append(rec)
    else:
        merge_rec_ls = new_ls


    ## Merge      # Merge local asm region and expand region
    start = -1
    end = -1
    rec = None
    new_rec = None
    ls = []
    merge_rec_ls.sort(key=lambda rec:rec.start)
    for rec in merge_rec_ls:
        if start >= 0:
            if rec.start - end < radius_2:
                end = rec.end
            else:
                if end - start > reg_len_3: reg_axpand = reg_expand_3
                elif end - start > reg_len_2: reg_axpand = reg_expand_2
                else: reg_axpand = reg_expand_1
                final_start = start - reg_axpand if start - reg_axpand > 0 else 0 
                final_end = end + reg_axpand if end + reg_axpand < chr_len else chr_len
                new_rec = Record(chr_id, final_start, final_end)
                new_rec.add_info(".")
                new_rec.add_operation("asm")
                new_rec.add_patch_id(".")
                ls.append(new_rec)
                ## 
                start = rec.start
                end = rec.end
        else:
            start = rec.start
            end = rec.end
    if start >= 0:
        if end - start > reg_len_3: reg_axpand = reg_expand_3
        elif end - start > reg_len_2: reg_axpand = reg_expand_2
        else: reg_axpand = reg_expand_1
        final_start = start - reg_axpand if start - reg_axpand > 0 else 0 
        final_end = end + reg_axpand if end + reg_axpand < chr_len else chr_len
        new_rec = Record(chr_id, final_start, final_end)
        new_rec.add_info(".")
        new_rec.add_operation("asm")
        new_rec.add_patch_id(".")
        ls.append(new_rec)
    write_rec_to_bed(bed4, ls)

    ## 进行最后的asm 区域之间的合并，将邻近的
    print("------------------Perform final denovo reg merge")
    asm_candidate_cluster_dis = config["get_candidate_op"]["asm_candidate_cluster_dis"]
    ls = rec_cluster_by_dis(ls, asm_candidate_cluster_dis)

    ls.extend(ls_in)    # 直接将初始的和新数组混在一起，后面在给过滤掉吧
    return ls

def remove_overlap_rec(rec_ls, chr_id):
    ls = []
    rec_ls.sort(key=lambda rec:rec.start)   # 
    start = -1
    end = -1
    merge_ls = []
    for rec in rec_ls:
        if start >= 0:
            if rec.start <= end:    # merge
                end = max(rec.end, end)
                need_to_merge.append(rec)
            else:
                if len(merge_ls) > 1:
                    new_rec =  Record(chr_id, start, end)
                    new_rec.add_info(".")
                    new_rec.add_operation("asm")
                    new_rec.add_patch_id(".")
                    ls.append(new_rec)
                else:
                    ls.append(merge_ls[0])
                need_to_merge = [rec]
                start, end = rec.start, rec.end
        else:
            start, end = rec.start, rec.end
            need_to_merge = [rec]
    if start >= 0:
        if len(merge_ls) > 1:
            new_rec =  Record(chr_id, start, end)
            new_rec.add_info(".")
            new_rec.add_operation("asm")
            new_rec.add_patch_id(".")
            ls.append(new_rec)
        else:
            ls.append(merge_ls[0])
    return ls



def get_candidate_op(bedinfo_ls_in:list, bam_in, bed_out_dir, Depthrec:DepthRec, config):  ## 返回Record ls
    '''处理通过特征得到的候选区间，进行裁剪合并，
    并筛选出哪些可以通过读数解决，哪些需要局部组装来解决
    single contig
    Format in: [chr1, 0, 100, "154478bp-cov_reg"] 
    '''
    if len(bedinfo_ls_in) < 1:
        return []
    ##
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    chr_id = bedinfo_ls_in[0][0]
    chr_len = bam_reader.get_reference_length(chr_id)
    bedinfo_ls_in.sort(key=lambda reg:reg[1])   # 进行排序

    ##
    # bedinfo_ls_in = clip_cov_reg_merge(bedinfo_ls_in)   # 
    candidate_op_ls_out = []
    clipinfo_ls = reg_expand(bedinfo_ls_in, chr_len, config)    # 全部都丢到一块进行聚类等处理吧，懒得搞了
    # print("clipinfo_ls:", clipinfo_ls)

    no_merge_dir = bed_out_dir + "/no_merge"
    no_merge_bed = no_merge_dir + "/" + chr_id + ".bed"
    with open(no_merge_bed, "w") as f:
        for rec in clipinfo_ls: f.write("{}\t{}\t{}\n".format(rec[0], rec[1], rec[2]))

    ##################################################################################
    ## 对clip_reg的作merge，1、just for del   ###容易出错的一步，将del的两侧合并起来
    ## 2、新增另一种判断模式，只要有一种判真，则合并区间
    ## 但是有很大可能是三个要合并，见鬼了。先不管了。
    ## 引入新的特征后这种聚合好像没啥卵用了，还毕竟费时
    logger.info("Step: clip merge !!!")
    clip_merge_dis = config["get_candidate_op"]["clip_merge_dis"]   # 500, 1000
    clipinfo_ls = reg_cluster_by_dis(clipinfo_ls, clip_merge_dis)  # 对合并的数据进行cluster
    logger.info("cluster by dis Done")
    # print("clustered_clipinfo_ls:", clipinfo_ls)
    if config["get_candidate_op"]["merge_by_clip_dp"]:
        merged_clipinfo_ls = merge_clip_regions.merge_by_clip_dp(clipinfo_ls, bam_in)
    else:
        merged_clipinfo_ls = clipinfo_ls
    # print("\nStep1_out:",merged_clipinfo_ls)

    ### 判断一下 哪些reg 是否可以使用读数进行填充，
    logger.info("Judge whether span by reads")
    solve_by_read_config = config["solve_by_read"]
    reads_span_candidate_op_ls = []
    asm_candidate_ls = []
    for rec in merged_clipinfo_ls:
        # res = judge_span_by_reads(rec[0], rec[1], rec[2], bam_in)     # 方法一
        can_solve, solve_seq, read_id = solve_region_by_reads.solve_by_reads(bam_in, rec[0], rec[1], rec[2], Depthrec, solve_by_read_config)
        if can_solve: # 需要判断是否返回的False
            record = Record(rec[0], rec[1], rec[2])
            record.add_info(solve_seq)
            record.add_operation("span_reads")  # merged_clip_reg / other_clip_reg
            record.add_patch_id(read_id)
            reads_span_candidate_op_ls.append(record)
        else:   # 不能横跨
            record = Record(rec[0], rec[1], rec[2])
            record.add_operation("reg_asm")
            record.add_info(".")
            record.add_patch_id(".")
            asm_candidate_ls.append(record)
    # print(merged_clipinfo_ls)
    # print("\nStep2:",asm_candidate_ls)


    ## 对asm_candidate 进行聚类 merge    //考虑将low_asm也一起聚类      这个聚类远远不够
    '''有一个问题：聚类后，原来位于二者之间能够span的区域仍会保留，所以在之后还要筛选一遍
    '''
    start = -1
    end = -1
    asm_candidate_cluster_dis = config["get_candidate_op"]["asm_candidate_cluster_dis"]    # # asm_candidates 聚类尺度 # 待调整
    asm_candidate_ls = rec_cluster_by_dis(asm_candidate_ls, asm_candidate_cluster_dis)
    candidate_op_ls_out.extend(asm_candidate_ls)    # 
    candidate_op_ls_out.extend(reads_span_candidate_op_ls)  ## 
    candidate_op_ls_out.sort(key=lambda rec:rec.start)

    ############################ 根据已有的 asm 区间对聚类 ##########################????
    ''' 根据asm 区域的大小等信息，使用asm区域合并周围小的asm区域，并将asm区域进行往两侧进行一定程度的延伸'''
    # write_rec_to_bed(bed_out_dir + "/debug.txt", candidate_op_ls_out)   # for debug
    candidate_op_ls_out = cluster_by_asm_candidate(candidate_op_ls_out, chr_id, chr_len, bed_out_dir, config)
    
    ###
    ############################# Remove 
    '''消除重叠的一些区间，主要是由于上面的合并会导致合并得到一些新的区间可能会将旧区间掩盖（特别是一些reads_span区间）'''
    print("candidate_op_ls_out:")
    candidate_op_ls_out.sort(key=lambda rec:rec.start)  # 聚类前先排序吧
    for rec in candidate_op_ls_out: print(rec.chr_id, rec.start, rec.end, rec.operation)
    final_candidate_op_ls = []
    start = -1
    end = -1
    need_to_merge = []
    for rec in candidate_op_ls_out:     # 消除交叠的区间 remove overlap
        if start > -1:
            if rec.start < end: # 与前面的区间有交叠，合并
                end = max(rec.end, end)     # 更新end
                need_to_merge.append(rec)
            else:   # 无交叠，将前面的区间加入 final_rec_ls 中  
                if len(need_to_merge) > 1:  # 要合并的话一定有一个是asm_region，直接
                    print("Merge {} overlap regions:".format(len(need_to_merge)))
                    has_asm = False
                    for REC in need_to_merge:   # 问题出在变量名的重复
                        print("{}:{}-{}".format(REC.chr_id, REC.start, REC.end, REC.operation), end=",")
                        if REC.operation.endswith("asm"):
                            has_asm = True
                            new_op = REC.operation
                    if not has_asm: # 
                        raise("ERROR region")
                    new_rec = Record(rec.chr_id, start, end)
                    new_rec.add_info(".")
                    new_rec.add_patch_id(".")
                    new_rec.add_operation(new_op)
                    final_candidate_op_ls.append(new_rec)
                    print("-> {}:{}-{}, {}".format(new_rec.chr_id, new_rec.start, new_rec.end, new_rec.operation))
                else:   # 不用合并的话保留原区间信息
                    final_candidate_op_ls.append(need_to_merge[0])  # 原区间存在了need_to_merge中
                    # print("不用合并：", need_to_merge[0].chr_id, need_to_merge[0].start, need_to_merge[0].end, need_to_merge[0].operation)
                # need_to_merge = [rec]   # 重置need_to_merge
                need_to_merge = []
                need_to_merge.append(rec)
                start = rec.start
                end = rec.end
        else:   # 首个区间
            start = rec.start
            end = rec.end
            need_to_merge = [rec]
    if start > -1:
        if len(need_to_merge) > 1:  # 要合并的话一定有一个是asm_region，直接
                    print("Merge {} overlap regions:".format(len(need_to_merge)))
                    has_asm = False
                    for rec in need_to_merge:
                        print("{}:{}-{}".format(rec.chr_id, rec.start, rec.end, rec.operation))
                        if rec.operation.endswith("asm"):
                            has_asm = True
                            new_op = rec.operation
                    print("\n")
                    if not has_asm:
                        raise("ERROR region")
                    new_rec = Record(rec.chr_id, start, end)
                    new_rec.add_info(".")
                    new_rec.add_patch_id(".")
                    new_rec.add_operation(new_op)
                    final_candidate_op_ls.append(new_rec)
                    print("-> {}:{}-{}, {}".format(new_rec.chr_id, new_rec.start, new_rec.end, new_rec.operation))
        else:   # 不用合并的话保留原区间信息
            final_candidate_op_ls.append(need_to_merge[0])  # 原区间存在了need_to_merge中

    '''对特殊情形的处理：如果染色体复杂区域过多，我们直接全部组装 '''
    if is_difficult_ctg(final_candidate_op_ls, chr_len, chr_id):
        new_rec = Record(chr_id, 0, chr_len)
        new_rec.add_info(".")
        new_rec.add_operation("whole_ctg_asm")  # 全染色体组装
        new_rec.add_patch_id(".")
        print("{} contig has too many asm reagions, will denovo for the whole contig!!!".format(chr_id))
        return [new_rec]    # 特殊情况


    print("final_candidate_op_ls:")
    for rec in final_candidate_op_ls: print(rec.chr_id, rec.start, rec.end, rec.operation)
    return final_candidate_op_ls


################################################# 用于get rec 和修改的函数 #################################################
def convert_reference_pos_to_raw_pos(read, candidate_pos, include_hardclip:bool):  # 输入read信息和要求的参考上的位置，将参考的位置转换为read上的位置
    '''将参考上某个坐标转换到原始读数上的坐标，所以需要将hardclip计算进去'''
    # candidate_pos = set(candidate_pos)
    # raw_ref_pos_map={}
    ref_pos = read.reference_start
    query_pos = 0
    for (ct,cl) in read.cigartuples:
        if ct==0:
            for i in range(cl):
                query_pos+=1
                ref_pos+=1
                if ref_pos == candidate_pos:
                    # raw_ref_pos_map[ref_pos] = query_pos
                    return query_pos
        elif ct==1:
            query_pos+=cl
        elif ct==2:
            for i in range(cl):
                ref_pos+=1
                if ref_pos == candidate_pos:
                    # raw_ref_pos_map[ref_pos] = query_pos
                    return query_pos
        elif ct==4:    
            query_pos+=cl
        elif ct==5:     # hard_clip
            if include_hardclip == True:    # 还需要考虑hardclip，
                query_pos += cl
        else:
            continue
    # return raw_ref_pos_map,read_length
    # return raw_ref_pos_map

def get_candidate_ctgs(bam_in, reg, config):
    '''
    min_align_length: 5000
    min_dis: 1000   # 5000
    '''
    ##
    left_candidate_ls = []
    right_candidate_ls = []
    candidate_dic = dict()
    candidate_ids_dic = dict()
    l_ids = set()
    r_ids = set()
    min_MQ = 20     # 太严格了，局部组装比对质量可能并不高
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    chr_id = reg.chr_id
    chr_len = bam_reader.get_reference_length(chr_id)   # 
    #
    min_dis = config["min_dis"]    # 只保留patch超过一定长度的
    min_align_length = config["min_align_length"]
    left_bound = max(reg.start - min_dis, 1)
    right_bound = min(reg.end + min_dis, chr_len - 1)
    for read in bam_reader.fetch(reg.chr_id, reg.start, reg.end):
        if read.is_secondary or read.mapping_quality < min_MQ or read.query_alignment_length < min_align_length:
            continue
        if read.reference_start < left_bound:    # patch左边界的点
            left_candidate_ls.append(read)
            l_ids.add(read.query_name)
        if read.reference_end > right_bound:        # patch右边界的点
            right_candidate_ls.append(read)
            r_ids.add(read.query_name)
    candidate_dic["left"] = left_candidate_ls
    candidate_dic["right"] = right_candidate_ls
    candidate_ids_dic["left"] = l_ids
    candidate_ids_dic["right"] = r_ids
    return candidate_dic

def get_info_op(candidate_ls, asm_fa_dic, direction, target_pos): # 
    '''
    target_pos为参考上跨过边界的点
    direction表示与参考良好match的那部分是左侧还是右侧
    '''
    # print(candidate_ls, direction, target_pos)
    ## get_best_ctg
    best_ctg = None
    for ctg in candidate_ls:    # 对候选ctg作过滤，加上对clip的限制，和长度的限制
        if best_ctg != None:
            if ctg.query_alignment_length > best_ctg.query_alignment_length:
                best_ctg = ctg
        else:
            best_ctg = ctg
    ## get info op
    if best_ctg != None:
        new_op = "asm_patch"
        ctg_seq = asm_fa_dic[best_ctg.query_name]
        # patch_id = best_ctg.query_name
        if best_ctg.is_reverse: # 判断是否反向链
            ctg_seq = Seq.reverse_complement(ctg_seq)
        if direction == "left":
            query_pos = convert_reference_pos_to_raw_pos(best_ctg, target_pos, True)
            new_info = ctg_seq[query_pos:]
            patch_id = best_ctg.query_name
        elif direction == "right":
            query_pos = convert_reference_pos_to_raw_pos(best_ctg, target_pos, True)
            new_info = ctg_seq[:query_pos]
            patch_id = best_ctg.query_name
        else:
            # raise("ERROR direction of {}:{}-{}".format(rec.chr_id, rec.start, rec.end))
            raise("ERROR direction")
    else:
        new_info = "."
        new_op = "N_fill"
        patch_id = "."
    return new_info, new_op, patch_id

def solve_t0(candidate_ls, rec, asm_fa_dic):    # 左端点的类型
    '''
    提供pass右边界的候选ctg, 初始的rec, 以及local_asm的fasta文件
    '''
    new_rec = Record(rec.chr_id, rec.start, rec.end)
    info, op, patch_id = get_info_op(candidate_ls, asm_fa_dic, "right", rec.end)  # 注意是pass右边界的
    new_rec.add_info(info)
    new_rec.add_operation(op)
    new_rec.add_patch_id(patch_id)
    if op == "asm_patch":
        print("{}:{}-{} left telomere can be asm_patch".format(rec.chr_id, rec.start, rec.end))
    else:
        print("{}:{}-{} left telomere N_fill".format(rec.chr_id, rec.start, rec.end))
    return new_rec


def solve_t1(candidate_ls, rec, asm_fa_dic):     # 右侧的类型
    '''
    提供pass左边界的候选ctg, 初始的rec, 以及local_asm的fasta文件
    '''
    new_rec = Record(rec.chr_id, rec.start, rec.end)
    info, op, patch_id = get_info_op(candidate_ls, asm_fa_dic, "left", rec.start)  # 注意是pass左边界的，找左侧匹配良好
    new_rec.add_info(info)
    new_rec.add_operation(op)
    new_rec.add_patch_id(patch_id)
    if op == "asm_patch":
        print("{}:{}-{} right telomere can be asm_patch".format(rec.chr_id, rec.start, rec.end))
        # logger
    else:
        print("{}:{}-{} right telomere N_fill".format(rec.chr_id, rec.start, rec.end))
    return new_rec


def solve_t2(left_candidate_ls, right_candidate_ls, rec, asm_fa_dic): # 非端点的类型，中间的区域
    new_rec = Record(rec.chr_id, rec.start, rec.end)

    left_ids = set()
    right_ids= set()
    for read in left_candidate_ls:
        # if read.
        left_ids.add(read.query_name)
    for read in right_candidate_ls:
        right_ids.add(read.query_name)
    common = left_ids.intersection(right_ids)   # 若某条读数同时跨过左右两侧，认为是最佳序列
    
    ## 
    if len(common) >= 1:
        print("{}:{}-{} can span by whole ctg".format(rec.chr_id, rec.start, rec.end))
        logger.info("{}:{}-{} can span by whole ctg".format(rec.chr_id, rec.start, rec.end))
        left_best_alignment_ls = []
        right_best_alignment_ls = []
        read_to_align_length = defaultdict(int)  # read_id: align_length
        if len(common) > 1:
            print("*********There are two span reads*********")
        for read in left_candidate_ls:
            if read.query_name in common:
                left_best_alignment_ls.append(read)
                read_to_align_length[read.query_name] += read.query_alignment_length
        for read in right_candidate_ls:
            if read.query_name in common:
                right_best_alignment_ls.append(read)
                read_to_align_length[read.query_name] += read.query_alignment_length
        ##
        best_ctg_id = max(read_to_align_length, key=read_to_align_length.get)  # 获取比对最长的ctg_id
        ## 
        '''有个问题是或许某个contig在某侧会有多个比对情况，极其复杂'''
        for read in left_best_alignment_ls:
            if read.query_name == best_ctg_id:
                left_query_pos = convert_reference_pos_to_raw_pos(read, rec.start, True)
                l_strand = read.is_reverse
        
        for read in right_best_alignment_ls:
            if read.query_name == best_ctg_id:
                right_query_pos = convert_reference_pos_to_raw_pos(read, rec.end, True)
                r_strand = read.is_reverse
        
        ## 得到最佳比对
        if l_strand != r_strand:
            print("{}:{}-{} has INV passed****".format(rec.chr_id, rec.start, rec.end))   ## need to improve！！！
            logger.info("{}:{}-{} has INV passed****".format(rec.chr_id, rec.start, rec.end))
            # new_rec.add_info(".")   
            new_rec.add_info("INV")     # ./INV
            new_rec.add_operation("N_fill")
            new_rec.add_patch_id(best_ctg_id)
        else:
            ctg_seq = asm_fa_dic[best_ctg_id]
            if l_strand: # 判断是否反向链
                ctg_seq = Seq.reverse_complement(ctg_seq)
            if left_query_pos < right_query_pos:
                info = ctg_seq[left_query_pos:right_query_pos]
                new_rec.add_info(info)
                new_rec.add_operation("asm_patch")
                new_rec.add_patch_id(best_ctg_id)
            elif left_query_pos == right_query_pos:     # 可能是del?? 如果直接有一个del跨过的话，会出问题，原因在于之前的区间选取的有问题
                print("{}:{}-{} may be del!!!".format(rec.chr_id, rec.start, rec.end))
                new_rec.add_info("")    # 
                new_rec.add_operation("asm_patch")      # 有待分析
                new_rec.add_patch_id(best_ctg_id)
            else:   # 有问题的pos
                print("{}:{}-{} ERROR query position!!!".format(rec.chr_id, rec.start, rec.end))
                # logger.info("{}:{}-{} ERROR query position!!!".format(rec.chr_id, rec.start, rec.end))
                # raise(ValueError)   # 调试阶段用    ???
                new_rec.add_info(".")
                new_rec.add_operation("N_fill")
                new_rec.add_patch_id(best_ctg_id)
        return new_rec
    else:   # common == 0
        print("{}:{}-{} can not span by whole ctg".format(rec.chr_id, rec.start, rec.end))
        logger.info("{}:{}-{} can not span by whole ctg".format(rec.chr_id, rec.start, rec.end))
        ## 
        left_info, left_op, left_patch_id = get_info_op(left_candidate_ls, asm_fa_dic, "left", rec.start)
        right_info, right_op, right_patch_id = get_info_op(right_candidate_ls, asm_fa_dic, "right", rec.end)
        if left_op == "N_fill" and right_op == "N_fill":    # 左右都是无穿过
            new_rec.add_info(".")
            new_rec.add_operation("N_fill")
            new_rec.add_patch_id(".")
        else:   # 可以部分解决，待完善，可以继续将部分读数收集  
            '''处理成 asm_patch,N_fill | N_fill,asm_patch | asm_patch,N_fill,asm_patch 这种形式？'''
            new_rec.add_info(left_info + "," + "." + "," + right_info)
            new_rec.add_operation(left_op + "," + "N_fill" + "," + right_op)
            new_rec.add_patch_id(left_patch_id + "," + "." + "," + right_patch_id)
        return new_rec

def solve_t3(chr_info:Chr_info, ):   # 处理低质量的contig
    new_rec = Record(chr_info.chr_id, 0, chr_info.chr_len)  # 
    # new_rec.add_info
    return new_rec


def is_difficult_ctg(rec_ls, chr_len, chr_id):     # 判断是否低质量contig
    threshold = 0.7
    asm_len = 0
    logger.info("asm_candidate_ls:->{}".format(chr_id))
    for rec in rec_ls:    # 
        if rec.operation.endswith("asm"):
            asm_len += rec.size
            logger.info("{} {} {}".format(rec.chr_id, rec.start, rec.end))
    asm_portion = asm_len / chr_len
    logger.info("{} contig asm_portion: {}, chr_len: {}".format(chr_id, asm_portion, chr_len))
    if asm_len / chr_len > threshold:
        logger.info("{} is difficult contig, {} portion is bad region!!!!".format(rec.chr_id, asm_len / chr_len))
        return True
    return False
def get_semi_rec(candidate_op_ls):
    '''
    1、process one ctg
    2、仅处理read patch区域
    '''
    final_op_ls = []
    for rec in candidate_op_ls:
        if rec.operation.endswith("reads"):
            rec.add_operation("reads_patch")
            final_op_ls.append(rec)
        elif rec.operation.endswith("asm"):
            final_op_ls.append(rec)
        else: raise ValueError

    return final_op_ls

def get_rec_by_asm(ctg, bam_in, candidate_op_ls, asm_fa_dic, config): # 传入asm_to_ref.bam, op_Ls
    '''
    process one ctg, 
    '''
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    final_op_ls = []
    asm_candidate_ls = []
    for rec in candidate_op_ls:
        # if rec.operation.endswith("reads"):
        #     rec.add_operation("reads_patch")
        #     final_op_ls.append(rec)
        #     logger.info("{}:{}-{} reads patch".format(rec.chr_id, rec.start, rec.end))
        # else:   # endswith("asm")
        #     asm_candidate_ls.append(rec)
        #     logger.info("{}:{}-{} asm".format(rec.chr_id, rec.start, rec.end))
        if rec.operation.endswith("reads"):
            rec.add_operation("reads_patch")
            final_op_ls.append(rec)
            # logger.info("{}:{}-{} reads_patch".format(rec.chr_id, rec.start, rec.end))
        elif rec.operation.endswith("skip"):
            final_op_ls.append(rec)
            print("Skip:{}:{}-{}".format(rec.chr_id, rec.start, rec.end))
        elif rec.operation.endswith("asm"):   # endswith("asm")  将所有asm区域收集
            # reg_id = rec.chr_id + ":" + str(rec.start) + "-" + str(rec.end)
            # if reg_id in asm_reg_ids:
            asm_candidate_ls.append(rec)
            # logger.info("{}:{}-{} asm".format(rec.chr_id, rec.start, rec.end))
        # 下面这些已经弃置
        # DeprecationWarning
        elif rec.operation.endswith("patch"):
            final_op_ls.append(rec)
            # logger.info("{}:{}-{}, {}".format(rec.chr_id, rec.start, rec.end, rec.operation))
        elif rec.operation == "replace_with_denovo":
            final_op_ls.append(rec)
        else:
            final_op_ls.append(rec)
            logger.info("ERROR rec {}:{}-{}, {}".format(rec.chr_id, rec.start, rec.end, rec.operation))
            raise("ERROR rec: {}:{}-{}, {}".format(rec.chr_id, rec.start, rec.end, rec.operation)) 
            # continue

    min_MQ = 20
    ##  for asm_candidate_ls
    for rec in asm_candidate_ls:
        ## 判断是否几乎whole contig 被认为是difficult region
        # ctg_len = ctg_len_dic[rec.chr_id]
        ctg_len = bam_reader.get_reference_length(rec.chr_id)
        if ctg_len == 0:
            print(rec.chr_id)
            exit(0)
        '''solve difficult contig'''
        if (rec.end - rec.start) / ctg_len > 0.8:    ## too difficult contig???一条染色体有一半是复杂区域，我们全部进行从头组装。使用从头组装结果来替代。
            info_ls = []
            patch_ls = []
            for read in bam_reader.fetch(rec.chr_id):
                if read.is_secondary or read.is_supplementary or read.mapping_quality < min_MQ:
                    continue
                info_ls.append(read.query_sequence)     # 收集所有的primary，将序列记录下来
                patch_ls.append(read.query_name)
            new_rec = Record(rec.chr_id, 0, ctg_len)
            new_info = ",".join(info_ls)    # 将局部组装的记录下来
            patch_ids = ",".join(patch_ls)
            for info in info_ls:
                new_info += str(info) + ","     # 将局部组装的记录下来
            ## 两类处理方式     由于部分染色体质量过差，选择全部从头组装。阈值设为0.7
            # new_rec.add_info(new_info)
            # new_rec.add_patch_id(patch_ids)     #
            new_rec.add_operation("replace_with_denovo")
            new_rec.add_info(".")
            new_rec.add_patch_id(".")   # 全部丢给从头组装
            print("{} replace_with_denovo: {}".format(rec.chr_id, patch_ls))
            logger.info("{} replace_with_denovo: {}".format(rec.chr_id, patch_ls))
            final_op_ls = [new_rec]
            return ctg, final_op_ls
        
        '''普通情形'''
        ## 找到pass左边和pass右边的最佳序列
        candidate_dic = get_candidate_ctgs(bam_in, Region(rec.chr_id, rec.start, rec.end), config)
        left_candidate_ls = candidate_dic["left"]
        right_candidate_ls = candidate_dic["right"]     # 
        ## 
        
        if rec.start < 1000: # 左端  取pass右端
            # print("{}:{}-{} left telomere".format(rec.chr_id, rec.start, rec.end))
            new_rec = solve_t0(right_candidate_ls, rec, asm_fa_dic)
        elif rec.end > ctg_len - 1000: # 右端，取pass左端
            # print("{}:{}-{} right telomere".format(rec.chr_id, rec.start, rec.end))
            new_rec = solve_t1(left_candidate_ls, rec, asm_fa_dic)
        else:
            # print("{}:{}-{} middle type".format(rec.chr_id, rec.start, rec.end))
            new_rec = solve_t2(left_candidate_ls, right_candidate_ls, rec, asm_fa_dic)
        final_op_ls.append(new_rec)

    return ctg, final_op_ls

# def SV_consensus_on_ref(candidate_op_ls, ref_seq, asm_fa_dic, bam_in, N_fill_size, final_rec_dic, final_seq_dic):
def SV_consensus_on_ref(chr_id, candidate_op_ls, ref_dic:dict, asm_fa_dic, bam_in, config):
    ''' 单线程的处理结果 '''
    ''' 注意提供的bam是asm_to_ref的bam, 提供新组装结果的consensus_fa_dic '''
    ''' 多线程操作字典修改可能会冲突，需要每个染色体单独返回一个字典 '''
    ctg, final_rec_ls = get_rec_by_asm(chr_id, bam_in, candidate_op_ls, asm_fa_dic, config["solve_by_denovo"])
    # new_seq_ls = apply_rec_on_ref(final_rec_ls, ref_seq, N_fill_size, consensus_fa_dic)
    connect_info, ctg_consensus_fa_dic = apply_rec.apply_rec_on_ref2(chr_id, final_rec_ls, ref_dic)
    # return final_rec_ls, new_seq_ls
    return final_rec_ls, connect_info, ctg_consensus_fa_dic
def SV_consensus_on_ref2(all_chrs, candidate_op_dic, ref_dic:dict, asm_fa_dic, bam_in, SVconsensus_bed_out, consensus_fasta_out, work_dir, threads, config):
    ''' 单线程的处理结果 '''
    ''' 注意提供的bam是asm_to_ref的bam, 提供新组装结果的consensus_fa_dic '''
    ''' 多线程操作字典修改可能会冲突，需要每个染色体单独返回一个字典 '''
    rec_dic = defaultdict(list)
    final_rec_ls = []
    pool = Pool(processes=threads)
    results = [pool.apply_async(Scaffolding.get_rec_by_asm, args=(ctg, bam_in, candidate_op_dic[ctg], asm_fa_dic, config["solve_by_denovo"])) for ctg in all_chrs]
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    for i, res in enumerate(results):
        ctg, ctg_final_rec_ls = res.get()
        rec_dic[ctg] = ctg_final_rec_ls
        final_rec_ls.extend(ctg_final_rec_ls)
    ## 
    connect_info_ls, patch_ids, consensus_dic = Scaffolding.consensus(bam_in, ref_dic, rec_dic, asm_fa_dic, consensus_fasta_out, work_dir, config, True)
    ## 
    # connect_info_ls = []
    # for ctg_id, ctg in consensus_dic.items():
    #     connect_info = Connect_info(ctg_id, [], []) # 
    #     connect_info.connect_ls.append(ctg_id)
    ## 
    Record.write_record(final_rec_ls, SVconsensus_bed_out)
    return connect_info_ls, patch_ids, consensus_dic

def run_SV_consensus_on_ref(threads, ctg_ls, SVconsensus_bed_out, consensus_fasta_out, semi_candidate_op_dic, asm_fa_dic, ref_dic, consensus_fa_dic, asm_to_ref_bam, Nfill_size):
    # run
    logger.info("Run SV_consensus_on_ref")
    pool = Pool(processes=threads)
    results = [pool.apply_async(SV_consensus_on_ref, args=(ctg, semi_candidate_op_dic[ctg], ref_dic, asm_fa_dic, asm_to_ref_bam, Nfill_size)) for ctg in ctg_ls]
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    connect_info_ls = []
    with open(SVconsensus_bed_out, "w") as fo2:     # 非并行写文件
        fo2.write("#chr_id\tstart\tend\topertion\tinfo\tpatch_id\n")
        for i, res in enumerate(results):
            final_rec_ls, connect_info, ctg_consensus_fa_dic = res.get()
            # ctg = all_chrs[i]
            connect_info_ls.append(connect_info)
            consensus_fa_dic.update(ctg_consensus_fa_dic)   # 将每个线程的字典 consensus 结果加入进来
            ## write record bed
            for rec in final_rec_ls:
                fo2.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(rec.chr_id, rec.start, rec.end, rec.operation, rec.info, rec.patch_id))
#################################################get rec的函数结束#################################################
def apply_reads_consensus():
    pass


####################



def candidate_process(ctg_bedinfo_ls, bam_in, bed_out, Depthrec, config):
    # print("Process {}".format(ctg_bedinfo_ls[0][0]))
    if ctg_bedinfo_ls:
        print("\nProcess {}".format(ctg_bedinfo_ls[0][0]))
        logger.info("\nProcess {}".format(ctg_bedinfo_ls[0][0]))
    ctg_candidate_op_out = get_candidate_op(ctg_bedinfo_ls, bam_in, os.path.dirname(bed_out), Depthrec, config)    # 比较费时，主要在于判断是否为del_pair
    with open(bed_out, "w") as fo:
        for rec in ctg_candidate_op_out:
            fo.write("\t".join([rec.chr_id, str(rec.start), str(rec.end), rec.operation, rec.info, rec.patch_id]))
            fo.write("\n")
    return ctg_candidate_op_out

def load_candidate_reg_bed(candidate_bed):  # 加载第一步通过特征找到的区间
    t0 = time.time()
    bedinfo_dic = {}
    with open(candidate_bed, 'r') as f:
        for line in f:
            fields = line.strip().split("\t")
            ctg = fields[0]
            if ctg not in bedinfo_dic:
                bedinfo_dic[ctg] = [[fields[0], int(fields[1]), int(fields[2]), fields[3]]]
            else:
                bedinfo_dic[ctg].append([fields[0], int(fields[1]), int(fields[2]), fields[3]])
    print("Load candidate bed, done cost:{}s".format(time.time() - t0))
    return bedinfo_dic

def get_frag_reads(fa_in, patch_ids):
    '''收集local denovo在后面的修补中用掉的碎片contig，返回剩余的序列片段'''
    '''chr     start   end     operation       info                patch_id'''
    ## get fragement sequence
    frag_dic = {}
    fa_dic = fasta_parser.read_sequence_dict(fa_in)
    for ctg,seq in fa_dic.items():
        if ctg not in patch_ids:  # exclude_ids -> patch_ids
            frag_dic[ctg] = seq
    return frag_dic

def get_depthrec(dp_file, chr_id, chr_len, win_size, block_size, whole_dp):
    dp_ls = DepthRec.read_mosdepth_dp_file(dp_file)
    ctg_depth_rec = DepthRec(chr_id, chr_len, dp_ls, win_size, block_size, whole_dp)
    return ctg_depth_rec

'''
改进：优化读数缺失问题,收集高剪切primary读数用于组装。从所有读数中搜/从所有区间搜

'''
def run_SVconsensus_parallel(out_dir, bam_in, process_ctg_ls,
                       candidate_bed, threads, fastq_in, data_type,
                       reference_fn, Nfill_size, genome_size,
                       ex_unmapped_denovo, out_denovo, keep_ls, args, config):
    ## 

    my_log = out_dir + "/get_fasta_consensus2.log"
    _enable_logging(my_log, debug=False, overwrite=True)
    ref_dic = fasta_parser.read_sequence_dict(reference_fn)
    t1 = time.time()
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai", threads=threads)
    all_chrs = bam_reader.references
    logger.info("Process contigs: {} !!!".format(process_ctg_ls))

    ## get depths_file
    dpinfo_f = os.path.dirname(out_dir) + "/step2_candidate_regions/dp_info/dpinfo.bed"     # 获取dpinfo_bed
    dpinfo_dic = Depth_info.read_dp_info(dpinfo_f)
    whole_dp = dpinfo_dic[process_ctg_ls[0]]["whole_dp"]    # get whole_dp
    task_ls = []
    win_size = config["depths"]["win_size"]
    block_size = config["depths"]["block_size"]  # config
    for chr_id in process_ctg_ls:
        dp_file = os.path.dirname(out_dir) + "/step2_candidate_regions/depths/" + chr_id + "/" + chr_id + ".regions.bed.gz"
        chr_len = bam_reader.get_reference_length(chr_id)
        task_ls.append([dp_file, chr_id, chr_len, win_size, block_size, whole_dp])
    pool = Pool(processes=threads)
    results = [pool.apply_async(get_depthrec, args=(task)) for task in task_ls]
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    depthrec_dic = {}
    for res in results:
        parse_res = res.get()
        depthrec_dic[parse_res.chr_id] = parse_res
    print("get depth done")

    ### ****************** [Candidate bed regions process] ******************
    ## load candidate bed file  
    bedinfo_dic = load_candidate_reg_bed(candidate_bed)
    ## parallel process ctgs   
    candidate_op_dir = out_dir + "/candidate_op"
    dirs = [candidate_op_dir, candidate_op_dir+"/no_merge", candidate_op_dir+"/merge0", candidate_op_dir+"/merge1", candidate_op_dir+"/merge2", candidate_op_dir+"/merge3", candidate_op_dir+"/merge4"]
    for dir in dirs:
        if not os.path.exists(dir):os.makedirs(dir)
    all_reference_ids = process_ctg_ls
    bed_to_merge = [candidate_op_dir + "/candidate_op_" + ctg + ".bed" for ctg in all_reference_ids]
    pool = Pool(processes=threads)
    results = [pool.apply_async(candidate_process, args=(bedinfo_dic.get(all_reference_ids[id], []), bam_in, bed_to_merge[id], depthrec_dic[all_reference_ids[id]], config)) for id in range(len(all_reference_ids))]
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    '''注意：apply_async的话需要使用get来得到返回值'''
    candidate_op_dic = defaultdict(list)   # 将所有染色体都加进去，便于统一
    # candidate_op_ls = []
    for res in results:     # res is ctg_candidate_op_ls
        parse_res = res.get()
        if len(parse_res) > 0:
            candidate_op_dic[parse_res[0].chr_id].extend(parse_res) # chr_id->rec_ls
            # candidate_op_ls.extend(res)
    # return
    ## merge bed  out: candidate_op.bed
    merge_out = candidate_op_dir + "/candidate_op.bed"
    merge_cmd_ls = ["cat", " ".join(bed_to_merge), "|", "sort -k 1,1 -k 2n,2", ">", merge_out]
    run_cmd_ls(merge_cmd_ls)
    t2 = time.time()
    # return

    if args.skip_denovo:    # correct mode, 跳过后面的denovo
        '''还需要完事correct mode'''
        SV_consensus_dir = out_dir
        make_dir(SV_consensus_dir)
        SVconsensus_bed_out = SV_consensus_dir + "/consensus.bed"
        consensus_fasta_out = SV_consensus_dir + "/consensus.fasta"
        consensus_fa_dic = {}
        # correct mode
        print("----------Skip denovo, correct mode----------")
        if args.cut:
            asm_opertion = "asm_cut"
            print("cut asm region")
        else:
            asm_opertion = "asm_skip"
            print("skip asm region")
        # 
        final_op_dic = defaultdict(list)
        for ctg in candidate_op_dic.keys(): # 获取final rec ls
            op_ls = candidate_op_dic[ctg]
            new_op_ls = []
            for rec in op_ls:
                if rec.operation.endswith("reads"): # ！！！理论上这种情况应该会很少，由读数填充的区域  ->改进
                    rec.add_operation("reads_patch")
                    new_op_ls.append(rec)
                elif rec.operation.endswith("asm"):   # endswith("asm")  将所有asm区域收集
                    rec.add_operation(asm_opertion)
                    new_op_ls.append(rec)
                else:
                    raise ValueError
            final_op_dic[ctg] = new_op_ls
        ## 
        connect_info_ls = []
        pool = Pool(processes=threads)
        results = [pool.apply_async(apply_rec.apply_rec_on_ref2, args=(ctg, final_op_dic[ctg], ref_dic)) for ctg in all_chrs]
        pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
        pool.join() # 等待进程池中的所有进程执行完毕
        for i, res in enumerate(results):
            connect_info, ctg_consensus_fa_dic = res.get()
            # ctg = all_chrs[i]
            connect_info_ls.append(connect_info)
            consensus_fa_dic.update(ctg_consensus_fa_dic)   # 将每个线程的字典 consensus 结果加入进来
        all_final_ls = []
        for ctg in final_op_dic.keys():
            all_final_ls.extend(final_op_dic[ctg])
        Record.write_record(all_final_ls, SVconsensus_bed_out)
        fasta_parser.write_fasta_dict(consensus_fa_dic, consensus_fasta_out)
        return
    else:   # 非correct mode，不跳过
        pass
    
    
    #### ****************** [local denovo assembly] ******************
    logger.info("Run local assembly module !!! ")
    denovo_asm_dir = out_dir +"/denovo_asm"     # 母目录, 存的是所有denovo区域的数据
    make_dir(denovo_asm_dir)
    ## 1、collect reads of local assembly regions and unmapped reads, 
    ## then extract fastq from reads_set()
    unmapped_fq = denovo_asm_dir + "/unmapped.fastq"
    get_unmapped_reads(threads, bam_in, unmapped_fq, "fastq")
    ### ----------------------------------new pipe---------------------------------- ###
    candidate_op_ls = []
    for op_ls in candidate_op_dic.values():
        candidate_op_ls.extend(op_ls)
    #  
    reg_read_ids_dic, clip_read_ids_dic = extract_read_fastq_utils.get_region_read_ids(out_dir, bam_in, threads, candidate_op_ls, whole_dp, config)   # asm_reg_id -> reg_read_ids
    asm_region_fq = denovo_asm_dir + "/region.fastq"
    region_ids_fn = denovo_asm_dir + "/region_ids.bed"
    extract_read_fastq_utils.write_read_ids(region_ids_fn, reg_read_ids_dic["all_read"])
    solve_by_denovo.select_reads_from_names(fastq_in, asm_region_fq, region_ids_fn, threads)
    # 
    print("Get raw fastq denoe!!!")
    t3 = time.time()

    ## ----------------------------------2、apply reads patch----------------------------------
    if config["apply_read_patch_first"]:
        print("------------------------------Apply reads patch------------------------------")
        mm_num = max(math.ceil(config["high_clip"]["min_high_clip_num"]), config["high_clip"]["min_portion"] * whole_dp)
        ##  抽取clip区域的读数
        clip_dir = out_dir + "/clip"
        make_dir(clip_dir)
        raw_clip_read_ids_fn = clip_dir + "/raw_clip_ids.bed"
        raw_clip_read_fq = clip_dir + "/raw_clip_ids.fastq"
        clip_bam = clip_dir + "/clip2ref.bam"
        clip_sam = clip_dir + "/clip2ref.sam"
        extract_read_fastq_utils.write_read_ids(raw_clip_read_ids_fn, reg_read_ids_dic["clip_read"])
        solve_by_denovo.select_reads_from_names(asm_region_fq, raw_clip_read_fq, raw_clip_read_ids_fn, threads)
        ##
        semi_consensus_bed = out_dir + "/semi_consensus.bed"
        semi_consensus_fa = out_dir + "/semi_consensus.fa"
        semi_candidate_op_dic = defaultdict(list)
        reg_trans = {}
        semi_consensus_fa_dic = {}
        results = []
        pool = Pool(processes=threads)
        for chr_id in all_chrs:
            rec_ls = apply_rec.get_semi_rec(candidate_op_dic[chr_id])
            results.append(pool.apply_async(apply_rec.apply_read_patch, args=(chr_id, rec_ls, ref_dic), error_callback=error_call_back))
        pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
        pool.join() # 等待进程池中的所有进程执行完毕
        for res in results:
            semi_ctg_op_dic, ctg_consensus_fa_dic, ctg_reg_trans = res.get()
            semi_candidate_op_dic.update(semi_ctg_op_dic)
            # print("ctg_reg_trans:", ctg_reg_trans)
            reg_trans.update(ctg_reg_trans)
            semi_consensus_fa_dic.update(ctg_consensus_fa_dic)
        fasta_parser.write_fasta_dict(semi_consensus_fa_dic, semi_consensus_fa)
        # 
        semi_op_ls = []
        for key in semi_candidate_op_dic.keys():
            semi_op_ls.extend(semi_candidate_op_dic[key])
        Record.write_record(semi_op_ls, semi_consensus_bed)
        # print(reg_trans)
        ## map clip reads to semi fasta 
        print("-------------------------map clip reads to semi fasta---------------------------")
        run_minimap2(semi_consensus_fa, raw_clip_read_fq, data_type, threads, clip_sam, ["-w16 -k24"])
        run_samtools(clip_sam, clip_bam, threads)
        ##
        clip_reg_fille = clip_dir + "/clip_reg.bed"
        final_clip_ids = set()
        results = []
        pool = Pool(processes=threads)
        for reg in clip_read_ids_dic.keys():    # 是否选择使用全部的读数填充缺口
            new_reg_id = reg_trans[reg]
            try:
                new_reg = id_to_reg(new_reg_id)
            except:
                print("Error id:", new_reg_id)
                exit()
            # flag, new_reg, ids_set = extract_read_fastq_utils.get_clip_read_ids(bam_in, threads, new_reg, mm_num, config)
            results.append(pool.apply_async(extract_read_fastq_utils.get_clip_read_ids, args=(clip_bam, threads, new_reg, mm_num, config), error_callback=error_call_back))
        pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
        pool.join() # 等待进程池中的所有进程执行完毕
        with open(clip_reg_fille, "w") as f:
            for res in results:
                flag, reg, ids_set = res.get()
                final_clip_ids.update(ids_set)
                if len(ids_set) >= mm_num:
                    f.write("{}\t{}\t{}\t{}\n".format(reg[0], reg[1], reg[2], len(ids_set)))
        print("Raw clip nuqm: {},-> Final clip read num:{}".format(len(reg_read_ids_dic["clip_read"]), len(final_clip_ids)))
        clip_read_ids = final_clip_ids
        fa_for_consensus = semi_consensus_fa
    else:
        fa_for_consensus = reference_fn
        clip_read_ids = reg_read_ids_dic["clip_read"]
        semi_candidate_op_dic = candidate_op_dic
    ## ----------------------------------3、Start assembly----------------------------------
    print("----------------------------------Start assembly----------------------------------")
    asm_regions, whole_ctg_asm_reg_ls = [], []
    for rec in candidate_op_ls:
        if rec.operation.endswith("asm"):
            asm_regions.append([rec.chr_id, rec.start, rec.end])
            if rec.operation == "whole_ctg_asm":
                whole_ctg_asm_reg_ls.append([rec.chr_id, rec.start, rec.end])
    
    # 2.1、ctg denovo
    ctg_denovo_dir = out_dir +"/ctg_denovo"
    make_dir(ctg_denovo_dir)
    print("Start ctg denovo")
    asm_ctg_ls = []
    if config["apply_ctg_denovo"]:
        whole_ctg_asm_reg_ls = solve_by_denovo.filter_empty_ctg(whole_ctg_asm_reg_ls, reg_read_ids_dic)
        reg_ls = []
        if config["denovo_by_ctg"]["apply_all"]:
            print("Apply on all denovo ctg")
            asm_ctg_ls = [reg[0] for reg in whole_ctg_asm_reg_ls]
            reg_ls = whole_ctg_asm_reg_ls
        else:
            for reg in whole_ctg_asm_reg_ls:
                if solve_by_denovo.is_abnormal_ctg(dpinfo_dic[reg[0]], config["denovo_by_ctg"]):
                    print("{} is abnormal ctg, dpinfo: {}".format(reg[0], dpinfo_dic[reg[0]]))
                    reg_ls.append(reg)
                    asm_ctg_ls.append(reg[0])
                else:
                    print("{} is normal ctg, dpinfo: {}".format(reg[0], dpinfo_dic[reg[0]]))
        print("Whole ctg denovo have: {}, Apply ctg denovo on: {}".format(whole_ctg_asm_reg_ls, ",".join(asm_ctg_ls)))
        if len(asm_ctg_ls) > 0:
            ctg_denovo_dic = solve_by_denovo.run_ctg_denovo(reg_ls, ctg_denovo_dir, reg_read_ids_dic, asm_region_fq, fa_for_consensus, threads, data_type, config)   # 后面会合并    
        else:
            print("No asm ctg denovo")
            ctg_denovo_dic = {}
    else:
        ctg_denovo_dic = {}
    t4 = time.time()
    
    # 2.2、merge_denovo
    merge_denovo_dir = out_dir +"/merge_denovo"
    make_dir(merge_denovo_dir)
    merge_denovo_fq = merge_denovo_dir + "/merge.fastq"
    merge_reg_fq = merge_denovo_dir + "/reg.fastq"
    merge_reg_read_ids_fn = merge_denovo_dir + "/asm_reg_ids.bed"
    merge_reg_read_ids = set()
    merge_local_reg_ls = []
    merge_reg_read_ids.update(clip_read_ids)
    print(reg_read_ids_dic.keys())
    if len(asm_ctg_ls) == 0:    # 
        merge_reg_read_ids = reg_read_ids_dic["all_read"]
        merge_local_reg_ls = asm_regions
    else:
        for reg in asm_regions:
            if reg[0] in asm_ctg_ls: continue
            reg_id = solve_by_denovo.reg_to_id(reg)
            merge_reg_read_ids.update(reg_read_ids_dic[reg_id]) # reads ids
            merge_local_reg_ls.append(reg)
            # print("read set:", reg_read_ids_dic[reg_id])
    # stats_file = merge_denovo_dir + "/" + "denovo_read.stats"
    # with open(stats_file, "w") as f:
    #     f.write("Asm read num:{}\n".format(len(asm_ids)))
    #     f.write("Clip read num:{}\n".format(len(clip_read_ids)))
    #     f.write("All read num:{}\n".format(len(all_ids)))
    #     f.write("Clip read add num:{}\n".format(len(all_ids)-len(asm_ids)))
    # exit(0)
    logger.info("Apply merge local denovo on: {}".format(merge_local_reg_ls))
    extract_read_fastq_utils.write_read_ids(merge_reg_read_ids_fn, merge_reg_read_ids)
    solve_by_denovo.select_reads_from_names(asm_region_fq, merge_reg_fq, merge_reg_read_ids_fn, threads)
    reads_merge_cmd = ["cat", unmapped_fq, merge_reg_fq, ">", merge_denovo_fq]
    logger.info("Running: %s", " ".join(reads_merge_cmd))
    subprocess.check_call(" ".join(reads_merge_cmd), shell=True)
    # 
    merge_denovo_asm_out = Run_for_denovo(merge_denovo_fq, merge_denovo_dir, threads, genome_size, data_type, config)
    if merge_denovo_asm_out == None:    # 组装失败
        print("Assembly failed, Done")
        exit(0)
    
    t5 = time.time()


    '''## get reads of asm_reagion     may be parallel??
    candidate_op_ls = []
    for op_ls in candidate_op_dic.values():
        candidate_op_ls.extend(op_ls)
    asm_regions, whole_ctg_asm_reg_ls = [], []
    for rec in candidate_op_ls:
        if rec.operation.endswith("asm"):
            asm_regions.append([rec.chr_id, rec.start, rec.end])
            if rec.operation == "whole_ctg_asm":
                whole_ctg_asm_reg_ls.append([rec.chr_id, rec.start, rec.end])
    
    logger.info("Candidate local denovo regions:{}, whole_ctg_asm_reg_ls: {}".format(asm_regions, whole_ctg_asm_reg_ls))
    asm_regions_ids = set()
    reg_read_ids_dic = defaultdict(set)     # asm_reg_id -> reg_read_ids
    # asm_regions_len = 0
    for reg in asm_regions:
        # asm_regions_len += reg[2] - reg[1]
        asm_reg_id = solve_by_denovo.reg_to_id(reg)
        reg_read_ids_dic[asm_reg_id] = set()
        if config["remove_secondary"]:
            for read in bam_reader.fetch(reg[0], reg[1], reg[2]):   # exclude secondary
                if read.is_secondary: continue
                asm_regions_ids.add(read.query_name)
                reg_read_ids_dic[asm_reg_id].add(read.query_name)
        else:
            for read in bam_reader.fetch(reg[0], reg[1], reg[2]):   # 加入了所有的序列
                asm_regions_ids.add(read.query_name)
                reg_read_ids_dic[asm_reg_id].add(read.query_name)
    
    

    t3 = time.time()
    logger.info(" ########################## Extract reads_id finished ########################## ")

    asm_region_fq = denovo_asm_dir + "/asm_region.fastq"
    asm_regions_ids_fn = denovo_asm_dir + "/asm_region_ids.bed"
    with open(asm_regions_ids_fn, "w") as f:
        for read_id in asm_regions_ids: f.write("{}\n".format(read_id))
    ## 读数提取
    select_reads_from_names(fastq_in, asm_region_fq, asm_regions_ids_fn, threads)  # 提取所有需要进行局部组装的区间读数
    
    t4 = time.time()
    logger.info(" ########################## Extract all fastq finished, cost {}s ########################## ".format(time.time() - t3))
    logger.info(" ########################## Extract reads_id cost {}s, Extract fastq cost {}s ########################## ".format(t3 - t2, t4 - t3))

##########################  Start assembly, 3steps  ##########################

    reg_denovo_dir = out_dir + "/" + "denovo_asm1"
    ctg_denovo_dir = out_dir + "/" + "denovo_asm2"
    merge_denovo_dir = out_dir + "/" + "denovo_asm3"
    make_dir(reg_denovo_dir)
    make_dir(ctg_denovo_dir)
    make_dir(merge_denovo_dir)
    merge_reg_read_ids_fn = merge_denovo_dir + "/asm_region_ids.bed"
    merge_reg_fq = merge_denovo_dir + "/reg.fastq"
    merge_denovo_fq = merge_denovo_dir + "/denovo.fastq"
    ## 2.1、run reg denovo
    if config["apply_denovo_by_reg"]:   # if choose the option
        semi_candidate_op_dic, failed_reg_ls, success_reg_ls = solve_by_denovo.run_reg_denovo(asm_regions, reg_denovo_dir, reg_read_ids_dic, candidate_op_dic, asm_region_fq, \
        bam_in, threads, data_type, fa_for_consensus, depthrec_dic, config)
    else:   # 跳过reg denovo
        semi_candidate_op_dic = candidate_op_dic
        failed_reg_ls = asm_regions
        success_reg_ls = []
    t4_1 = time.time()
    print("Reg denovo assembly Done, cost {}s !!!".format(time.time() - t4))
    ## 2.2 run ctg denovo for abnormal whole ctg asm
    print("Start ctg denovo")
    asm_ctg_ls = []
    if config["apply_ctg_denovo"]:
        whole_ctg_asm_reg_ls = solve_by_denovo.filter_empty_ctg(whole_ctg_asm_reg_ls, reg_read_ids_dic)
        reg_ls = []
        if config["denovo_by_ctg"]["apply_all"]:
            print("Apply on all denovo ctg")
            asm_ctg_ls = [reg[0] for reg in whole_ctg_asm_reg_ls]
            reg_ls = whole_ctg_asm_reg_ls
        else:
            for reg in whole_ctg_asm_reg_ls:
                if solve_by_denovo.is_abnormal_ctg(dpinfo_dic[reg[0]], config["denovo_by_ctg"]):
                    print("{} is abnormal ctg, dpinfo: {}".format(reg[0], dpinfo_dic[reg[0]]))
                    reg_ls.append(reg)
                    asm_ctg_ls.append(reg[0])
                else:
                    print("{} is normal ctg, dpinfo: {}".format(reg[0], dpinfo_dic[reg[0]]))
        print("Whole ctg denovo have: {}, Apply ctg denovo on: {}".format(whole_ctg_asm_reg_ls, ",".join(asm_ctg_ls)))
        if len(asm_ctg_ls) > 0:
            ctg_denovo_dic = solve_by_denovo.run_ctg_denovo(reg_ls, ctg_denovo_dir, reg_read_ids_dic, asm_region_fq, fa_for_consensus, threads, data_type, config)   # 后面会合并    
        else:
            print("No asm ctg denovo")
            ctg_denovo_dic = {}
    else:
        ctg_denovo_dic = {}
    t4_2 = time.time()
    ## 2.3 run merge local denovo assembly
    logger.info("Start merge local denovo")
    merge_reg_read_ids = set()
    merge_local_reg_ls = []
    for reg in failed_reg_ls:
        if reg[0] in asm_ctg_ls: continue
        reg_id = solve_by_denovo.reg_to_id(reg)
        merge_reg_read_ids.update(reg_read_ids_dic[reg_id]) # reads ids
        merge_local_reg_ls.append(reg)     # 添加数据
    with open(merge_reg_read_ids_fn, "w") as f:
        for read_id in merge_reg_read_ids: 
            f.write("{}\n".format(read_id))
    select_reads_from_names(asm_region_fq, merge_reg_fq, merge_reg_read_ids_fn, threads)
    logger.info("Apply merge local denovo on: {}".format(merge_local_reg_ls))
    # 读数合并，unmapped_fq + merge_reg_fq -> final_denovo_fq
    reads_merge_cmd = ["cat", unmapped_fq, merge_reg_fq, ">", merge_denovo_fq]
    logger.info("Running: %s", " ".join(reads_merge_cmd))
    subprocess.check_call(" ".join(reads_merge_cmd), shell=True)
    # asm
    if out_denovo:  # For test
        print("Skip merge denovo, provide asm from out:", out_denovo)
        merge_denovo_asm_out = out_denovo
    else:
        merge_denovo_asm_out = Run_for_denovo(merge_denovo_fq, merge_denovo_dir, threads, genome_size, data_type, config)
        if merge_denovo_asm_out == None:    # 组装失败
            pass
    t5 = time.time()
    logger.info("Merge denovo assembly Done, cost {}s !!!".format(time.time() - t4_1))
    logger.info("Local assembly done, Time stats:\nreg cost {}s\nctg cost {}s\nmerge cost {}s\nsum: {}s".format(t4_1 - t4, t4_2 - t4_1, t5 - t4_2, t5 - t4))
    '''
##########################   Assembly Done   ##########################
    ## 3、map asm fasta to reference  ## 命令还要改一下
    logger.info("****************** map denovo asm to reference start ******************")
    fa_for_consensus_dic = fasta_parser.read_sequence_dict(fa_for_consensus)
    query_asm = merge_denovo_asm_out    # denovo asm out
    sorted_bam_out = merge_denovo_dir + "/denovoasm_to_ref.sorted.bam"
    minimap2_cmd_ls = ["minimap2", "-ax", "asm20", "-t", str(threads), fa_for_consensus, query_asm, "|", "samtools", "sort", "-O", "BAM", "-@", str(threads), "-o", sorted_bam_out] 
    index_cmd_ls = ["samtools index -@", str(threads), sorted_bam_out]
    run_cmd_ls(minimap2_cmd_ls)
    run_cmd_ls(index_cmd_ls)
    logger.info("******************map denovo asm to reference Done!!!******************")
    t6 = time.time()
    
    ### ****************** [SV consensus on reference] ******************
    ## scan regions to determine specific opertion and consensus on reference
    # SV_consensus_dir = out_dir + "/SV_consensus"
    consensus_fa_dic = {}
    SV_consensus_dir = out_dir
    make_dir(SV_consensus_dir)
    SVconsensus_bed_out = SV_consensus_dir + "/consensus.bed"
    consensus_fasta_out = SV_consensus_dir + "/consensus.fasta"
    asm_to_ref_bam = sorted_bam_out
    denovo_fa = merge_denovo_asm_out    # denovo asm out
    asm_fa_dic = fasta_parser.read_sequence_dict(denovo_fa)

    # run
    logger.info("Run SV_consensus_on_ref")
    if config["apply_scaffod"]:
        print("-----------------------apply_scaffod mode-----------------------")
        connect_info_ls, patch_ids, consensus_fa_dic = SV_consensus_on_ref2(all_chrs, semi_candidate_op_dic, fa_for_consensus_dic, asm_fa_dic, asm_to_ref_bam, SVconsensus_bed_out, consensus_fasta_out, out_dir, threads, config)

    else:
        # pool = Pool(processes=threads)
        # results = [pool.apply_async(SV_consensus_on_ref, args=(ctg, semi_candidate_op_dic[ctg], fa_for_consensus_dic, asm_fa_dic, asm_to_ref_bam, config)) for ctg in all_chrs]
        # # results = [pool.apply_async(SV_consensus_on_ref, args=(ctg, semi_candidate_op_dic[ctg], ref_dic, asm_fa_dic, asm_to_ref_bam, Nfill_size)) for ctg in all_chrs]
        # pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
        # pool.join() # 等待进程池中的所有进程执行完毕
        # connect_info_ls = []
        # final_rec_ls = []
        # for i, res in enumerate(results):
        #     ctg_final_rec_ls, connect_info, ctg_consensus_fa_dic = res.get()
        #     # ctg = all_chrs[i]
        #     connect_info_ls.append(connect_info)
        #     consensus_fa_dic.update(ctg_consensus_fa_dic)   # 将每个线程的字典 consensus 结果加入进来
        #     final_rec_ls.extend(ctg_final_rec_ls)
        # Record.write_record(final_rec_ls, SVconsensus_bed_out)
        # fasta_parser.write_fasta_dict(consensus_fa_dic, consensus_fasta_out)
        # 
        # ## patch_ids for next step
        # patch_ids = set()
        # with open(SVconsensus_bed_out, "r") as f:
        #     for line in f:
        #         if line.startswith("#"): continue
        #         fields = line.strip().split("\t")   # 注意加上\t，兼容del
        #         op_ls = fields[3].split(",")
        #         patch_ls = fields[5].split(",")     # 
        #         for i, op in enumerate(op_ls):
        #             if op == "replace_with_denovo":
        #                 for id in patch_ls:patch_ids.add(id)    #   replace这种也要排除
        #             if op == "asm_patch":   # 用掉的碎片ctg
        #                 patch_ids.add(patch_ls[i])
        pass
    logger.info("consensus fasta write Done !!! ")
    logger.info("Run SV_consensus_on_ref Done !!! ")
    t7 = time.time()

    ## collect denovo unmapped assembly 收集组装未比对上的读数，这部分读数(asm)一般是参考上所缺失的
    # denovo_unmapped_fa = SV_consensus_dir + "/denovo_unmapped.fasta"
    # get_unmapped_reads(threads, asm_to_ref_bam, denovo_unmapped_fa, "fasta")    # 通过工具将未比对上的组装结果提取出来
    # denovo_unmapped_dic = fasta_parser.read_sequence_dict(denovo_unmapped_fa)   # 将其读取到字典中吧
    # logger.info("denovo unmapped fasta write Done !!!")

    ## collect denovo fragemants
    denovo_fragments_fa = SV_consensus_dir + "/denovo_fragments.fasta"  # 组装碎片，主要指的是未用于填充修补的denovo结果
    print(patch_ids)
    denovo_fragments_dic = get_frag_reads(denovo_fa, patch_ids) # 收集组装碎片
    fasta_parser.write_fasta_dict(denovo_fragments_dic, denovo_fragments_fa)
    logger.info("denovo fragements write Done !!!")

    ## write merge fasta
    merge_fa_out = SV_consensus_dir + "/merge.fasta"
    merge_fa_dic = consensus_fa_dic
    if not ex_unmapped_denovo:  # 将碎片保留下来, default
        # merge_fa_dic = copy.deepcopy(consensus_fa_dic)  # 先拷贝过来把
        denovo_merge_num = 0
        min_denovo_to_keep = config["min_denovo_to_keep"]  # 设置组装碎片最小保留长度
        for ctg,seq in denovo_fragments_dic.items():
            if len(seq) > min_denovo_to_keep:
                merge_fa_dic[ctg + "_fragemant"] = seq
                denovo_merge_num += 1 
                ## 记录连接信息吧
                seq_id = ctg + "_fragemant"
                connect_info_ls.append(Connect_info(seq_id, [seq_id], []))
        logger.info("merge {} fragements from denovo".format(denovo_merge_num))
        fasta_parser.write_fasta_dict(merge_fa_dic, merge_fa_out)
    else:   # 不需要合并，直接保留consensus结果
        # merge_fa_dic = consensus_fa_dic
        logger.info("keep same with consensus fasta")
        # shutil.copy(consensus_fasta_out, merge_fa_out)
    ## 保持不变的contig
    # for chr in keep_ls:
    #     merge_fa_dic[chr] = ref_dic[chr]

    ##
    merge_fa_dic.update(ctg_denovo_dic)     # 合并来自ctg asm的数据
    fasta_parser.write_fasta_dict(merge_fa_dic, merge_fa_out)
    connect_info_file = os.path.join(SV_consensus_dir, "connect.txt")
    Connect_info.write_connect_info(connect_info_ls, connect_info_file)
    logger.info("final merge fasta write Done !!!")
    t8 = time.time()

    print("Time stastic of run_SVconsensus_parallel:")
    time_ls = [t1, t2, t3, t4, t5, t6, t7, t8]
    for i in range(len(time_ls) - 1):
        print("t{}-t{}: {}s".format(str(i+1), str(i+2), str(time_ls[i+1]-time_ls[i])))

def main():
    ## for debug
    # out_dir = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi_ctgs/my_pipe/step3_SV_consensus"
    out_dir = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/create_consensus_by_bed/test_consensus"
    bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi_ctgs/NC_19/my_pipe/step1_mapping/aln.sorted.bam"
    process_ctg_ls = ["NC_000019.10"]
    # candidate_bed = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi_ctgs/NC_19/my_pipe/step2_candidate_regions/candidate.bed"
    candidate_bed = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/find_candidate_regions/test_find/candidate.bed"
    threads = 20
    fastq_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/chm13/hifi/NC_060943.1.fastq"
    data_type = "hifi"
    reference_fn = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/chm13/REF/ref1/chrs/NC_000019.10.fasta"
    Nfill_size = 20
    genome_size = 5800000
    ex_unmapped_denovo = False
    config = {}
    run_SVconsensus_parallel(out_dir, bam_in, process_ctg_ls, candidate_bed, 
                             threads, fastq_in, data_type, reference_fn, 
                             Nfill_size, genome_size, ex_unmapped_denovo, config)
    return
    
    ### ****************** [Argument Parser] ******************
    parser = argparse.ArgumentParser(description="get_SV_consensus")
    parser.add_argument("-t", "--threads", dest="threads", type=int, default=1)
    parser.add_argument("--ref", dest="ref", required=True)
    parser.add_argument("--fastq", dest="fastq_in", required=True, help="fastq reads file")
    parser.add_argument("--bam", dest="bam_in", required=True, help="bam file")
    parser.add_argument("--out-dir", dest="out_dir", required=True)
    parser.add_argument("--contig_ls", dest="ctg_ls", required=True, help="format like:chr1,chr2")
    parser.add_argument("--bed", dest="candidate_bed", required=True, help="provide candidate bed, will get record from this")
    parser.add_argument("--min-contig", dest="min_contig", help="skip process contig shorter than this, keep with raw", default=200000, type=int)   # 调 1,000,000     5,000,000
    parser.add_argument("--data-type", dest="data_type", required=True, choices=["ont", "hifi"], help="fastq file type")
    parser.add_argument("--fill-size", dest="Nfill_size", type=int, default=20, help="N_fill_size")
    parser.add_argument("-a", dest="all_chrs", action='store_true')
    # group = parser.add_mutually_exclusive_group()
    # group.add_argument('--preset1', action='store_true')   
    # group.add_argument('--preset2', action='store_true')

    args = parser.parse_args()
    bam_reader = pysam.AlignmentFile(args.bam_in, "rb")
    if args.all_chrs:
        process_ctg_ls = list(bam_reader.references)
    else:
        process_ctg_ls = [ctg for ctg in args.ctg_ls.split(",") if bam_reader.get_reference_length(ctg) > args.min_contig]
    logger.info("Process contigs: {} !!!".format(process_ctg_ls))
    # if args.preset1:    ## 配置默认参数
    #     parser.set_defaults(data_type="ont")
    # elif args.preset2:
    #     parser.set_defaults()
    run_SVconsensus_parallel(args)

if __name__ == "__main__":

    # from create_consensus_by_bed import solve_by_denovo
    # dpinfo_file = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/yeast_ont/my_pipe3/step2_candidate_regions/dp_info/dpinfo.bed"
    # dpinfo_dic = Depth_info.read_dp_info(dpinfo_file)
    # print(dpinfo_dic)
    # main()
    # return
    # candidate_bed = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi/my_pipe/step2_candidate_regions/candidate.bed"
    # out_dir = "/public/home/hpc214712170/test/te"
    # process_ctg_ls = ['NC_000001.11', 'NC_000002.12', 'NC_000003.12', 'NC_000004.12', 'NC_000005.10', 'NC_000006.12', 'NC_000007.14', 'NC_000008.11', 'NC_000009.12', 'NC_000010.11', 'NC_000011.10', 'NC_000012.12', 'NC_000013.11', 'NC_000014.9', 'NC_000015.10', 'NC_000016.10', 'NC_000017.11', 'NC_000018.10', 'NC_000019.10', 'NC_000020.11', 'NC_000021.9', 'NC_000022.11', 'NC_000023.11', 'NC_000024.10']
    # threads = 40
    # bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi/my_pipe/step1_mapping/aln.sorted.bam"
    # t1 = time.time()

    # bedinfo_dic = load_candidate_reg_bed(candidate_bed)
    # ## parallel process ctgs   
    # candidate_op_dir = out_dir + "/candidate_op"
    # if not os.path.exists(candidate_op_dir):
    #     os.makedirs(candidate_op_dir)
    # all_reference_ids = process_ctg_ls
    # random.shuffle(all_reference_ids)
    # bed_to_merge = [candidate_op_dir + "/candidate_op_" + ctg + ".bed" for ctg in all_reference_ids]
    # # pool = Pool(processes=threads)
    # # results = [pool.apply_async(candidate_process, args=(bedinfo_dic.get(all_reference_ids[id], []), bam_in, bed_to_merge[id])) for id in range(len(all_reference_ids))]
    # # pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    # # pool.join() # 等待进程池中的所有进程执行完毕
    # # '''注意：apply_async的话需要使用get来得到返回值'''
    # res = []
    # for id in range(len(all_reference_ids)):
    #     tmp_t = time.time()
    #     candidate_process(bedinfo_dic.get(all_reference_ids[id], []), bam_in, bed_to_merge[id])
    #     print(all_reference_ids[id], "Cost: ", time.time() - tmp_t)
    # candidate_op_dic = defaultdict(list)   # 将所有染色体都加进去，便于统一
    # print("All Cost: ", time.time()-t1)
    pass
