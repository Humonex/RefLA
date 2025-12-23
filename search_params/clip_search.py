import pysam
import re
from collections import namedtuple
import numpy as np
from multiprocessing import Pool
import search_pileup_params
Region = namedtuple('Region', ["chr_id", "start", "end"])

def cal_sum(arry, l, r):
    return sum(arry[l:r])
def stats_ls(data_ls):
    
    return

def find_from_clip(bam_in, ctg, MIN_CLIP_LEN, win_size):  # find candidate regions by clips
    # region_ls_merge = []    # 保存所有ctg的region
    bam_index = bam_in + ".bai"
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_index)
    ctg_len = bam_reader.get_reference_length(ctg)  # 172126628
    reg_start = 0
    reg_end = ctg_len

    print("Find {}:{}-{}".format(ctg, reg_start, reg_end))
    clip_ls = [0] * ctg_len # 记录contig每个位置>min_clip_len的clip数目``
    MIN_MAPPING_QUALITY = 10    # 10
    for read in bam_reader.fetch(ctg, start=reg_start, end=reg_end):
        if read.is_unmapped or read.is_secondary or read.mapping_quality < MIN_MAPPING_QUALITY:
            continue
        ref_pos = read.reference_start
        cigar = read.cigarstring
        tokens = re.findall("[\d]{0,}[A-Z]{1}", cigar)
        # left
        left = tokens[0]
        if left[-1] in "HS":
            if int(left[:-1]) > MIN_CLIP_LEN:
                # clip_ls.append([read.reference_start, int(left[:-1])]) # 记录下clip的位置，长度
                clip_ls[read.reference_start] += 1
        
        # right
        right = tokens[-1]
        if right[-1] in "HS":
            if int(right[:-1]) > MIN_CLIP_LEN:
                # clip_ls.append([read.reference_end, int(left[:-1])]) # 记录下clip的位置，长度
                # if read.reference_end > len(clip_ls)-1:
                #     print(read.query_name, "\t"+ctg+":"+str(read.reference_start)+"-"+str(read.reference_end))
                clip_ls[read.reference_end-1] += 1

    print(ctg+":"+"FIND clip_regions")
    win_size = 1000  # 200 400
    stride = 500
    # MIN_CLIP_NUM = 5
    win_num = ctg_len // stride + 1 # 
    win_clip_num = [0] * win_num    # clip num of a window

    for i in range(win_num):    # i*stride, i*stride+win_size
        win_clip_num[i] = cal_sum(clip_ls, i*stride, i*stride+win_size) if i*stride+win_size < ctg_len else cal_sum(clip_ls, i*stride, ctg_len)
    return ctg, win_clip_num # 一个数组记录


def stats_clip(bam_in, ctg_ls, MIN_CLIP_LEN, out_file, threads=20):
    print("Stats clip info")
    win_clip_num = []
    dic = {}
    '''for ctg in ctg_ls:
        ctg_win_clip_num_ls = find_from_clip(bam_in, ctg, MIN_CLIP_LEN, 1000)
        win_clip_num.extend(ctg_win_clip_num_ls)
        dic[ctg] = ctg_win_clip_num_ls'''
    pool = Pool(processes=threads)
    results = [pool.apply_async(find_from_clip, args=(bam_in, ctg, MIN_CLIP_LEN, 1000)) for ctg in ctg_ls]
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    
    for res in results:
        ctg, ctg_win_clip_num_ls = res.get()
        dic[ctg] = ctg_win_clip_num_ls
        win_clip_num.extend(ctg_win_clip_num_ls)
    ## 
    if len(win_clip_num) == 0:
        print("----------null win_clip_num ls, check again!!!----------")
    win_clip_num.sort()
    # print(win_clip_num)
    win_clip_num = [float(a) for a in win_clip_num]
    win_clip_num = np.array(win_clip_num)

    ### stats
    search_pileup_params.stats_dic({'clip':win_clip_num}, {'clip':False}, out_file)
    '''print("median: ", np.median(win_clip_num))
    step = 0.01
    quantile_ls = [i * step for i in range(int(1 / step))]
    quantile_num_ls = np.quantile(win_clip_num, quantile_ls)
    for i in range(len(quantile_ls)):
        print("{:3f}:{:.4f}".format(quantile_ls[i], quantile_num_ls[i]))
        print(quantile_ls[i], ":", quantile_num_ls[i])
'''
if __name__ == "__main__":
    import sys
    '''bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/thaliana_ont/tests/my_pipe4/step1_mapping/aln.sorted.bam"
    ctg_ls = ["LR699760.1","LR699761.1","LR699762.1","LR699763.1","LR699764.1"]'''
    bam_in = sys.argv[1]
    out_file = sys.argv[2]
    MIN_CLIP_LEN = 500
    bam_reader = pysam.AlignmentFile(bam_in, "rb")
    ctg_ls = list(bam_reader.references)
    stats_clip(bam_in, ctg_ls, MIN_CLIP_LEN, out_file)
# [0.05, 0.25, 0.5, 0.75, 0.8, 0.85, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99]:  [ 0.  0.  0.  1.  1.  2.  2.  2.  3.  3.  4.  6.  8. 11. 18. 28.]
