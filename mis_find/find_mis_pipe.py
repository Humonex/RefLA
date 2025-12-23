
import glob
import logging
import math
import signal
import sys
sys.path.append("/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe")
import os
from multiprocessing import Pool
import re
import subprocess
import numpy as np
import yaml
import pysam
from collections import defaultdict, Counter
from find_candidate_regions.find_reg_by_depth import Depth_info, get_dp, read_mosdepth_dp_file, get_dp_info_parallel
# import collect_signature, cluster_sigs    # vscode
import mis_find.collect_signature, mis_find.cluster_sigs
def make_dir(dir):
    if not os.path.isdir(dir):
        os.makedirs(dir)
def cal_sum(arry, l, r):
    return sum(arry[l:r])

class Mis_Info():   # mis_info
    def __init__(self, ctg, start, end, ) -> None:
        self.ctg = ctg
        self.start = start
        self.end = end
        # self.avg_dp = 
        # self.correct_portion
        # self.disagree_portion = 
        # self.differ_portion = 
        pass

class Alarm(Exception):
    pass

def alarm_handler(signum, frame):
    raise Alarm

def exe(cmd, timeout=-1):
    """
    Executes a command through the shell.
    timeout in minutes! so 1440 mean is 24 hours.
    -1 means never
    """
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, \
                            stderr=subprocess.STDOUT, close_fds=True,\
                            preexec_fn=os.setsid)
    signal.signal(signal.SIGALRM, alarm_handler)
    if timeout > 0:
        signal.alarm(int(timeout*60))  
    try:
        stdoutVal, stderrVal =  proc.communicate()
        signal.alarm(0)  # reset the alarm
    except Alarm:
        logging.error(("Command was taking too long. "
                       "Automatic Timeout Initiated after %d" % (timeout)))
        os.killpg(proc.pid, signal.SIGTERM)
        proc.kill()
        return 214,None,None
    
    retCode = proc.returncode
    return retCode,stdoutVal,stderrVal

def call_back(res):
    print(res)

def error_call_back(error_code):
    print(error_code)

def clear_bed_files(directory):
    # 构建文件匹配模式
    pattern = os.path.join(directory, '*.bed')
    # 使用glob找出所有匹配的文件
    files = glob.glob(pattern)
    if not files:
        print(f"No '.bed' files found in directory {directory}")
        return
    for file in files:
        try:
            os.remove(file)
            print(f"Removed file: {file}")
        except OSError as e:
            print(f"Error: {e.strerror}, while trying to remove file: {e.filename}")

def cluster_by_dis(reg_ls_in, dis): # 根据距离进行聚类
    if len(reg_ls_in) <= 1:
        return reg_ls_in
    reg_ls_in = sorted(reg_ls_in, key=lambda x: x[1])   # 聚类之前先排序
    reg_start = -1
    reg_end = -1
    chr_id = reg_ls_in[0][0]
    reg_ls_out = []
    for reg in reg_ls_in:
        if reg_start > -1:  # 非首次
            if reg[1] - reg_end <= dis:
                reg_end = reg[2]
                need_to_cluster.append(reg)
            else:   # new_reg
                reg_ls_out.append([chr_id, reg_start, reg_end])
                reg_start = reg[1]
                reg_end = reg[2]
                need_to_cluster = [reg]
        else:
            reg_start = reg[1]
            reg_end = reg[2]
            need_to_cluster = [reg]
    if reg_start > -1:
        reg_ls_out.append([chr_id, reg_start, reg_end])
    return reg_ls_out

def cal_idels(read:pysam.AlignedRead, reg_start, reg_end):
    '''计算ref_start-ref_end之间的indels'''
    ref_pos = read.reference_start
    cigar = read.cigartuples
    ins_num = 0
    del_num = 0
    sum_ins_num = 0
    sum_del_num = 0
    query_pos = 0
    match_num = 0
    min_size = 30
    
    for op, op_len in cigar:
        if ref_pos > reg_end: break
        if op == 0:
            if op_len + ref_pos < reg_start or ref_pos > reg_end:
                pass
            elif op_len + ref_pos < reg_end:
                match_num += op_len
            else:
                match_num += reg_end - ref_pos
            query_pos += op_len
            ref_pos += op_len
        elif op == 1:   # ins
            query_pos += 1
            if ref_pos > reg_start and ref_pos < reg_end:
                # print("Add ins")
                if op_len > min_size:
                    ins_num += op_len
            sum_ins_num += op_len
        elif op == 2:   # del
            if op_len + ref_pos < reg_start:
                pass
            elif op_len + ref_pos < reg_end:    # []
                # print("Add del")
                if op_len > min_size:
                    del_num += op_len
            else:   # > reg_end
                # print("Add del")
                if op_len > min_size:
                    del_num += op_len
            ref_pos += op_len
            sum_del_num += op_len
        elif op == 4:   # soft clip
            query_pos += op_len
        else:
            continue
    # print(ins_num, del_num)
    return ins_num, del_num

def cal_cov_fraction(dp_ls, threshold):
    size = len(dp_ls)
    if size < 1: return 1
    num = 0
    for dp in dp_ls:
        if dp < threshold:num += 1
    return num / size


'''min_span_num = 2
    min_supp_portion = 0.4  # 有多少成的读数表明没有mis
    MIN_MAPPING_QUALITY = 20    # 10    map_Q太严了，导致大量的比对过滤掉了->3
    MIN_ALIGN_LENGTH = 10000
    MIN_ALIGN_RATE = 0.95
    ins_threshold = 100
    del_threshold = 100
'''
def win_check(ctg, start, end, ctg_len, bam_reader, params):
    print("Check:{}:{}-{}".format(ctg, start, end))
    '''
    True: mis
    False: no mis
    1、对深度的限制，一些重复区不一定能从cigar检测出mis。或者存在干扰
    2、对
    hifi使用不同的参数
    '''
    min_span_num = params['min_span_num']
    min_supp_portion = params['min_supp_portion']  # 有多少成的读数表明没有mis
    min_MQ = params['min_MQ']    # 10    map_Q太严了，导致大量的比对过滤掉了->3
    MIN_ALIGN_LENGTH = params['MIN_ALIGN_LENGTH']
    MIN_ALIGN_RATE = params['MIN_ALIGN_RATE']
    ins_threshold = params['ins_threshold']
    del_threshold = params['del_threshold']
    indel_threshold = ins_threshold + del_threshold
    min_clip_len = params['min_clip_len']
    bound_span = params['bound_span']   # 5000/10000
    # 

    span_ls = []
    for read in bam_reader.fetch(ctg, start, end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < min_MQ:
            continue
        # for bound
        is_bound = False
        if read.reference_start <= max(500, bound_span) or read.reference_end >= min(ctg_len - 500, ctg_len - bound_span):
            MIN_ALIGN_LENGTH = MIN_ALIGN_LENGTH - 2000
            MIN_ALIGN_RATE = MIN_ALIGN_RATE - 0.1
            is_bound = True
        ## 
        if read.query_alignment_length < MIN_ALIGN_LENGTH or (read.query_alignment_length / read.query_length) < MIN_ALIGN_RATE:
            continue
        cigar = read.cigartuples
        left = cigar[0]
        if left[0] == 4 or left[0] == 5:
            if left[1] > min_clip_len and not is_bound:    # left clip非左边界，且其长度大于阈值
                continue
        right = cigar[-1]
        if right[0] == 4 or right[0] == 5:
            if right[1] > min_clip_len and not is_bound: # right clip非右边界，且其长度大于阈值
                continue
        if read.reference_start <= max(0, start - 500) and read.reference_end >= min(ctg_len, end + 500):
            span_ls.append(read)
    
    #  
    if len(span_ls) < min_span_num and start >= 5000 and end <= ctg_len - 5000:
        # filtered_ls.append([reg[0], reg[1], reg[2], "low_supp_mis"])
        return True, "low_supp_mis"
    else:
        min_span_num = max(min_span_num, math.ceil(min_supp_portion*len(span_ls)))
        # cal indel of the span
        Inss_ls = []
        Dels_ls = []
        Indels_ls = []
        for read in span_ls:   # 收集span_ls读数的indels
            read_l = max(start - 5000, read.reference_start + 100)
            read_r = min(end + 5000, read.reference_end - 100)
            indels = cal_idels(read, read_l, read_r)
            # indels = cal_idels(read, start, end)    # max(0, start - 500), min(ctg_len, end + 500)
            Indels_ls.append(indels[0] + indels[1])
            Inss_ls.append(indels[0])
            Dels_ls.append(indels[1])
            # if Indels_ls[-1] == 0:
            #     print(read.query_name)
        # 
        Inss_ls.sort()
        Dels_ls.sort()
        Indels_ls.sort()
        # check indels
        avg_ins = sum(Inss_ls[:min_span_num]) // min_span_num
        avg_del = sum(Dels_ls[:min_span_num]) // min_span_num
        avg_indels = sum(Indels_ls[:min_span_num]) // min_span_num
        # print(Inss_ls, Dels_ls, Indels_ls)
        if avg_ins < ins_threshold and avg_del < del_threshold and avg_indels < indel_threshold: # no mis
            return False, "No_mis"
        else:
            return True, "Reads_mis"

def filter(bam, ctg, reg_start, reg_end, ctg_len, clu_ls, out_dir):
    bam_reader = pysam.AlignmentFile(bam, "rb", index_filename=bam + ".bai")
    reg_id = ctg + ":" + str(reg_start) + "-" + str(reg_end)
    print("Start classify and filter: {}:{}-{}".format(ctg, reg_start, reg_end))
    filtered_dir = os.path.join(out_dir, "filtered")
    filtered_bed = os.path.join(filtered_dir, reg_id + ".bed")
    # 
    # low_dp_threshold = 3
    # low_cov_fraction = 0.1
    # min_clip_ratio = 0.4
    MIN_MAPPING_QUALITY = 20    # 10
    boundary = 5000
    min_span_num = 3
    min_supp_portion = 0.2
    MIN_ALIGN_LENGTH = 10000
    MIN_ALIGN_RATE = 0.95
    ins_threshold = 100
    del_threshold = 100
    filtered_ls = []
    out_ls = []
    # F1:
    for reg in clu_ls:
        ctg, start, end = reg
        # filter boundary
        # if (start < 1000 and end - start < boundary) or (end > ctg_len - 1000 and end - start < boundary): 
        #     out_ls.append([reg, "boundary"])
        #     continue
        if start < 5000:
            start = 5000
        if end > ctg_len - 5000:
            end = ctg_len - 5000
        if end <= start or start > ctg_len or end <= 0:
            out_ls.append([reg, "bounder"])
            continue
        # 
        '''start_idx = start // stride
        end_idx = (end - win_size) // stride if end < reg_end else win_num - 1
        reg_dp_ls = [round(win_dp_ls[i], 4) for i in range(start_idx, end_idx + 1)]
        reg_clip_ls = [round(win_clip_num[i], 4) for i in range(start_idx, end_idx + 1)]
        avg_dp = np.mean(reg_dp_ls)
        cov_fraction = cal_cov_fraction(reg_dp_ls, low_dp_threshold)'''
        # Type 1
        if reg[2] - reg[1] > 50000:
            filtered_ls.append([reg[0], reg[1], reg[2], "Large_mis"])
            continue
        # 
        '''left_idx = max(start - 4000, 0) // stride
        right_idx = (end + 4000 - win_size) // stride if (end + 4000) < reg_end else win_num - 1
        surround_dp_ls = []
        for i in range(left_idx, start_idx):
            surround_dp_ls.append(win_dp_ls[i])
        for i in range(end_idx, right_idx):
            surround_dp_ls.append(win_dp_ls[i])
        surround_dp = np.median(surround_dp_ls)
        clip_ratio = max(reg_clip_ls) // surround_dp if surround_dp != 0 else 0'''
        span_ls = []
        for read in bam_reader.fetch(ctg, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < MIN_MAPPING_QUALITY or read.query_alignment_length < MIN_ALIGN_LENGTH or (read.query_alignment_length / read.query_length) < MIN_ALIGN_RATE:
                continue
            if read.reference_start < max(0, start - 500) and read.reference_end > min(ctg_len, end + 500):
                span_ls.append(read)
        #  
        if len(span_ls) < min_span_num:
            filtered_ls.append([reg[0], reg[1], reg[2], "low_supp_mis"])
        else:
            min_span_num = max(min_span_num, math.ceil(min_supp_portion*len(span_ls)))
            # cal indel of the span
            Inss_ls = []
            Dels_ls = []
            for i in range(len(span_ls)):   # 收集span_ls读数的indels
                read = span_ls[i]
                indels = cal_idels(read, start, end)    # max(0, start - 500), min(ctg_len, end + 500)
                Inss_ls.append(indels[0])
                Dels_ls.append(indels[1])
            # 
            Inss_ls.sort()
            Dels_ls.sort()
            # check indels
            avg_ins = sum(Inss_ls[:min_span_num]) // min_span_num
            avg_del = sum(Dels_ls[:min_span_num]) // min_span_num
            if avg_ins < ins_threshold and avg_del < del_threshold: # no mis
                out_ls.append([reg, "No_mis"])
            else:
                filtered_ls.append([reg[0], reg[1], reg[2], "Medium_mis"])
    
    ##
    print("filtered_ls:", filtered_ls)
    with open(filtered_bed, "w") as f:
        for ele in filtered_ls:
            f.write("{}\t{}\t{}\t{}\n".format(ele[0], ele[1], ele[2], ele[3]))
    print("{}:{}-{} out_ls:{}".format(ctg, reg_start, reg_end, out_ls))
    print("Stop classify and filter")

def filter2(bam, ctg, reg_start, reg_end, ctg_len, clu_ls, out_dir, config):
    ## classify and filter the reg
    '''
    过滤：
    1、边界区域，跳过5000bp
    2、过大的区域。采用滑动窗口的过滤方式
    3、过滤出错区域 //使用dp、sim、diff、dis指标来跳过这类区域
    1) 高倍重复区
    2) 嘈杂的区域
    '''
    bam_reader = pysam.AlignmentFile(bam, "rb", index_filename=bam + ".bai")
    reg_id = ctg + ":" + str(reg_start) + "-" + str(reg_end)
    print("Start classify and filter2: {}:{}-{}".format(ctg, reg_start, reg_end))
    filtered_dir = os.path.join(out_dir, "filtered2")
    filtered_bed = os.path.join(filtered_dir, reg_id + ".bed")
    # 
    params = config["filter2"]
    # low_dp_threshold = 3
    # low_cov_fraction = 0.1
    # min_clip_ratio = 0.4
    # min_span_num = 3
    # MIN_ALIGN_LENGTH = 10000
    # MIN_ALIGN_RATE = 0.95
    # ins_threshold = 100
    # del_threshold = 100
    filtered_ls = []
    out_ls = []
    # F2：滑动窗口式的过滤策略
    # check_win_size = 10000   # 5000
    check_win_size = params['check_win_size']
    bound_span = params['bound_span']   # 5000/10000
    for reg in clu_ls:
        ctg, start, end = reg
        flag = False
        start = bound_span if start < bound_span else start
        end = ctg_len - bound_span if end > ctg_len - bound_span else end
        if start > end: 
            print("Skip:", reg)
            continue
        if end - start <= check_win_size:
            res = win_check(ctg, start, end, ctg_len, bam_reader, params)
            if res[0]:
                filtered_ls.append([ctg, start, end, res[1], [start, end]])
                flag = True # has mis
        else:
            # slide check
            win_start = start
            win_end = start + check_win_size
            while win_start <= end and win_end - win_start > 0:
                res = win_check(ctg, win_start, win_end, ctg_len, bam_reader, params)
                if res[0]:
                    filtered_ls.append([ctg, start, end, res[1], [win_start, win_end]])
                    flag = True # has mis
                    break
                win_start += check_win_size
                win_end = min(end, win_end + check_win_size)
        if not flag:    # no mis
            out_ls.append([reg, "No_mis"])
    # F1:
    # for reg in clu_ls:
    #     ctg, start, end = reg
    #     # filter boundary
    #     if (start < 1000 and end - start < boundary) or (end > ctg_len - 1000 and end - start < boundary): 
    #         out_ls.append([reg, "boundary"])
    #         continue
        
    #     # 
    #     '''start_idx = start // stride
    #     end_idx = (end - win_size) // stride if end < reg_end else win_num - 1
    #     reg_dp_ls = [round(win_dp_ls[i], 4) for i in range(start_idx, end_idx + 1)]
    #     reg_clip_ls = [round(win_clip_num[i], 4) for i in range(start_idx, end_idx + 1)]
    #     avg_dp = np.mean(reg_dp_ls)
    #     cov_fraction = cal_cov_fraction(reg_dp_ls, low_dp_threshold)'''
    #     # Type 1
    #     if reg[2] - reg[1] > 50000:
    #         filtered_ls.append([reg[0], reg[1], reg[2], "Large_mis"])
    #         continue
    #     # 
    #     '''left_idx = max(start - 4000, 0) // stride
    #     right_idx = (end + 4000 - win_size) // stride if (end + 4000) < reg_end else win_num - 1
    #     surround_dp_ls = []
    #     for i in range(left_idx, start_idx):
    #         surround_dp_ls.append(win_dp_ls[i])
    #     for i in range(end_idx, right_idx):
    #         surround_dp_ls.append(win_dp_ls[i])
    #     surround_dp = np.median(surround_dp_ls)
    #     clip_ratio = max(reg_clip_ls) // surround_dp if surround_dp != 0 else 0'''
    #     span_ls = []
    #     for read in bam_reader.fetch(ctg, start, end):
    #         if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < MIN_MAPPING_QUALITY or read.query_alignment_length < MIN_ALIGN_LENGTH or (read.query_alignment_length / read.query_length) < MIN_ALIGN_RATE:
    #             continue
    #         if read.reference_start < max(0, start - 500) and read.reference_end > min(ctg_len, end + 500):
    #             span_ls.append(read)
    #     #  
    #     if len(span_ls) < min_span_num:
    #         filtered_ls.append([reg[0], reg[1], reg[2], "low_supp_mis"])
    #     else:
    #         min_span_num = max(min_span_num, 0.2*len(span_ls))
    #         # cal indel of the span
    #         Inss_ls = []
    #         Dels_ls = []
    #         for i in range(len(span_ls)):   # 收集span_ls读数的indels
    #             read = span_ls[i]
    #             indels = cal_idels(read, start, end)    # max(0, start - 500), min(ctg_len, end + 500)
    #             Inss_ls.append(indels[0])
    #             Dels_ls.append(indels[1])
    #         # 
    #         Inss_ls.sort()
    #         Dels_ls.sort()
    #         # check indels
    #         avg_ins = sum(Inss_ls[:min_span_num]) // min_span_num
    #         avg_del = sum(Dels_ls[:min_span_num]) // min_span_num
    #         if avg_ins < ins_threshold and avg_del < del_threshold: # no mis
    #             out_ls.append([reg, "No_mis"])
    #         else:
    #             filtered_ls.append([reg[0], reg[1], reg[2], "Medium_mis"])
    
    ##
    print("filtered_ls:", filtered_ls)
    with open(filtered_bed, "w") as f:
        for ele in filtered_ls:
            f.write("{}\t{}\t{}\t{}\t{}\n".format(ele[0], ele[1], ele[2], ele[3], ele[4]))
    print("{}:{}-{} out_ls:{}".format(ctg, reg_start, reg_end, out_ls))
    print("Stop classify and filter2")

def find_candidate(ctg, reg_start, reg_end, ref, bam, out_dir, dp_info:Depth_info, config):
    
    bam_reader = pysam.AlignmentFile(bam, "rb", index_filename=bam + ".bai")
    ctg_len = bam_reader.get_reference_length(ctg)
    print("Find mis for:{}:{}-{}".format(ctg, reg_start, reg_end))
    # OUT: 
    # 
    reg_len = reg_end - reg_start
    win_size = config["win_size"]
    stride = win_size // 2
    win_num = (reg_len - win_size) // stride + 1 # 末端？？
    # dp params
    dp_params = config["dp_params"]
    dp_lower_bound = dp_params["dp_lower_bound"]
    dp_upper_bound = dp_params["dp_upper_bound"]
    lower_dp = dp_info.whole_dp * dp_lower_bound
    upper_dp = dp_info.whole_dp * dp_upper_bound
    dp_ls = dp_info.dp_ls
    dp_win_size = dp_info.win_size
    
    # print("block_batch:", block_batch, "block_size:", block_size, 'block_num:', block_num)
    # clip params
    clip_params = config["clip_params"]
    min_clip_portion = clip_params["min_clip_portion"]
    min_clip_num = min_clip_portion * dp_info.whole_dp
    min_clip_len = clip_params["min_clip_len"]
    # pileup params
    pileup_params = config["pileup_params"]
    win_size, step_size, min_correct_portion, max_differ_portion, max_disagree_portion, cluster_dis \
        = pileup_params["win_size"], pileup_params["step_size"], pileup_params["min_correct_portion"], \
        pileup_params["max_differ_portion"], pileup_params["max_disagree_portion"], pileup_params["cluster_dis"]
    
    ## -----------------------一、cal info-----------------------
    # 1、get pileup info
    min_MQ = pileup_params['min_MQ']
    print("Get pileup info")
    region = ctg + ":" + str(reg_start) + "-" + str(reg_end)
    pileup_dir = os.path.join(out_dir, "pileup")
    pileup_file = os.path.join(pileup_dir, region + ".pileup.txt")
    pileup_cmd = ["samtools mpileup -B", "-q", str(min_MQ), "-aa", "-d 200","-r", region, "-f", ref, bam, "-o", pileup_file]
    subprocess.check_call(" ".join(pileup_cmd), shell=True)
    # pileup_stream = pysam.mpileup("-B", "-q", "20", "-aa", "-d 100","-r", region, "-f", ref, bam)    # 计算pileup信息
    print("Get pileup done")
    pileup_dict={"contig":[],"correct":[],"ambiguous":[],"insert":[],"deletion":[],
                 "disagree":[],"depth":[],"differ":[]}
    # for line in pileup_stream.split("\n"):
    with open(pileup_file, "r") as f:
        for line in f:
            record = line.strip().split()
            # print(record)
            if not record:
                continue
            ## 
            pileup_dict['contig'].append(record[0])
            match_detail=record[4]
            pileup_dict['correct'].append(match_detail.count('.')+match_detail.count(','))
            pileup_dict["depth"].append(int(record[3]))     # 反映的是pileup计算的depth，与实际计算的可以不一致
            st = ''.join(re.split('[\+|\-][0-9]+[ATCGatcg]+', match_detail))
            st_counter = Counter(st)
            disagree_numd = st_counter['a'] + st_counter['A'] + st_counter['g'] + st_counter['G'] + st_counter['c'] + st_counter['C'] + st_counter['t'] + st_counter['T']   # 单碱基的不一致，snp
            pileup_dict["disagree"].append(disagree_numd)
            match_counter = Counter(match_detail)
            differ_numd = match_counter['a'] + match_counter['A'] + match_counter['g'] + match_counter['G'] + match_counter['c'] + match_counter['C'] + match_counter['t'] + match_counter['T'] - disagree_numd # 包括插入删除的大小信息
            # ins_numd = 
            pileup_dict["differ"].append(differ_numd)
    
    window_pileup_dict={"contig":[],"start_pos":[],"end_pos":[],"correct_portion":[],
                        "disagree_portion":[],"differ_portion":[]}
    for i in range(0, len(pileup_dict["correct"]), step_size):
        start = i
        end = i + win_size if i + win_size <= len(pileup_dict["correct"]) else len(pileup_dict["correct"])
        ## 
        win_start = start + reg_start
        win_end = win_start + win_size if win_start + win_size < reg_end else reg_end
        total = np.sum(pileup_dict["depth"][start:end])
        window_pileup_dict["contig"].append(ctg)
        window_pileup_dict["start_pos"].append(win_start)
        window_pileup_dict["end_pos"].append(win_end)
        if total == 0:
            window_pileup_dict["correct_portion"].append(np.nan)
            window_pileup_dict["differ_portion"].append(np.nan)
            window_pileup_dict["disagree_portion"].append(np.nan)
            '''ls.append([ctg, start + reg_start, end + reg_start])'''
            continue
        window_pileup_dict["correct_portion"].append(np.sum(pileup_dict['correct'][start:end])/total)
        window_pileup_dict["disagree_portion"].append(np.sum(pileup_dict["disagree"][start:end])/total)
        window_pileup_dict["differ_portion"].append(np.sum(pileup_dict["differ"][start:end])/total) 
        # print("{}-{}: correct:{},differ:{}".format(start + reg_start, end + reg_start, window_pileup_dict["correct_portion"][-1], window_pileup_dict["differ_portion"][-1]))
        '''if window_pileup_dict["correct_portion"][-1] < min_correct_portion \
            or window_pileup_dict["differ_portion"][-1] > max_differ_portion \
            or window_pileup_dict["disagree_portion"][-1] > max_disagree_portion:
            ls.append([ctg, start + reg_start, end + reg_start])'''
    # pileup_out = out_dir + "/" + ctg + ":" + str(reg_start) + "-" + str(reg_end) + "_pileup_feature.txt"
    # write_pileup_dic(window_pileup_dict, pileup_out)
    win_correct_portion = window_pileup_dict["correct_portion"]
    win_differ_portion = window_pileup_dict["differ_portion"]
    win_disagree_portion = window_pileup_dict["disagree_portion"]
    print("Get pileup info done")
    # 2、get clip info
    print("Get clip info")
    clip_ls = [0] * reg_len # 记录contig每个位置>min_clip_len的clip数目
    min_clip_len = 500      # 500 300
    # MIN_MAPPING_QUALITY = 20    # 10
    min_MQ = clip_params['min_MQ']
    for read in bam_reader.fetch(ctg, start=reg_start, end=reg_end):
        if read.is_unmapped or read.is_secondary or read.mapping_quality < min_MQ:
            continue
        # ref_pos = read.reference_start
        cigar = read.cigartuples

        # left
        left = cigar[0]
        if left[0] == 4 or left[0] == 5:
            if left[1] > min_clip_len:
                clip_ls[read.reference_start] += 1
        # right
        right = cigar[-1]
        if right[0] == 4 or right[0] == 5:
            if right[1] > min_clip_len:
                clip_ls[read.reference_end-1] += 1
    win_clip_num = [0] * win_num    # clip num of a window
    for i in range(win_num):    # i*stride, i*stride+win_size -> reg_start+i*stride, reg_start+i*stride+win_size
        win_clip_num[i] = cal_sum(clip_ls, reg_start+i*stride, reg_start+i*stride+win_size) if reg_start+i*stride+win_size < reg_len else cal_sum(clip_ls, reg_start+i*stride, reg_len)
    # 3、get dp info, skipped
    print("Get dp info")
    win_dp_ls = [0] * win_num
    norm_win_size = win_size // dp_win_size
    norm_stride = stride // dp_win_size
    bias = reg_start // dp_win_size
    for i in range(win_num):
        l = i * norm_stride + bias
        r = i * norm_stride + norm_win_size + bias if i * norm_stride + norm_win_size + bias <= len(dp_ls) else len(dp_ls)
        win_dp_ls[i] = cal_sum(dp_ls, l, r) / (r - l)
        # print("{}, Win {}:".format(ctg, i), win_dp_ls[i], dp_ls[l:r])
    
    ## -----------------------二、find regions-----------------------
    print("dp_lower_bound:{}, dp_upper_bound:{}, lower_dp:{}, upper_dp:{}".format(dp_lower_bound, dp_upper_bound, lower_dp, upper_dp))
    print("min_clip_num:", min_clip_num)
    print()
    ''' ctg:reg_start-reg_end'''
    '''
    win_clip_num
    win_correct_portion
    win_differ_portion
    win_disagree_portion
    win_dp_ls
    '''
    reg_id = ctg + ":" + str(reg_start) + "-" + str(reg_end)
    info_bed_dir = os.path.join(out_dir, "info")
    if not os.path.isdir(info_bed_dir):os.makedirs(info_bed_dir)
    info_bed = os.path.join(info_bed_dir, reg_id + ".info.bed")
    f1 = open(info_bed, "w")
    f1.write("#contig\tstart\tend\tdp\tclip_num\tcorrect_portion\tdiffer_portion\tdisagree_portion\n")
    candidate_dir = os.path.join(out_dir, "candidate")
    if not os.path.isdir(candidate_dir):os.makedirs(candidate_dir)
    candidate_bed = os.path.join(candidate_dir, reg_id + ".bed")
    f2 = open(candidate_bed, "w")
    candidate_ls = []
    # if config['']
    apply_dp = config['apply_dp']
    apply_clip = config['apply_clip']
    apply_pileup = config['apply_pileup']
    print("apply_dp:{}, apply_cip:{}, apply_pileup:{}".format(apply_dp, apply_clip, apply_pileup))
    for i in range(win_num):
        win_start = reg_start + i * stride
        win_end = reg_start + i * stride + win_size if reg_start + i * stride + win_size <= reg_end else reg_end
        # print(win_start, "-", win_end)
        f1.write("{}\t{}\t{}\t{:.2f}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\n".format(ctg, win_start, win_end, win_dp_ls[i], win_clip_num[i], win_correct_portion[i], win_differ_portion[i], win_disagree_portion[i]))
        # if win_dp_ls[i] == 0: 
        #     print(win_correct_portion[i] <= min_correct_portion)
        #     print(win_correct_portion[i])
        '''筛选条件：dp、clip、pileup
        消融实验：
        '''
        '''if win_dp_ls[i] <= lower_dp \
            or win_dp_ls[i] >= upper_dp \
            or win_clip_num[i] >= min_clip_num \
            or win_correct_portion[i] <= min_correct_portion \
            or win_differ_portion[i] >= max_differ_portion \
            or win_disagree_portion[i] >= max_disagree_portion:
            candidate_ls.append([ctg, win_start, win_end])
            f2.write("{}\t{}\t{}\t{:.2f}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\n".format(ctg, win_start, win_end, win_dp_ls[i], win_clip_num[i], win_correct_portion[i], win_differ_portion[i], win_disagree_portion[i]))
            # print(win_start, "-", win_end)'''

        check = [False, False, False]
        if apply_dp:
            if win_dp_ls[i] <= lower_dp \
                or win_dp_ls[i] >= upper_dp:
                check[0] = True
        if apply_clip:
            if win_clip_num[i] >= min_clip_num:
                check[1] = True
        if apply_pileup:
            if win_correct_portion[i] <= min_correct_portion \
                or win_differ_portion[i] >= max_differ_portion \
                or win_disagree_portion[i] >= max_disagree_portion:
                check[2] = True
        if check[0] or check[1] or check[2]:
            candidate_ls.append([ctg, win_start, win_end])
            f2.write("{}\t{}\t{}\t{:.2f}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{}\n".format(ctg, win_start, win_end, win_dp_ls[i], win_clip_num[i], win_correct_portion[i], win_differ_portion[i], win_disagree_portion[i], check))
        # if win_dp_ls[i] <= lower_dp:
        #     print(win_start, "-", win_end, "type1")
        # elif win_dp_ls[i] >= upper_dp:
        #     print(win_start, "-", win_end, "type2")
        # elif win_clip_num[i] >= min_clip_num:
        #     print(win_start, "-", win_end, "type3")
        # elif win_correct_portion[i] <= min_correct_portion:
        #     print(win_start, "-", win_end, "type4")
        # elif win_differ_portion[i] >= max_differ_portion:
        #     print(win_start, "-", win_end, "type5")
        # elif win_disagree_portion[i] >= max_disagree_portion:
        #     print(win_start, "-", win_end, "type6")
    
    ## cluster by dis
    print("len(candidate_ls):", len(candidate_ls))
    clu_ls = cluster_by_dis(candidate_ls, 5000)
    print("len(clu_ls)", len(clu_ls))
    
    ##
    print("clu_ls:", clu_ls)
    clu1_bed = os.path.join(candidate_dir, reg_id + ".clu1.bed")
    f3 = open(clu1_bed, "w")
    # boundary = 10000
    for reg in clu_ls:
        chr, start, end = reg
        start_idx = start // stride
        end_idx = (end - win_size) // stride if end < reg_end else win_num - 1
        # print(start_idx, end_idx)
        ls1 = [round(win_dp_ls[i], 4) for i in range(start_idx, end_idx + 1)]
        ls2 = [round(win_clip_num[i], 4) for i in range(start_idx, end_idx + 1)]
        # 
        ls3 = [round(win_correct_portion[i], 4) for i in range(start_idx, end_idx + 1)]
        ls4 = [round(win_differ_portion[i], 4) for i in range(start_idx, end_idx + 1)]
        ls5 = [round(win_disagree_portion[i], 4) for i in range(start_idx, end_idx + 1)]
        f3.write("{}\t{}\t{}\t{}\t{}\n".format(chr, start, end, ls1, ls2))

    ### -----------------------三、Filter-----------------------
    filter(bam, ctg, reg_start, reg_end, ctg_len, clu_ls, out_dir)
    filter2(bam, ctg, reg_start, reg_end, ctg_len, clu_ls, out_dir, config)
    
    f1.close()
    f2.close()
    f3.close()


def run_find_pipe(ref, bam, ctg_ls, out_dir, threads, config):
    print("----------------find mis pipe----------------")
    bam_reader = pysam.AlignmentFile(bam, "rb", index_filename=bam + ".bai")
    ## ----------------0、get dp info----------------
    print("----------------Start get_dp_info!!!----------------")
    dp_win_size = config["dp_params"]["dp_win_size"]
    block_size = config["dp_params"]["block_size"]
    # min_MQ = 10     # for hifi assembly
    min_MQ = config['dp_params']['min_MQ']
    dp_info_dic = get_dp_info_parallel(bam, threads, out_dir, dp_win_size, block_size, min_MQ)# get dpinfo from mosdepth，include avg_depth
    print("----------------get_dp_info done!!!----------------")
    
    ## ----------------1、find candidate----------------
    # dirs
    candidate_dir = os.path.join(out_dir, "candidate")
    if not os.path.isdir(candidate_dir):os.mkdir(candidate_dir)
    info_bed_dir = os.path.join(out_dir, "info")
    if not os.path.isdir(info_bed_dir):os.mkdir(info_bed_dir)
    filtered_dir = os.path.join(out_dir, "filtered")
    if not os.path.isdir(filtered_dir):os.mkdir(filtered_dir)
    filtered2_dir = os.path.join(out_dir, "filtered2")
    if not os.path.isdir(filtered2_dir):os.mkdir(filtered2_dir)
    pileup_dir = os.path.join(out_dir, "pileup")
    if not os.path.isdir(pileup_dir):os.mkdir(pileup_dir)
    # rm_cmd1 = ["rm", filtered_dir+"/*.bed"]
    # rm_cmd2 = ["rm", filtered2_dir+"/*.bed"]
    # subprocess.check_call(" ".join(rm_cmd1), shell=True)
    # subprocess.check_call(" ".join(rm_cmd2), shell=True)    
    # clear filtered dir
    clear_bed_files(filtered_dir)
    clear_bed_files(filtered2_dir)
    # run find
    pool = Pool(processes=threads)
    for ctg in ctg_ls:
        ctg_len = bam_reader.get_reference_length(ctg)
        reg_start = 0
        reg_end = ctg_len
        # find_candidate(ctg, 0, ctg_len, ref, bam, out_dir, dp_info_dic[ctg], config)
        pool.apply_async(find_candidate, args=(ctg, reg_start, reg_end, ref, bam, out_dir, dp_info_dic[ctg], config), callback=call_back, error_callback=error_call_back)
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    print("----------------------find candidate done----------------------")
    # 
    merge_dir = filtered_dir + "/merge"
    merge_dir2 = filtered2_dir + "/merge"
    merge_bed = merge_dir + "/merge.bed"
    merge_bed2 = merge_dir2 + "/merge.bed"
    make_dir(merge_dir)
    make_dir(merge_dir2)
    merge_cmd1 = ["cat", filtered_dir+"/*.bed", ">", merge_bed]
    merge_cmd2 = ["cat", filtered2_dir+"/*.bed", ">", merge_bed2]
    subprocess.check_call(" ".join(merge_cmd1), shell=True)
    subprocess.check_call(" ".join(merge_cmd2), shell=True)
    # post process
    '''merge_bed = candidate_dir+"/merge.clu1.bed"
    merge_cmd = ["cat", candidate_dir+"/*clu1.bed", "| cut -f1-3", ">", merge_bed]
    subprocess.check_call(" ".join(merge_cmd), shell=True)
    candidate_reg_dic = defaultdict(list)
    with open(merge_bed, "r") as f:
        for line in f:
            if line.startswith("#"): continue
            ctg, start, end = line.strip().split("\t")
            start, end = int(start), int(end)
            candidate_reg_dic[ctg].append([start, end])'''
    ## ----------------2、collect sigs from cigar parallel----------------
    '''## ----------------2、collect sigs from cigar parallel----------------
    # params
    min_length = 30
    min_mapq = 20
    max_split_parts = 7
    min_read_len = 500
    min_siglength = 10
    merge_del_threshold = 100     # 500   # sigs的合并距离
    merge_ins_threshold = 100   # 500
    MaxSize = 100000
    # run
    sigs_dir = out_dir + "/signatures"
    if not os.path.isdir(sigs_dir):os.mkdir(sigs_dir)
    pool = Pool(processes=threads)
    for ctg in ctg_ls:
        task = [ctg, 0, bam_reader.get_reference_length(ctg)]
        pool.apply_async(collect_signature.single_pipe, args=(bam, min_length, min_mapq, max_split_parts, min_read_len, out_dir + "/", 
                task, min_siglength, merge_del_threshold, merge_ins_threshold, MaxSize, candidate_reg_dic[ctg]), callback=call_back, error_callback=error_call_back)
        # collect_signature.single_pipe(bam, min_length, min_mapq, max_split_parts, min_read_len, sigs_dir, 
        #         task, min_siglength, merge_del_threshold, merge_ins_threshold, MaxSize, candidate_reg_dic[ctg])
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    
    ## ----------------3、cluster、filter----------------
    pool = Pool(processes=threads)
    for reg in candidate_reg_ls:
        pool.apply_async(cluster_sigs.clustar_all, args=(bam, reg, out_dir, dp_info_dic))
        pool.apply_async(cluster_sigs.clustar_inv, args=())
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    ## ----------------4、----------------'''
    return



if __name__ == "__main__":
    # print("Start")
    # ctg = 'NC_001133.9'
    # ctg = 'NC_001134.8'
    # ctg = 'NC_001136.10'
    # reg_start = 0
    # reg_end = 230218
    # ref = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/GCF_000146045.fna"
    # bam = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/aln2alt.sort.bam"
    # out_dir = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/my_pipe/step2_candidate_regions"
    
    # bam_reader = pysam.AlignmentFile(bam, "rb", index_filename=bam + ".bai")
    # chr_len = bam_reader.get_reference_length(ctg)
    # reg_end = chr_len
    # dp_win_size = 100
    # # dp_out_file = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/my_pipe/step2_candidate_regions/depths/NC_001133.9/NC_001133.9.regions.bed.gz"
    # dp_dir = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/my_pipe/step2_candidate_regions/depths"
    # dp_out_file = os.path.join(dp_dir, ctg, ctg + ".regions.bed.gz")
    # dp_ls = read_mosdepth_dp_file(dp_out_file)
    # block_size = 10000000
    # whole_dp = 37.75
    # # dp_info = Depth_info(ctg, chr_len, dp_ls, dp_win_size, block_size, whole_dp)
    # config_file = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/find_mis_configs/Config.yaml"
    # with open(config_file, "r") as f: # config参数获取
    #     config = yaml.safe_load(f.read())     # 获取部分参数
    # # find_candidate(ctg, reg_start, reg_end, ref, bam, out_dir, dp_info, config["step2"])
    
    # # Test
    # '''pool = Pool(processes=20)
    # bam_reader = pysam.AlignmentFile(bam, "rb", index_filename=bam + ".bai")
    # for ctg in bam_reader.references:
    #     chr_len = bam_reader.get_reference_length(ctg)
    #     dp_dir = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/my_pipe/step2_candidate_regions/depths"
    #     dp_out_file = os.path.join(dp_dir, ctg, ctg + ".regions.bed.gz")
    #     dp_ls = read_mosdepth_dp_file(dp_out_file)
    #     dp_win_size = 100
    #     block_size = 10000000
    #     whole_dp = 37.75
    #     dp_info = Depth_info(ctg, chr_len, dp_ls, dp_win_size, block_size, whole_dp)
    #     # find_candidate(ctg, 0, bam_reader.get_reference_length(ctg), ref, bam, out_dir, dp_info, config["step2"])
    #     pool.apply_async(find_candidate, args=(ctg, 0, bam_reader.get_reference_length(ctg), ref, bam, out_dir, dp_info, config["step2"]))
    #     # pass
    # pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    # pool.join() # 等待进程池中的所有进程执行完毕'''
    # # ref = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/chm13/simu/template_simu.fasta"
    # # bam = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/chm13/simu/aln2simu.sort.bam"
    # # out_dir = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/chm13/my_pipe/step2_candidate_regions"
    # # bam_reader = pysam.AlignmentFile(bam, "rb", index_filename=bam + ".bai")
    
    # for ctg in bam_reader.references:
    #     chr_len = bam_reader.get_reference_length(ctg)
    #     dp_dir = out_dir + "/depths"
    #     dp_out_file = os.path.join(dp_dir, ctg, ctg + ".regions.bed.gz")
    #     dp_ls = read_mosdepth_dp_file(dp_out_file)
    #     dp_win_size = 100
    #     block_size = 10000000
    #     whole_dp = 37.75
    #     dp_info = Depth_info(ctg, chr_len, dp_ls, dp_win_size, block_size, whole_dp)
    #     find_candidate(ctg, 0, bam_reader.get_reference_length(ctg), ref, bam, out_dir, dp_info, config["step2"])
    # exit(1)


    ##  ------------------------Test filter---------------------------
    # out_dir = "/public/home/hpc214712170/Test/mis_detect/asm/thaliana_hifi/my_pipe/test"
    # bam = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/mis_detect/asm/thaliana_hifi/data/aln2asm.sort.bam"
    # clu_dir = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/mis_detect/asm/thaliana_hifi/my_pipe/step2_candidate_regions/candidate"
    # 
    config_file = "/public/home/hpc214712170/Test/mis_detect/asm/chm13_ont/flye/Config.yaml"
    with open(config_file, "r") as f: # config参数获取
        config = yaml.safe_load(f.read())     # 获取部分参数
    bam = "/public/home/hpc214712170/Test/mis_detect/asm/chm13_ont/flye/data/aln2asm.sort.bam"
    filtered_dir = "/public/home/hpc214712170/Test/mis_detect/asm/chm13_ont/flye/my_pipe/step2_candidate_regions/filtered2"
    out_dir = "/public/home/hpc214712170/Test/mis_detect/asm/chm13_ont/flye/my_pipe/step2_candidate_regions"
    clu_dir = "/public/home/hpc214712170/Test/mis_detect/asm/chm13_ont/flye/my_pipe/step2_candidate_regions/candidate"
    bam_reader = pysam.AlignmentFile(bam, "rb", index_filename=bam + ".bai")
    file_ls = os.listdir(clu_dir)
    pool = Pool(processes=40)
    for file in file_ls:
        if file.endswith("clu1.bed"):
            file_path = clu_dir + "/" + file
            ctg = file.split(":")[0]
            reg_start = 0
            ctg_len = bam_reader.get_reference_length(ctg)
            reg_end = ctg_len
            clu_ls = []
            with open(file_path, "r") as f:
                for line in f:
                    if not line: continue
                    fields = line.strip().split("\t")
                    ctg, start, end = fields[0:3]
                    start, end = int(start), int(end)
                    clu_ls.append([ctg, start, end])
            # filter2(bam, ctg, reg_start, reg_end, ctg_len, clu_ls, out_dir, config['step2'])
            pool.apply_async(filter2, args=(bam, ctg, reg_start, reg_end, ctg_len, clu_ls, out_dir, config['step2']), callback=call_back, error_callback=error_call_back)
    # pool = Pool(processes=40)
    # for reg in candidate_reg_ls:
    #     pool.apply_async(filter, args=(), callback=call_back, error_callback=error_call_back)
    #     # collect_signature.single_pipe(bam, min_length, min_mapq, max_split_parts, min_read_len, sigs_dir, 
    #     #         task, min_siglength, merge_del_threshold, merge_ins_threshold, MaxSize, candidate_reg_dic[ctg])
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    merge_bed = filtered_dir + "/merge/merge.bed"
    merge_cmd = ["cat", filtered_dir+"/*.bed", ">", merge_bed]
    subprocess.check_call(" ".join(merge_cmd), shell=True)
    exit(1)
    ## win_check()
    # params = config['step2']["filter2"]
    # ctg = "contig_140"
    # start = 71084200
    # end = 71094200
    # ctg_len = 72820461
    # bam_reader = pysam.AlignmentFile(bam, "rb", index_filename=bam + ".bai")
    # res = win_check(ctg, start, end, ctg_len, bam_reader, params)
    # print(res)
    # exit(1)
    ## 
    candidate_dir = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/my_pipe/step2_candidate_regions/candidate"
    merge_bed = candidate_dir+"/merge.clu1.bed"
    merge_cmd = ["cat", candidate_dir+"/*clu1.bed", "| cut -f1-3", ">", merge_bed]
    subprocess.check_call(" ".join(merge_cmd), shell=True)
    candidate_reg_dic = defaultdict(list)
    candidate_reg_ls = []
    with open(merge_bed, "r") as f:
        for line in f:
            if line.startswith("#"): continue
            ctg, start, end = line.strip().split("\t")
            start, end = int(start), int(end)
            candidate_reg_dic[ctg].append([start, end])
            candidate_reg_ls.append([ctg, start, end])
    min_length = 30
    min_mapq = 20
    max_split_parts = 7
    min_read_len = 500
    min_siglength = 10
    merge_del_threshold = 100     # 500   # sigs的合并距离
    merge_ins_threshold = 100   # 500
    MaxSize = 100000
    sigs_dir = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/my_pipe/step2_candidate_regions/signatures"
    out_dir = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/my_pipe/step2_candidate_regions/"
    if not os.path.isdir(sigs_dir):os.mkdir(sigs_dir)
    pool = Pool(processes=20)
    for reg in candidate_reg_ls:
        task = reg
        pool.apply_async(collect_signature.single_pipe, args=(bam, min_length, min_mapq, max_split_parts, min_read_len, out_dir, 
                task, min_siglength, merge_del_threshold, merge_ins_threshold, MaxSize, None), callback=call_back, error_callback=error_call_back)
        # collect_signature.single_pipe(bam, min_length, min_mapq, max_split_parts, min_read_len, sigs_dir, 
        #         task, min_siglength, merge_del_threshold, merge_ins_threshold, MaxSize, candidate_reg_dic[ctg])
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    # for reg in candidate_reg_ls:
    #     collect_signature.single_pipe(bam, min_length, min_mapq, max_split_parts, min_read_len, out_dir, 
    #             reg, min_siglength, merge_del_threshold, merge_ins_threshold, MaxSize, None)
    # temporary_dir = out_dir
    # analysis_pools = Pool(processes=int(20))
    # cmd_del = ("cat %ssignatures/*.bed | grep -w DEL | sort -u -T %s | sort -k 2,2 -k 3,3n -T %s > %sDEL.sigs"%(temporary_dir, temporary_dir, temporary_dir, temporary_dir))
    # cmd_ins = ("cat %ssignatures/*.bed | grep -w INS | sort -u -T %s | sort -k 2,2 -k 3,3n -T %s > %sINS.sigs"%(temporary_dir, temporary_dir, temporary_dir, temporary_dir))
    # cmd_inv = ("cat %ssignatures/*.bed | grep -w INV | sort -u -T %s | sort -k 2,2 -k 3,3 -k 4,4n -T %s > %sINV.sigs"%(temporary_dir, temporary_dir, temporary_dir, temporary_dir))
    # cmd_tra = ("cat %ssignatures/*.bed | grep -w TRA | sort -u -T %s | sort -k 2,2 -k 5,5 -k 3,3 -k 4,4n -T %s > %sTRA.sigs"%(temporary_dir, temporary_dir, temporary_dir, temporary_dir))
    # cmd_dup = ("cat %ssignatures/*.bed | grep -w DUP | sort -u -T %s | sort -k 1,1r -k 2,2 -k 3,4n -T %s > %sDUP.sigs"%(temporary_dir, temporary_dir, temporary_dir, temporary_dir))
    # cmd_reads = ("cat %ssignatures/*.reads > %sreads.sigs"%(temporary_dir, temporary_dir))
    # for i in [cmd_ins, cmd_del, cmd_dup, cmd_tra, cmd_inv, cmd_reads]:
    #     analysis_pools.map_async(exe, (i,))
    # analysis_pools.close()
    # analysis_pools.join()
    # run_pipe(ref, bam, bam_reader.references, out_dir, 20, config["step2"])
    pass

'''def find_candidate_parallel(ref, bam_in, ctg_ls, out_dir, threads, config, args):
    print("get_dp_info!!!")
    dp_win_size = config["dp_params"]["dp_win_size"]
    block_size = config["dp_params"]["block_size"]
    ## get dpinfo from mosdepth，include avg_depth
    dp_info_dic = get_dp_info_parallel(bam_in, threads, out_dir, dp_win_size, block_size)
    print("get_dp_info done!!!")

    ## 
    dirs = [os.path.join(out_dir, "clip"), os.path.join(out_dir, "cov")]
    for dir in dirs:
        if not os.path.exists(dir): os.makedirs(dir)
    pool = Pool(processes=threads)
    for ctg in ctg_ls:
        pool.apply_async(run_find_candidate2, args=(bam_in, ctg, out_dir, config, dp_info_dic), callback=call_back, error_callback=error_call_back) 
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    # for ctg in ctg_ls:
    #     run_find_candidate2(bam_in, ctg, out_dir, config, dp_info_dic)

    ## find from pileup
    if config["apply_pileup"]:
        print("apply pileup")
        pileup_dir = os.path.join(out_dir, "pileup")
        if not os.path.exists(pileup_dir):os.makedirs(pileup_dir)
        # parse_pileup_parallel(ctg_ls, ref, bam_in, threads, config["pileup_params"], pileup_dir)
    
    print("find_candidate done!!!")


    ## bed merge
    bed_ls = []
    if config["apply_pileup"]:
        bed_ls.append(os.path.join(pileup_dir, "candidate_pileup.bed"))
    for ctg in ctg_ls:
        bed_ls.append(os.path.join(out_dir, "clip", ctg+"_clip.bed"))
        bed_ls.append(os.path.join(out_dir, "cov", ctg+"_cov.bed"))
    bed_out_merge = os.path.join(out_dir, "candidate.bed")
    # bed_merge(bed_ls, bed_out_merge)
    merge_cmd = ["cat", " ".join(bed_ls), "| sort -k 1,1 -k 2n,2", ">", bed_out_merge]
    subprocess.check_call(" ".join(merge_cmd), shell=True)'''
# /public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/mis_detect/asm/chm13/my_pipe2/step2_candidate_regions/pileup/contig_598:0-110030683.pileup.txt