from multiprocessing import Pool
import os
import math
import pysam
from create_consensus_by_bed import fasta_parser
from collections import defaultdict
'''Ne, Na = 0, 0
    ctg_len = bam_reader.get_reference_length(ctg)
    reg_start = 0
    reg_end = ctg_len
    region = ctg + ":" + str(reg_start) + "-" + str(reg_end)
    pileup_file = os.path.join(work_dir, "pileup", region + ".pileup.txt")
    pileup_file = os.path.join(work_dir, "pileup", ctg + ":" + str(reg_start) + "-" + str(reg_end - 1) + ".pileup.txt")  
    mis_bed = os.path.join(work_dir, "filtered", region + ".bed")
    mis_ls = []
    with open(mis_bed, "r") as f:
        for line in f:
            fields = line.strip().split()
            mis_ls.append([fields[0], int(fields[1]), int(fields[2])])        
    if len(mis_ls) > 0:
        mis = mis_ls[0]
    else:            
        mis = [ctg, -1, -1]
    mis_length = 0
    mis_length += mis[2] - mis[1]
    mis_idx = 0
    print(mis)
    with open(pileup_file, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            pos = int(fields[1])
            if pos >= mis[1] and pos <= mis[2]: 
                # print(pos, end="\t")
                continue # in mis region
            if pos > mis[2]:  # over this mis
                if mis_idx + 1 < len(mis_ls):
                    mis_idx += 1
                    mis = mis_ls[mis_idx]
                    mis_length += mis[2] - mis[1]
                    print(mis)
            # print(pos, end="\t")
            depth = int(fields[3])  # 不计算del为深度，不合理
            match_detail = fields[4]
            quality_detail = fields[5]
            # 
            read_nums = len(quality_detail)
            if read_nums < min_dp: continue
            Na += 1
            if match_detail.count('.') + match_detail.count(',') / read_nums < mm_rate:
                Ne += 1
                        
            if pos < mis[1]:
                depth = int(fields[3])  # 不计算del为深度，不合理
                match_detail = fields[4]
                quality_detail = fields[5]
            elif pos >= mis[1] and pos < mis[2]:
                continue
            elif pos > mis[2]:  # 超过当期mis
                if mis_idx + 1 < len(mis_ls):
                    mis_idx += 1
                    mis = mis_ls[mis_idx]'''
        # print(Ne, Na)

# def call_back(res):
#     print(res)

def error_call_back(error_code):
    print(error_code)
def simple_stats(ref, min_contig, work_dir):
    print("Simple stats")
    simple_file = os.path.join(work_dir, "simple.ststs")
    ref_dic = fasta_parser.read_sequence_dict(ref)
    f = open(simple_file, "w")
    ## stats
    length_ls = []
    ctg_ls = [] # 有效的ctg，执行了Mis检测的
    for ctg_id, ctg_seq in ref_dic.items():
        length_ls.append(len(ctg_seq))
        if len(ctg_seq) > min_contig:
            ctg_ls.append(ctg_id)
    # 
    length_ls = sorted(length_ls, reverse=True)
    total_length = sum(length_ls)
    f.write("Contigs: {}\n".format(len(length_ls)))
    f.write("Length: {}\n".format(total_length))
    f.write("Length > {}bp: {}\n".format(min_contig, len(ctg_ls)))
    f.write("Longest contig: {}\n".format(max(length_ls)))
    N50_target = total_length // 2
    tmp = 0
    for length in length_ls:
        tmp += length
        if tmp >= N50_target:
            f.write("N50: {}\n".format(length))
            break
    f.close()
    return ctg_ls
def cal_qv(Ne, Na):
    qv = -10 * math.log10(Ne / Na)
    return qv
def cal_Ne(work_dir, bam, ctg, min_dp, mm_rate):
    bam_reader = pysam.AlignmentFile(bam, "rb", index_filename=bam + ".bai")
    '''该函数通过解析pileup文件，计算qv计算所需要的Ne和Na'''
    Ne, Na = 0, 0
    ctg_len = bam_reader.get_reference_length(ctg)
    reg_start = 0
    reg_end = ctg_len
    region = ctg + ":" + str(reg_start) + "-" + str(reg_end)
    pileup_file = os.path.join(work_dir, "pileup", region + ".pileup.txt")
    pileup_file = os.path.join(work_dir, "pileup", ctg + ":" + str(reg_start) + "-" + str(reg_end) + ".pileup.txt")  
    mis_bed = os.path.join(work_dir, "filtered2", region + ".bed")
    mis_ls = []
    with open(mis_bed, "r") as f:
        for line in f:
            fields = line.strip().split()
            mis_ls.append([fields[0], int(fields[1]), int(fields[2])])        
    if len(mis_ls) > 0:
        mis = mis_ls[0]
    else:            
        mis = [ctg, -1, -1]
    mis_length = 0
    mis_length += mis[2] - mis[1]
    mis_idx = 0
    # print(mis)
    with open(pileup_file, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            pos = int(fields[1])
            if pos >= mis[1] and pos <= mis[2]: 
                # print(pos, end="\t")
                continue # in mis region
            if pos > mis[2]:  # over this mis
                if mis_idx + 1 < len(mis_ls):
                    mis_idx += 1
                    mis = mis_ls[mis_idx]
                    mis_length += mis[2] - mis[1]
                    # print(mis)
            # print(pos, end="\t")
            # depth = int(fields[3])  # 不计算del为深度，不合理
            match_detail = fields[4]
            quality_detail = fields[5]
            # 
            read_nums = len(quality_detail)
            if read_nums < min_dp: continue
            Na += 1
            if match_detail.count('.') + match_detail.count(',') / read_nums < mm_rate:
                Ne += 1  
            '''if pos < mis[1]:
                depth = int(fields[3])  # 不计算del为深度，不合理
                match_detail = fields[4]
                quality_detail = fields[5]
            elif pos >= mis[1] and pos < mis[2]:
                continue
            elif pos > mis[2]:  # 超过当期mis
                if mis_idx + 1 < len(mis_ls):
                    mis_idx += 1
                    mis = mis_ls[mis_idx]'''
    return Ne, Na
def get_qv(ref, bam, work_dir, ctg_ls, threads):
    print("Start cal QV")
    simple_file = os.path.join(work_dir, "simple.ststs")
    f = open(simple_file, "a+")
    # bam_reader = pysam.AlignmentFile(bam, "rb", index_filename=bam + ".bai")
    min_dp = 5
    mm_rate = 0.7   # 根据读数
    # 
    results = []
    pool = Pool(processes=threads)
    for ctg in ctg_ls:
        results.append(pool.apply_async(cal_Ne, args=(work_dir, bam, ctg, min_dp, mm_rate), error_callback=error_call_back))
    # results = [pool.apply_async(cal_Ne, args=(work_dir, bam_reader, ctg, min_dp, mm_rate)), callback=call_back, error_callback=error_call_back]
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    Ne = 0
    Na = 0
    for res in results:
        ctg_Ne, ctg_Na = res.get()
        Ne += ctg_Ne
        Na += ctg_Na
    qv = cal_qv(Ne, Na)
    f.write("QV: {}\n".format(qv))
    f.close()


if __name__ == "__main__":
    print("Start")
    min_ctg = 30000
    # ref = "/public/home/hpc214712170/Test/mis_detect/asm/sativa_hifi/flye/data/flye_hifi.fa"
    # work_dir = "/public/home/hpc214712170/Test/mis_detect/asm/sativa_hifi/flye/my_pipe/step2_candidate_regions"
    # bam = "/public/home/hpc214712170/Test/mis_detect/asm/sativa_hifi/flye/data/aln2asm.sort.bam"
    ref = "/public/home/hpc214712170/Test/mis_detect/asm/sativa_hifi/hifiasm/data/hifiasm_hifi.fa"
    work_dir = "/public/home/hpc214712170/Test/mis_detect/asm/sativa_hifi/hifiasm/my_pipe/step2_candidate_regions"
    bam = "/public/home/hpc214712170/Test/mis_detect/asm/sativa_hifi/hifiasm/data/aln2asm.sort.bam"
    ctg_ls = simple_stats(ref, min_ctg, work_dir)
    # ctg_ls = ["contig_424"]
    fa_in = ""
    threads = 40
    # exit(0)
    get_qv(fa_in, bam, work_dir, ctg_ls, threads)   # ref, bam, work_dir, ctg_ls, threads
    # print(cal_qv(10, 1000))
    pass
