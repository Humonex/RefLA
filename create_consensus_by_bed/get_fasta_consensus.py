# from Bio import SeqIO
# import pandas as pd
import pysam
import sys
import os
import random
import signal
import multiprocessing as mp
from collections import namedtuple
import argparse
'''
bed format1(candidate regions found by features)     info: low_dep,
# chr start   end   info
chr1    0   1000    low_dep
chr1    10001   2000    clip_reg

bed format2:
# chr start   end   info        operation
chr1    0   10    .             DEL     
chr1    12  13    AGGGGGGG      INS   
chr1    15  20    NNNNNNNNNN    N_fill
chr1    33  44    .             INV   
chr1    46  47    GGGAAATT      DUP   

region format :[ctg, start, end]
'''
Record = namedtuple('Record', ["chr_id", "start", "end", "info", "operation"])

def fasta_write(seq_id, seq, out_stream):
    out_stream.write(">"+seq_id+"\n")
    chars_per_line = 80
    for i in range(0, len(seq), chars_per_line):
        out_stream.write(f'{seq[i:i + chars_per_line]}\n')

def can_solve(vcf_reader:pysam.AlignmentFile, region)->bool:   # 判断低覆盖区域是否为del，若能检测出是del，则直接return True根据vcf来处理
    ctg, start, end = region
    cover_len = 0
    for vcf in vcf_reader.fetch(ctg, start, end):
        if vcf.info["SVTYPE"] == "DEL":
            cover_len += vcf.rlen   # rlen:record length on chrom/contig (aka rec.stop - rec.start)
    if cover_len > 0.7:
        return True
    return False

def get_record_in_vcf(vcf_reader:pysam.AlignmentFile, region):  # [ctg, start, end]
    record_ls = []    
    ctg, start, end = region
    for vcf in vcf_reader.fetch(ctg, start, end):
        record_ls.append([vcf.contig, vcf.start, vcf.stop, vcf.alts[0], vcf.info["SVTYPE"]])
    return record_ls    # record_ls format: [ctg, start, end, info, operation]  info use alt

def has_overlap(reg1, reg2):    # [1, r]  判断两个区间是否有交叠
    l1, r1 = reg1
    l2, r2 = reg2
    if r1 < l2 or r2 < l1:
        return False
    return True # has overlap

def get_record_size(record):    # record format:[chr_id, start, end, info, operation]
    chr_id, start, end, info, operation = record
    if operation == "INS":
        return len(info) - 1
    if operation == "DEL" or operation == "DUP" or operation == "INV":
        return int(end) - (int(start) + 1)
    if operation == "BND":
        return 1000
    return 0

def filter_record(record_ls, exclude_region_ls):    # region format :[ctg, start, end]
    if len(record_ls) < 1:
        return []
    
    ## filter by exclude region
    tmp_record_ls = []
    for record in record_ls:
        record_start, record_end = record[1], record[2]
        flag = 0
        for reg in exclude_region_ls:
            reg_start, reg_end = reg[1], reg[2]
            if record_start > reg_end:
                break
            if has_overlap([record_start, record_end], [reg_start, reg_end]):
                flag = 1
        if not flag:
            tmp_record_ls.append(record)
    record_ls = tmp_record_ls
    
    ## filter by record size
    # if len(record_ls) < 1:
    #     return []
    # tmp_record_ls = []


    ## filter by record, remove overlap record, remove the second one   按变异大小好像更好一点（filter小的）
    if len(record_ls) < 1:
        return []
    filtered_record_ls = []
    pre_record = record_ls[0]
    for i in range(1, len(record_ls)):
        record = record_ls[i]
        pre_reg = [pre_record[1], pre_record[2]]
        now_reg = [record[1], record[2]]
        # if not has_overlap(pre_reg, now_reg):     # remove the second one
        #     filtered_record_ls.append(pre_record)
        #     pre_record = record
        if not has_overlap(pre_reg, now_reg):
            filtered_record_ls.append(pre_record)
            pre_record = record
        else:
            if get_record_size(pre_record) < get_record_size(record):   # remove the small one
                pre_record = record
    filtered_record_ls.append(pre_record)

    filtered_record_ls = list(set([tuple(record) for record in filtered_record_ls]))
    return filtered_record_ls

def bad_region_cluster(bad_region_ls:list): # format: [chr_id, pre_start, pre_end]
    bad_region_ls.sort(key=lambda region:region[1])
    new_ls = []
    max_len = 3000
    pre_start = -1
    pre_end = -1
    for region in bad_region_ls:
        chr_id, start, end = region
        if start - pre_end < max_len:   # bam region可能存在连接
            pre_end = end
        else:
            if pre_start > -1:  # 非首次，要将其加入到ls中
                new_ls.append([chr_id, pre_start, pre_end])
            pre_start = start
            pre_end = end
    if pre_start > -1:
        new_ls.append([chr_id, pre_start, pre_end])
    return new_ls

def low_dep_cluster(low_dep_region:list):
    low_dep_region.sort(key=lambda region:region[1])
    new_ls = []
    max_len = 500
    pre_start = -1
    pre_end = -1
    for region in low_dep_region:
        chr_id, start, end = region
        if start - pre_end < max_len:   # bam region可能存在连接
            pre_end = end
        else:
            if pre_start > -1:  # 非首次，要将其加入到ls中
                new_ls.append([chr_id, pre_start, pre_end])
            pre_start = start
            pre_end = end
    if pre_start > -1:
        new_ls.append([chr_id, pre_start, pre_end])
    return new_ls

def prase_record_file(bed_in):
    record_ls = []
    with open(bed_in, "r") as bed_reader:     # 根据bed修改，
        for line in bed_reader:
            if not line:
                continue
            chr_id, start, end, info, operation = line.strip().split("\t")[:5]
            record_ls.append([chr_id, int(start), int(end), operation, info])
    return record_ls

def get_operation_record(ctg_ls, bed_in, vcf_in, bed_out, MIN_N_FILL):    # bed_in + vcf_in -> bed_out     record format:[chr_id, start, end, info, operation]
    print("Get_operation_record!!!")
    bed_in_region_ls = []
    low_dep_ls = []
    clip_reg_ls = []
    with open(bed_in, "r") as bed_reader:
        for line in bed_reader:
            if line.startswith("#"):
                continue
            chr_id, start, end, info = line.strip().split()
            start, end = int(start), int(end)
            if info == "low_dep":
                # low_dep_ls.append([chr_id, start, end, info])
                low_dep_ls.append([chr_id, start, end])
            elif info == "clip_reg":
                # clip_reg_ls.append([chr_id, start, end, info])
                clip_reg_ls.append([chr_id, start, end])
            bed_in_region_ls.append([chr_id, start, end, info])
    bed_in_region_ls.sort(key=lambda reg:(reg[0], reg[1]))    # sort by chr_id and start
    low_dep_ls.sort(key=lambda reg:(reg[0], reg[1]))
    clip_reg_ls.sort(key=lambda reg:(reg[0], reg[1]))

    record_ls_merge = []
    fout = open(bed_out, "w")   # path/consensus.bed
    vcf_reader = pysam.VariantFile(vcf_in, "r")
    for ctg in ctg_ls:
        ## get region of ctg
        ctg_low_dep_ls = []
        ctg_clip_reg_ls = []
        for reg in low_dep_ls:
            if reg[0] == ctg:
                ctg_low_dep_ls.append(reg)
        for reg in clip_reg_ls:
            if reg[0] == ctg:
                ctg_clip_reg_ls.append(reg)

        ## select regions
        ctg_bad_region_ls = []
        ctg_record_ls = []  # 记录需要进行的修改

        ### select from lowdep region
        ctg_low_dep_ls = low_dep_cluster(ctg_low_dep_ls)
        for reg in ctg_low_dep_ls:  # 
            chr_id, start, end = reg
            # if chr_id != ctg:
            #     continue
            if end - start > 20000:     # N_fill   20000
                ctg_bad_region_ls.append([chr_id, start, end])
            else:
                if can_solve(vcf_reader, [chr_id, start, end]):
                    res_ls = get_record_in_vcf(vcf_reader, reg)
                    ctg_record_ls.extend(res_ls)
                else:   # 
                    ctg_bad_region_ls.append([chr_id, start, end])
        ctg_bad_region_ls = bad_region_cluster(ctg_bad_region_ls)
        tmp_ls = []
        for region in ctg_bad_region_ls:    # filter small low depth region
            if region[2] - region[1] > MIN_N_FILL:
                tmp_ls.append(region)
        ctg_bad_region_ls = tmp_ls


        ## select from clip region
        for reg in ctg_clip_reg_ls:
            # chr_id, start, end = reg
            res_ls = get_record_in_vcf(vcf_reader, reg)
            ctg_record_ls.extend(res_ls)

        ## filter ctg_record_ls by bad region and remove duplicate operation
        ctg_record_ls = list(set([tuple(record) for record in ctg_record_ls]))
        ctg_record_ls.sort(key=lambda record:int(record[1]))
        ctg_record_ls = filter_record(ctg_record_ls, ctg_bad_region_ls)

        ## add N_fill region to ctg_record_ls
        for region in ctg_bad_region_ls:
            ctg_record_ls.append([region[0], region[1], region[2], ".", "N_fill"])
        ctg_record_ls.sort(key=lambda record:int(record[1]))

        ## write record to bed_out
        for record in ctg_record_ls:
            chr_id, start, end, info, operation = record[:5]
            fout.write("{}\t{}\t{}\t{}\t{}\n".format(chr_id, start, end, info, operation))
        record_ls_merge.extend(ctg_record_ls)
    fout.close()        
    vcf_reader.close()
    print("get_operation_record Done!!!")
    return record_ls_merge

def solve_region_by_record(ctg_ls, fasta_in, record_ls, fasta_out, MIN_CONTIG):     # 
    '''
    根据操作信息对参考序列进行修改
    '''
    fout = open(fasta_out, "w")
    for ctg in ctg_ls:
        # record_ls = []
        # with open(bed_in, "r") as bed_reader:     # 根据bed修改，
        #     for line in bed_reader:
        #         chr_id, start, end, info, operation = line.strip().split("\t")[:5]
        #         if chr_id == ctg:
        #             record_ls.append([chr_id, int(start), int(end), operation, info])

        ref_seq = pysam.FastaFile(fasta_in).fetch(ctg)
        new_ref_seq = ref_seq
        offect = 0  # 记录偏移量
        if len(ref_seq) < MIN_CONTIG:
            fasta_write(ctg, new_ref_seq, fout)
            continue

        ctg_record_ls = []
        for record in record_ls:
            if record[0] == ctg:
                ctg_record_ls.append(record)
        ctg_record_ls.sort(key=lambda record:int(record[1]))
        for record in ctg_record_ls:
            chr_id, start, end, info, operation = record
            start += offect
            end += offect
            if operation == "INS":
                new_ref_seq = new_ref_seq[:start] + info + new_ref_seq[end:]
                offect += len(info) - 1     # info(alt)记录了变异，但第一位是参考上的
            elif operation == "DEL":
                new_ref_seq = new_ref_seq[:start+1] + new_ref_seq[end:] # 包括start处
                offect -= end - (start + 1)          # del的INFO记录只有一位字符，通过start end求SV长度
            elif operation == "INV":    # alt:('<INV>',)
                new_ref_seq = new_ref_seq[:start+1] + new_ref_seq[start+1:end][::-1] + new_ref_seq[end:]    # reverse[start+1, end]
            elif operation == "DUP":    # DUP的alt记录为('<DUP>',)
                dup_length = end - (start + 1)
                new_ref_seq = new_ref_seq[:start] + ref_seq[(record[1]+1):(record[1]+1)+dup_length] + new_ref_seq[start:]   # 获取原序列上对应序列。在start处插入一个新序列
                offect += dup_length
            elif operation == "N_fill": # 
                N_FILL = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
                new_ref_seq = new_ref_seq = new_ref_seq[:start] + N_FILL + new_ref_seq[end:]    # 填充N，删除了[start, end)的元素
                offect = offect + len(N_FILL) - (end-start)
            else:
                continue
        fasta_write(ctg, new_ref_seq, fout)
    fout.close()

def run_consensus(candidate_bed, candidate_vcf, ctg_ls, out_dir, fasta_in):
    # for ctg in ctg_ls:
    bed_out = os.path.join(out_dir, "consensus.bed")
    fasta_out = os.path.join(out_dir, "consensus.fasta")
    record_ls = get_operation_record(ctg_ls, candidate_bed, candidate_vcf, bed_out)
    solve_region_by_record(ctg_ls, fasta_in, record_ls, fasta_out)


def main():
    # import time
    # start_t = time.time()
    # candidate_bed = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/find_candidate_regions/test_find/find_candidate.bed"
    # candidate_vcf = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/HG002/SV/cuteSV/NC30_hap1/NC_060930.1_hap1.vcf.gz"
    # bed_out = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/create_consensus_by_bed/test_consensus/operation.bed"
    # contig_ls = ["NC_060930.1"]
    # record_ls = get_operation_record(contig_ls, candidate_bed, candidate_vcf, bed_out)

    # fasta_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/T2T_CHM13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
    # fasta_out = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/create_consensus_by_bed/test_consensus/NC_060930.1_hap1.fasta"
    # solve_region_by_record(contig_ls, fasta_in, record_ls, fasta_out)

    # end_t = time.time()
    # print('cost %.3f second' % (end_t - start_t))   # 11s

    ## single pipe
    # candidate_bed="/public/home/hpc214712170/shixf/projects/ref-guided-assembly/test_flye_hapdup/NC_060930.1/my_pipe/hap1/candidate_regions/candidate.bed"
    # candidate_vcf="/public/home/hpc214712170/shixf/projects/ref-guided-assembly/test_flye_hapdup/NC_060930.1/my_pipe/hap1/cute_SV/cuteSV.sorted.vcf.gz"
    # ctg_ls="contig_1,contig_11,contig_12,contig_14,contig_15,contig_16,contig_17,contig_19,contig_2,contig_20,contig_21,contig_23,contig_27,contig_28,contig_29,contig_3,contig_32,contig_35,contig_41,contig_43,contig_45,contig_52,contig_53,contig_54,contig_55,contig_58,contig_60,contig_63,contig_64,contig_66,contig_7,contig_8,contig_9".split(",")
    # fasta_in="/public/home/hpc214712170/shixf/projects/ref-guided-assembly/test_flye_hapdup/NC_060930.1/assembly.fasta"
    # out_dir="/public/home/hpc214712170/shixf/projects/ref-guided-assembly/test_flye_hapdup/NC_060930.1/my_pipe/hap1/SV_consensus"
    # candidate_bed=sys.argv[1]
    # candidate_vcf=sys.argv[2]
    # ctg_ls=sys.argv[3].split(",")
    # fasta_in=sys.argv[4]
    # out_dir=sys.argv[5]

    
    # bed_out = os.path.join(out_dir, "consensus.bed")
    # fasta_out = os.path.join(out_dir, "consensus.fasta")
    # record_ls = get_operation_record(ctg_ls, candidate_bed, candidate_vcf, bed_out)
    # print("SV consensus on: {}".format(ctg_ls))
    # solve_region_by_record(ctg_ls, fasta_in, record_ls, fasta_out)

    ### add parser
    parser = argparse.ArgumentParser(description="get_SV_consensus")
    parser.add_argument("--fasta", dest="fasta_in", required=True)
    parser.add_argument("--vcf", dest="candidate_vcf", required=True)
    parser.add_argument("--out-dir", dest="out_dir", required=True)
    parser.add_argument("--contig_ls", dest="ctg_ls", required=True, help="format like:chr1,chr2")
    parser.add_argument("--bed", dest="candidate_bed", required=True, help="provide candidate bed, will get record from this")
    parser.add_argument("--record_bed", dest="record_ls_bed", default=False, help="provide operation record file, if provide it, will skip get record step and follow this record ls")
    parser.add_argument("--min-contig", dest="min_contig", help="skip contig shorter than this, keep with raw", default=2000000, type=int)   # 调 1,000,000     5,000,000
    parser.add_argument("--min-fill", dest="min_N_fill", default=5000, help="min region to perform N_fill", type=int)

    args = parser.parse_args()
    bed_out = os.path.join(args.out_dir, "consensus.bed")
    fasta_out = os.path.join(args.out_dir, "consensus.fasta")
    
    fasta_reader = pysam.FastaFile(args.fasta_in)
    process_ctg_ls = []
    for ctg in args.ctg_ls.split(","):
        if fasta_reader.get_reference_length(ctg) > int(args.min_contig):
            process_ctg_ls.append(ctg)
    ##
    if args.record_ls_bed:  # 提供了record_ls，直接使用
        print("Use record ls by provide!!")
        record_ls = prase_record_file(args.record_ls_bed)
    else:
        record_ls = get_operation_record(process_ctg_ls, args.candidate_bed, args.candidate_vcf, bed_out, int(args.min_N_fill))
    ##
    print("SV consensus on: {} !!!".format(process_ctg_ls))
    solve_region_by_record(args.ctg_ls.split(","), args.fasta_in, record_ls, fasta_out, int(args.min_contig))    # 提供初始contig，因为其他contig也要输出
    
    
if __name__ == "__main__":
    main()


'''
def run_consensus_parallel(candidate_bed, candidate_vcf, ctg_ls, out_dir, fasta_in, num_threads):
    all_reference_ids = ctg_ls
    random.shuffle(all_reference_ids)
    chunk_size = len(all_reference_ids) // num_threads + 1
    threads = []
    orig_sigint = signal.signal(signal.SIGINT, signal.SIG_IGN)
    for i in range(num_threads):
        contigs_list = all_reference_ids[i*chunk_size:(i+1)*chunk_size]
        if not contigs_list:
            continue
        # threads.append(mp.Process(target= , args=()))
    signal.signal(signal.SIGINT, orig_sigint)
    for t in threads:
        t.start()
    try:
        for t in threads:
            t.join()
            if t.exitcode != 0:
                raise Exception("One of the processes exited with code: {0}".format(t.exitcode))
    except KeyboardInterrupt:
        for t in threads:
            t.terminate()
        raise
'''
