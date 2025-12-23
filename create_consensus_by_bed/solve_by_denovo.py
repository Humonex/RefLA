import logging
import subprocess
import pysam
import os
import shutil
from Bio import Seq
import time
from multiprocessing import Pool
from collections import namedtuple, defaultdict
from create_consensus_by_bed.Utils import Run_for_denovo, Record, write_reg_to_bed, make_dir, run_cmd_ls, DepthRec, Asm_Region
from create_consensus_by_bed import fasta_parser

'''
特殊情况：supp比primary比对的更好，应该选supp作为填充读数
'''
# from Utils import Run_for_denovo, Record, write_reg_to_bed, make_dir, run_cmd_ls
# import fasta_parser
## 
logger = logging.getLogger()
Region = namedtuple('Region', ["chr_id", "start", "end"])
Chr_info = namedtuple('Chr_info', ["chr_id", "chr_len"])
class denovo_info():
    def __init__(self, chr_id, start, end, reg_ctg_ls) -> None:
        self.chr_id = chr_id
        self.start = start
        self.end = end
        self.reg = [chr_id, start, end]
        self.reg_ctg_ls = reg_ctg_ls
        # self.solve = False
    def write_info(info_ls, fout):
        with open(fout, "w") as f:
            for info in info_ls:
                if len(info.reg_ctg_ls) == 0:
                    f.write("{}\t{}\t{}\tfailed_asm\n".format(info.chr_id, info.start, info.end, ",".join(info.reg_ctg_ls)))
                    pass
                else:
                    f.write("{}\t{}\t{}\t{}\n".format(info.chr_id, info.start, info.end, ",".join(info.reg_ctg_ls)))

def can_solve_bound(rec):
    if rec.operation != "asm_patch":
        return False
    else:
        patch_len = len(rec.info)
        reg_len = rec.end - rec.start
        diff =  abs(patch_len - reg_len)
        print("{}:{}-{}, diff:{}, diff_porotion:{}".format(rec.chr_id, rec.start, rec.end, diff, diff/reg_len))
        if reg_len > 100000:
            if diff / reg_len <= 0.2:
                return True
        elif reg_len <= 100000:
            if diff <= 20000:
                return True
        return False

#### get rec by asm
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

def get_candidate_ctgs(bam_in, reg):
    ##
    left_candidate_ls = []
    right_candidate_ls = []
    candidate_dic = dict()
    candidate_ids_dic = dict()
    l_ids = set()
    r_ids = set()
    min_MQ = 20     # 太严格了，局部组装比对质量可能并不高
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    # chr_id = reg.chr_id
    # chr_len = bam_reader.get_reference_length(chr_id)   # 
    min_patch_len = 5000    # 只保留patch超过一定长度的
    for read in bam_reader.fetch(reg.chr_id, reg.start, reg.end):
        if read.is_secondary or read.mapping_quality < min_MQ:
            continue
        if read.reference_start < reg.start - min_patch_len:    # patch左边界的点
            left_candidate_ls.append(read)
            l_ids.add(read.query_name)
        if read.reference_end > reg.end + min_patch_len:        # patch右边界的点
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
    '''
    选出最佳填充contig
    1、与原始ref交叠长度
    2、与原始ref交叠部分相似度，好像多此一举
    '''
    best_ctg = None
    for ctg in candidate_ls:    # 对候选ctg作过滤，加上对clip的限制，和长度的限制
        if best_ctg != None:
            if ctg.query_alignment_length > best_ctg.query_alignment_length:
                best_ctg/= ctg
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


def get_rec_by_asm(ctg, bam_in, candidate_op_ls, asm_fa_dic, process_id_set, apply_bound): # 传入asm_to_ref.bam, op_Ls
    '''
    注意只能process one ctg, 
    功能：
    对于候选的区域，筛选出所有需要局部组装处理的区域，
    对每个区域进行判定：
    1、能够完全进行跨过->使用局部组装的结果；
    2、否则->将无法组装的区域记录下来failed，送入下一轮的全局混合局部组装
    '''
    print("Process {}".format(ctg))
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    final_op_ls = []
    failed_reg_ls = []
    success_reg_ls = [] # [[chr1, 0, 1000]]
    asm_candidate_ls = []
    res_rec_ls = []
    for rec in candidate_op_ls:
        reg_id = reg_to_id([rec.chr_id, rec.start, rec.end])
        if rec.operation.endswith("reads"):
            rec.add_operation("reads_patch")
            final_op_ls.append(rec)
            # logger.info("{}:{}-{} reads_patch".format(rec.chr_id, rec.start, rec.end))
        elif rec.operation.endswith("skip"):
            final_op_ls.append(rec)
            print("Skip:{}:{}-{}".format(rec.chr_id, rec.start, rec.end))
        elif rec.operation.endswith("patch"):
            final_op_ls.append(rec)
            logger.info("{}:{}-{}, {}".format(rec.chr_id, rec.start, rec.end, rec.operation))
        elif rec.operation.endswith("asm"):   # endswith("asm")  将所有asm区域收集
            # reg_id = rec.chr_id + ":" + str(rec.start) + "-" + str(rec.end)
            # if reg_id in asm_reg_ids:
            if reg_id in process_id_set:
                asm_candidate_ls.append(rec)
            else:
                final_op_ls.append(rec)
            # logger.info("{}:{}-{} asm".format(rec.chr_id, rec.start, rec.end))
        else:
            final_op_ls.append(rec)
            logger.info("{}:{}-{}, {}".format(rec.chr_id, rec.start, rec.end, rec.operation))
            raise("ERROR rec: {}:{}-{}, {}".format(rec.chr_id, rec.start, rec.end, rec.operation)) 
            # continue

    ##  for asm_candidate_ls
    for rec in asm_candidate_ls:
        ## 判断是否几乎 whole contig 被认为是difficult region
        # ctg_len = ctg_len_dic[rec.chr_id]
        ctg_len = bam_reader.get_reference_length(rec.chr_id)
        '''solve difficult contig'''
        if (rec.end - rec.start) / ctg_len > 0.8 or rec.operation == "whole_ctg_asm":    ## 一条染色体绝大部分是复杂区域，整条进行从头组装
            new_rec = Record(rec.chr_id, 0, ctg_len)
            new_rec.add_operation("replace_with_denovo")
            new_rec.add_info(".")
            new_rec.add_patch_id(".")   # 全部丢给从头组装
            logger.info("{} replace_with_denovo".format(rec.chr_id))
            final_op_ls = [new_rec]
            failed_reg_ls.append([rec.chr_id, rec.start, rec.end])
            res_rec_ls.append(new_rec)
            break
        
        '''普通情形'''
        ## 找到pass左边和pass右边的最佳序列
        candidate_dic = get_candidate_ctgs(bam_in, Region(rec.chr_id, rec.start, rec.end))
        left_candidate_ls = candidate_dic["left"]
        right_candidate_ls = candidate_dic["right"]     # 
        ## 暂不将端点考虑进去
        if rec.start < 1000: # 左端  取pass右端
            # print("{}:{}-{} left telomere".format(rec.chr_id, rec.start, rec.end))
            if apply_bound:
                res_rec = solve_t0(right_candidate_ls, rec, asm_fa_dic)
                res_rec.add_loca("left")
                res_rec_ls.append(res_rec)
                if can_solve_bound(res_rec):
                    new_rec = res_rec
                else:
                    new_rec = rec
            else:
                new_rec = rec
        elif rec.end > ctg_len - 1000: # 右端，取pass左端
            # print("{}:{}-{} right telomere".format(rec.chr_id, rec.start, rec.end))
            if apply_bound:
                res_rec = solve_t1(left_candidate_ls, rec, asm_fa_dic)
                res_rec.add_loca("right")
                res_rec_ls.append(res_rec)
                if can_solve_bound(res_rec):
                    new_rec = res_rec
                else:
                    new_rec = rec
            else:
                new_rec = rec
        else:
            # print("{}:{}-{} middle type".format(rec.chr_id, rec.start, rec.end))
            res_rec = solve_t2(left_candidate_ls, right_candidate_ls, rec, asm_fa_dic)
            res_rec.add_loca("mid")
            res_rec_ls.append(res_rec)
            if res_rec.operation == "asm_patch":
                new_rec = res_rec
            else:
                new_rec = rec
        final_op_ls.append(new_rec)
        if new_rec.operation == "asm_patch":
            success_reg_ls.append([new_rec.chr_id, new_rec.start, new_rec.end])
        else:
            failed_reg_ls.append([new_rec.chr_id, new_rec.start, new_rec.end])

    return final_op_ls, failed_reg_ls, success_reg_ls, res_rec_ls

###########
## 
@DeprecationWarning
def select_reads_from_names2(fastq_in, out_dir, read_ids, threads):
    ## 使用seqkit实现 
    '''seqkit grep -f id.txt -j threads seqs.fq.gz -o result.fq.gz'''
    # seqkit_cmd = ["seqkit", "grep", "-f", read_ids, "-j", str(threads), fastq_in, "-o", fastq_out]
    # logger.info("Running: %s", " ".join(seqkit_cmd))
    reg_reads_ids_fn = os.path.join(out_dir, "reg_ids.bed")
    fastq_out = os.path.join(out_dir, "reg.fastq")
    with open(reg_reads_ids_fn, "w") as f:
        for read_id in read_ids:
            f.write("{}\n".format(read_id))
    cmd = ["seqtk", "subseq", fastq_in, reg_reads_ids_fn, ">", fastq_out] 
    subprocess.check_call(" ".join(cmd), shell=True)
    return fastq_out

def select_reads_from_names(fastq_in, fastq_out, read_ids, threads):
    ## 使用seqkit实现 
    '''seqkit grep -f id.txt -j threads seqs.fq.gz -o result.fq.gz'''
    # seqkit_cmd = ["seqkit", "grep", "-f", read_ids, "-j", str(threads), fastq_in, "-o", fastq_out]
    # logger.info("Running: %s", " ".join(seqkit_cmd))
    t0 = time.time()
    cmd = ["seqtk", "subseq", fastq_in, read_ids, ">", fastq_out] 
    subprocess.check_call(" ".join(cmd), shell=True)
    print("Get fastq done, cost:{}s".format(time.time() - t0))
    return fastq_out

def select_reads_from_names3(fastq_in, out_dir, asm_reg:Asm_Region, part, threads):
    ## 使用seqkit实现 
    '''seqkit grep -f id.txt -j threads seqs.fq.gz -o result.fq.gz'''
    # seqkit_cmd = ["seqkit", "grep", "-f", read_ids, "-j", str(threads), fastq_in, "-o", fastq_out]
    # logger.info("Running: %s", " ".join(seqkit_cmd))
    read_ids = asm_reg.get_reads(part)
    reg_id = reg_to_id([asm_reg.chr_id, asm_reg.start, asm_reg.end])
    reg_reads_ids_fn = out_dir + "/" + reg_id + "_" + part + ".ids"
    fastq_out = out_dir + "/" + reg_id + "_" + part + ".fq"
    Asm_Region.write_reads(reg_reads_ids_fn, read_ids)
    cmd = ["seqtk", "subseq", fastq_in, reg_reads_ids_fn, ">", fastq_out] 
    subprocess.check_call(" ".join(cmd), shell=True)
    return fastq_out

def get_reg_ids(bam_in, reg_ls):
    read_ids = set()
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    for reg in reg_ls:
        for read in bam_reader.fetch(reg[0], reg[1], reg[2]):
            read_ids.add(read.query_name)
    return read_ids

def reg_to_id(reg):
    return reg[0] + ":" + str(reg[1]) + "-" + str(reg[2])

def id_to_reg(reg_str:str):
    ctg = reg_str.split(":")[0]
    start, end = reg_str.split(":")[1].split("-")
    start, end = int(start), int(end)
    return [ctg, start, end]

def cal_max_reg_len(config):
    base_max_reg_len = config["denovo_by_reg"]["max_reg_len"]
    res = base_max_reg_len
    reg_len_2 = config["cluster_by_asm_candidate"]["reg_len_2"]    # 20000
    reg_len_3 = config["cluster_by_asm_candidate"]["reg_len_3"]    # 1000000
    reg_expand_1 = config["cluster_by_asm_candidate"]["reg_expand_1"] # 10000   # 对所有候选区域进行的扩充
    reg_expand_2 = config["cluster_by_asm_candidate"]["reg_expand_2"] # 20000
    reg_expand_3 = config["cluster_by_asm_candidate"]["reg_expand_3"] # 50000
    if base_max_reg_len < reg_len_2:
        res += 2 * reg_expand_1
    elif base_max_reg_len < reg_len_3:
        res += 2 * reg_expand_2
    else:
        res += 2 * reg_expand_3
    return res

def contig_pool(samfile):
    contig_len={}
    for (ref,lens) in zip(samfile.references,samfile.lengths):
        contig_len[ref]=lens
    return contig_len

def is_bound(reg, contig_len_dic):
    if reg[1] < 1000 or contig_len_dic[reg[0]] - reg[2] < 1000:
        print("{}:{}-{} is bound".format(reg[0], reg[1], reg[2]))
        return True
    return False

def is_abnormal_ctg(dpinfo, params):
    min_cov_ratio, dp_upper_bound, dp_lower_bound = params["min_cov_ratio"], params["dp_upper_bound"], params["dp_lower_bound"]
    if dpinfo["cov_ratio"] < min_cov_ratio or dpinfo["chr_avg_dp"] / dpinfo["whole_dp"] < dp_lower_bound or dpinfo["chr_avg_dp"] / dpinfo["whole_dp"] > dp_upper_bound:
        return True
    return False

def filter_empty_ctg(whole_ctg_asm_reg_ls, reg_read_ids_dic):
    ls = []
    for reg in whole_ctg_asm_reg_ls:
        if is_empty_reg(reg, reg_read_ids_dic):
            print("Empty ctg: {}".format(reg))
        else:
            ls.append(reg)
    return ls
def is_empty_reg(reg, reg_read_id_dic):
    reg_id = reg_to_id(reg)
    reg_reads_ids = reg_read_id_dic[reg_id]
    if len(reg_reads_ids) == 0: return True
    return False

def judge_(rec_ls, read: pysam.AlignedSegment):
    '''  rec: left or right half solve'''
    if len(rec_ls) == 0:
        return False
    if rec_ls[0] == "replace_with_denovo": return True
    ## 
    span_reg = [read.reference_start, read.reference_end]
    for rec in rec_ls:
        # if read.
        pass

    return

def contig_filter(filter_type, fa_dic, bam_in, rec_dic, config): # rec_dic
    '''
    contig 过滤器，根据组装得到的contig在参考上的位置，判断是否为多余的组装碎片。
    扫描asm_to_ref文件，保留primary进行判断。不知道是不是可行
    根据比对位置进行判断；若：1、good region，丢弃contig; 2、bad region; 判断现有的方法能够做到什么程度，若基本完成，丢弃，否则保留。
    '''
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    ls = [] # record contig type
    new_fa_dic = {}
    unmapped_ls = []
    unsolve_dic = {}
    for chr_id in bam_reader.references:
        ctg_unsolve_ls = []
        for rec in ctg_rec_ls:
            if rec.operation == "read_patch":continue
            elif rec.operation == "asm_patch":  # left or right or mid, keep left and right
                if rec.loca == "left" or rec.loca == "right":
                    ctg_unsolve_ls.append(rec)
                else: continue
            else:   # replace_with_denovo / N_fill / asm_patch,N_fill,N_fill.....
                ctg_unsolve_ls.append(rec)
        unsolve_dic[chr_id] = ctg_unsolve_ls

    for read in bam_reader.fetch(until_eof=True):
        if read.is_supplementary or read.is_secondary: continue
        if read.is_unmapped:
            unmapped_ls.append([read.query_name, "unmapped"])
            continue
        ctg_rec_ls = rec_dic[chr_id]
        ctg_unsolve_ls = unsolve_dic[chr_id] 
        ## primary read
        
        # 
    
    pass

def filter2(ctg_ls, fa_dic, bam_in, config):
    ''' 只保留比对到该条contig上的序列，别的不保留了'''
    new_fa_dic = {}
    min_len = config["min_len"]
    # unmapped_ls = []
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    for read in bam_reader.fetch(until_eof=True):
        if read.is_supplementary or read.is_secondary: continue
        if read.is_unmapped:
            new_fa_dic[read.query_name] = fa_dic[read.query_name]
        if read.reference_name in ctg_ls and len(fa_dic[read.query_name]) > min_len:    # 
            new_fa_dic[read.query_name] = fa_dic[read.query_name]
    return new_fa_dic

def extract_reg_ids(bam_in, asm_reg_ls):
    t0 = time.time()
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    new_asm_reg_ls = []
    for asm_reg in asm_reg_ls:
        all_reads = set()
        l_reads = set()
        r_reads = set()
        for read in bam_reader.fetch(asm_reg.chr_id, asm_reg.start, asm_reg.end):   # 加入了所有的序列
            all_reads.add(read.query_name)
        if asm_reg.loca == "mid":
            for read in bam_reader.fetch(asm_reg.chr_id, asm_reg.start, asm_reg.start + 10):    # 跨过左端点的读数
                l_reads.add(read.query_name)
            for read in bam_reader.fetch(asm_reg.chr_id, asm_reg.start, asm_reg.start + 10):    # 跨过右端点的读数
                r_reads.add(read.query_name)
        ## 
        if asm_reg.loca == "mid":
            asm_reg.add_l_reads(l_reads)
            asm_reg.add_r_reads(r_reads)
            asm_reg.add_all_reads(all_reads)
        else:
            asm_reg.add_all_reads(all_reads)
        new_asm_reg_ls.append(new_asm_reg_ls)
    print("Extract read ids cost: {}".format(time.time() - t0))
    return new_asm_reg_ls

def extract_reg_reads(fq_in, bam_in, asm_reg_ls, whole_ctg_asm_reg_ls, threads, out_dir):
    '''提取子区域的读数
    1、get all reg reads
    '''
    # get reg read ids
    asm_reg_ls = extract_reg_ids(bam_in, asm_reg_ls)
    # get reg reads
    task_ls = []
    pool = Pool(processes=threads)
    for task in task_ls:
        pool.apply_async(select_reads_from_names3, args=())    # fastq_in, out_dir, asm_reg:Asm_Region, part, threads
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕

def run_reg_denovo(region_ls, out_dir, reg_read_id_dic, candidate_op_dic, fastq_in, bam_in, threads, data_type, reference, Depthrec_dic, config):
    t0 = time.time()
    asm_dir = out_dir + "/asm"
    make_dir(out_dir)
    make_dir(asm_dir)
    ls1 = []
    semi_op_ls, failed_reg_ls, success_reg_ls, res_rec_ls = [], [], [], []
    process_regid_set = set()   # 进行这一步处理的reg_id
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    contig_len_dic = contig_pool(bam_reader)
    max_reg_len = cal_max_reg_len(config)
    dp_upper_bound = config["denovo_by_reg"]["dp_upper_bound"]
    block_high_dp_bound = config["denovo_by_reg"]["block_high_dp_bound"]
    print("max_reg_len: {}".format(max_reg_len))
    apply_bound = config["denovo_by_reg"]["apply_bound"]
    if apply_bound: print("Will reg denovo apply on bound!!!")
    dp_ls = []
    for idx, reg in enumerate(region_ls): 
        dp_ls.append([Depthrec_dic[reg[0]].get_block_dp(reg), Depthrec_dic[reg[0]].get_block_high_dp(reg)])
        print("{}, avg_dp:{}, high_block_dp:{}".format(reg, str(dp_ls[idx][0]), str(dp_ls[idx][1])))
    for idx, reg in enumerate(region_ls):   # all denovo reg
        if reg[2] - reg[1] > max_reg_len:   # 超过大小，排除
            print("Longger than the maxreglen: {}".format(max_reg_len))
            failed_reg_ls.append(reg)
        elif (reg[2] - reg[1]) / contig_len_dic[reg[0]] > 0.7:  # whole ctg
            print("Too long denovo reg: {}, chr_len: {}".format(reg, contig_len_dic[reg[0]]))
            failed_reg_ls.append(reg)
        elif not apply_bound and is_bound(reg, contig_len_dic):    # 是边界且不需要边界，排除 
            print("{} is bound and donot need bound".format(reg))
            failed_reg_ls.append(reg)
        elif dp_ls[idx][0] / Depthrec_dic[reg[0]].whole_dp > dp_upper_bound or dp_ls[idx][1] / Depthrec_dic[reg[0]].whole_dp > block_high_dp_bound:
            print("Too high dp reg: {}, avg_dp: {}, block_high_dp: {}".format(reg, dp_ls[idx][0], dp_ls[idx][1]))
            failed_reg_ls.append(reg)
        elif is_empty_reg(reg, reg_read_id_dic):
            print("Empty reg: {}".format(reg))
            failed_reg_ls.append(reg)
        else:   # 剩下的区间，用于reg denovo
            ls1.append(reg)
            process_regid_set.add(reg_to_id(reg))
    
    print("asm_regions: {}".format(region_ls))
    print("ls1: {}".format(ls1))
    print("Do not perform on: {}".format(failed_reg_ls))
    ## choose for denovo
    task_ls = []    # [[]]
    for reg in ls1:
        reg_id = reg_to_id(reg)
        reg_work_dir = os.path.join(asm_dir, reg_id)
        make_dir(reg_work_dir)
        reg_reads_ids = reg_read_id_dic[reg_id]
        reg_reads_ids_fn = os.path.join(reg_work_dir, "reg_denovo_ids.bed")
        reg_fastq = os.path.join(reg_work_dir, "reg_denovo.fastq")
        with open(reg_reads_ids_fn, "w") as f:
            for read_id in reg_reads_ids:f.write("{}\n".format(read_id))
        task_ls.append([reg, reg_work_dir, reg_fastq, reg_reads_ids_fn])
        
    ## extracrt_fq_parallel
    t1 = time.time()
    pool = Pool(processes=threads)
    for task in task_ls:
        pool.apply_async(select_reads_from_names, args=(fastq_in, task[2], task[3], 1))
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    print("Extract reg read done, cost {}s".format(time.time() - t1))
    ### reg denovo Assembly
    denovo_info_f = os.path.join(out_dir, "denovo_info.bed")
    reg_denovo_info_ls = []
    t2 = time.time()
    denovo_asm_dic = {}
    # parallel
    if config["denovo_by_reg"]["parallel"]:
        results = []
        pool_threads = threads // 4 + 2
        pool = Pool(processes=pool_threads)
        for task in task_ls:
            reg = task[0]
            genome_size = min((reg[2] - reg[1]) * 2, reg[2] - reg[1] + 40000)
            print("Add {} to Pool".format(task[-1]))
            res = pool.apply_async(Run_for_denovo, args=(task[2], task[1], threads, genome_size, data_type, config))    # Args: fq_in, work_dir, threads, genome_size, data_type, config
            results.append(res)
        pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
        pool.join() # 等待进程池中的所有进程执行完毕
        print("Reg assembly finished, cost {}s".format(time.time() - t2))
        for idx, res in enumerate(results):     # Get results
            task = task_ls[idx]
            reg = task[0]
            reg_denovo_asm_out = res.get()
            reg_denovo_asm_dic = fasta_parser.read_sequence_dict(reg_denovo_asm_out)
            reg_ctg_ls = []
            for ctg,seq in reg_denovo_asm_dic.items():
                new_ctg = "reg" + str(idx) + "_" + ctg
                denovo_asm_dic[new_ctg] = seq
                reg_ctg_ls.append(new_ctg)
            reg_denovo_info_ls.append(denovo_info(reg[0], reg[1], reg[2], reg_ctg_ls))
    # 串行
    else:
        for idx, task in enumerate(task_ls):
            reg = task[0]
            genome_size = min((reg[2] - reg[1]) * 2, reg[2] - reg[1] + 40000)
            reg_denovo_asm_out = Run_for_denovo(task[2], task[1], threads, genome_size, data_type, config)  # Args: fq_in, work_dir, threads, genome_size, data_type, config
            reg_denovo_asm_dic = fasta_parser.read_sequence_dict(reg_denovo_asm_out)
            reg_ctg_ls = []
            for ctg,seq in reg_denovo_asm_dic.items():
                new_ctg = "reg" + str(idx) + "_" + ctg
                denovo_asm_dic[new_ctg] = seq
                reg_ctg_ls.append(new_ctg)
            reg_denovo_info_ls.append(denovo_info(reg[0], reg[1], reg[2], reg_ctg_ls))
            print("Run {} assembly done".format(reg_to_id(reg)))
    ## 
    denovo_info.write_info(reg_denovo_info_ls, denovo_info_f)
    print("Reg assembly finished, cost {}s".format(time.time() - t2))
    ## 
    denovo_asm_fa = os.path.join(out_dir, "denovo.fa")
    fasta_parser.write_fasta_dict(denovo_asm_dic, denovo_asm_fa)

    ### mapping to ref
    sorted_bam_out = os.path.join(out_dir, "denovoasm_to_ref.sorted.bam")
    mapping_cmd_ls = ["minimap2", "-ax", "asm20", "-t", str(threads), reference, denovo_asm_fa, "|", "samtools", "sort", "-O", "BAM", "-@", str(threads), "-o", sorted_bam_out, "&&", "samtools index -@", str(threads), sorted_bam_out] 
    subprocess.check_call(" ".join(mapping_cmd_ls), shell=True)
    ### apply SV consensus by asm, just get rec in this process
    print("Start reg denovo patch")
    semi_candidate_op_bed = out_dir + "/semi_candidate_op.bed"
    res_rec_bed = out_dir + "/denovo_rec.bed"
    failed_reg_file = out_dir + "/failed.bed"
    success_reg_file = out_dir + "/success.bed"
    all_chrs = bam_reader.references
    pool = Pool(processes=threads)
    results = [pool.apply_async(get_rec_by_asm, args=(ctg, sorted_bam_out, candidate_op_dic[ctg], denovo_asm_dic, process_regid_set, apply_bound)) for ctg in all_chrs]
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    for res in results:
        ctg_semi_op_ls, ctg_failed_reg_ls, ctg_success_reg_ls, ctg_res_rec_ls = res.get()
        semi_op_ls.extend(ctg_semi_op_ls)
        failed_reg_ls.extend(ctg_failed_reg_ls) 
        success_reg_ls.extend(ctg_success_reg_ls)
        res_rec_ls.extend(ctg_res_rec_ls)
    
    Record.write_record(semi_op_ls, semi_candidate_op_bed)
    Record.write_rec_short(res_rec_ls, res_rec_bed)
    write_reg_to_bed(failed_reg_file, failed_reg_ls)
    write_reg_to_bed(success_reg_file, success_reg_ls)
    ## 
    semi_candidate_rec_dic = defaultdict(list)
    for rec in semi_op_ls:
        semi_candidate_rec_dic[rec.chr_id].append(rec)
    
    print("run_reg_denovo done, cost {}s".format(time.time() - t0))
    return semi_candidate_rec_dic, failed_reg_ls, success_reg_ls

def run_ctg_denovo(reg_ls, out_dir, reg_read_id_dic, fastq_in, reference, threads, data_type, config):
    ''' '''
    t0 = time.time()
    ## 
    task_ls = []
    for reg in reg_ls:  # get task_ls
        ctg = reg[0]
        ctg_workdir_dir = out_dir + "/" + ctg
        make_dir(ctg_workdir_dir)
        reg_id = reg_to_id(reg)
        ctg_reads_ids = reg_read_id_dic[reg_id]
        ctg_reads_ids_fn = os.path.join(ctg_workdir_dir, "reg_denovo_ids.bed")
        ctg_fastq = os.path.join(ctg_workdir_dir, "reg_denovo.fastq")
        with open(ctg_reads_ids_fn, "w") as f:
            for read_id in ctg_reads_ids:f.write("{}\n".format(read_id))
        task_ls.append([reg, ctg_workdir_dir, ctg_fastq, ctg_reads_ids_fn])
    
    ### get reads
    t1 = time.time()
    try:
        if len(reg_ls) != 0:
            pool = Pool(processes=threads)
            for task in task_ls:
                pool.apply_async(select_reads_from_names, args=(fastq_in, task[2], task[3], 1))
            pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
            pool.join() # 等待进程池中的所有进程执行完毕
    except:
        raise("run_ctg_denovo get reads failed!!!")
    print("Extract ctg read done, cost {}s".format(time.time() - t1))

    ### assembly
    fa_dic = {}
    ctg_denovo_info_ls = []
    ctg_denovo_info_f = out_dir + "/denovo_info.bed"
    ctg_denovo_fa = out_dir + "/ctg_asm.fa"
    # parallel perform
    results = []
    pool = Pool(processes=threads)
    for task in task_ls:
        task_threads = threads
        ctg = task[0][0]
        genome_size = task[0][2] - task[0][1]
        print("Add {} to Pool".format(task[-1]))
        results.append(pool.apply_async(Run_for_denovo, args=(task[2], task[1], task_threads, genome_size, data_type, config)))
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    for idx, res in enumerate(results):
        task = task_ls[idx]     # 与添加的时候顺序一致
        ctg = task[0][0]
        ctg_denovo_asm_out = res.get()
        ctg_denovo_asm_dic = fasta_parser.read_sequence_dict(ctg_denovo_asm_out)
        new_ctg_ls = []
        for asm_ctg,seq in ctg_denovo_asm_dic.items():
            new_ctg = ctg + "_" + asm_ctg   # 
            fa_dic[new_ctg] = seq
            new_ctg_ls.append(new_ctg)
        ctg_denovo_info_ls.append(denovo_info(task[0][0], task[0][1], task[0][2], new_ctg_ls))

    # 普通串行
    # for task in task_ls:
    #     ctg = task[0][0]
    #     genome_size = task[0][2] - task[0][1]
    #     ctg_denovo_asm_out = Run_for_denovo(task[2], task[1], threads, genome_size, data_type, config)
    #     ctg_denovo_asm_dic = fasta_parser.read_sequence_dict(ctg_denovo_asm_out)
    #     new_ctg_ls = []
    #     for asm_ctg,seq in ctg_denovo_asm_dic.items():
    #         new_ctg = ctg + "_" + asm_ctg   # 
    #         fa_dic[new_ctg] = seq
    #         new_ctg_ls.append(new_ctg)
    #     ctg_denovo_info_ls.append(denovo_info(task[0][0], task[0][1], task[0][2], new_ctg_ls))
    ##
    denovo_info.write_info(ctg_denovo_info_ls, ctg_denovo_info_f)
    fasta_parser.write_fasta_dict(fa_dic, ctg_denovo_fa)
    print("Ctg denovo asm Done, cost {}s".format(time.time() - t0))
    
    ### filter
    
    if config["apply_contig_filter"]:
        ctg_ls = [reg[0] for reg in reg_ls]
        ### mapping
        bam_out = out_dir + "/ctg_to_ref.sort.bam"
        minimap2_cmd_ls = ["minimap2", "-ax", "asm20", "-t", str(threads), reference, ctg_denovo_fa, \
            "|", "samtools", "sort", "-O", "BAM", "-@", str(threads), "-o", bam_out, \
            "&&", "samtools index -@", str(threads), bam_out] 
        run_cmd_ls(minimap2_cmd_ls)
        filtered_fa_dic = filter2(ctg_ls, fa_dic, bam_out, config["contig_filter"])
        filter_ls = filtered_fa_dic.keys()
        print("Filter from: {} \n-> {}".format(fa_dic.keys(), filter_ls))
    else:
        filtered_fa_dic = fa_dic
    
    return filtered_fa_dic

def run_merge_local_denovo(reg_ls, asm_ctg_ls, out_dir, reg_read_id_dic, fastq_in, unmapped_fq, reference, threads, data_type, config):
    # reg_ls, out_dir, reg_read_id_dic, fastq_in, reference, threads, data_type, config
    merge_reg_read_ids_fn = out_dir + "/asm_region_ids.bed"
    merge_reg_fq = out_dir + "/reg.fastq"
    merge_denovo_fq = out_dir + "/denovo.fastq"
    merge_reg_read_ids = set()
    merge_local_reg_ls = []
    denovo_info_ls = []
    genome_size = 0
    for reg in reg_ls:
        if reg[0] in asm_ctg_ls: continue
        reg_id = reg_to_id(reg)
        merge_reg_read_ids.update(reg_read_id_dic[reg_id])
        merge_local_reg_ls.append(reg)     # 添加数据
        genome_size += min((reg[2] - reg[1]) * 2, reg[2] - reg[1] + 40000)
    with open(merge_reg_read_ids_fn, "w") as f:
        for read_id in merge_reg_read_ids: 
            f.write("{}\n".format(read_id))
    select_reads_from_names(fastq_in, merge_reg_fq, merge_reg_read_ids_fn, threads)
    # 读数合并，unmapped_fq + merge_reg_fq -> final_denovo_fq
    reads_merge_cmd = ["cat", unmapped_fq, merge_reg_fq, ">", merge_denovo_fq]
    logger.info("Running: %s", " ".join(reads_merge_cmd))
    subprocess.check_call(" ".join(reads_merge_cmd), shell=True)
    ## Assembly
    logger.info("Apply merge local denovo on: {}".format(merge_local_reg_ls))
    
    merge_denovo_asm_out = Run_for_denovo(merge_denovo_fq, out_dir, threads, genome_size, data_type, config)
    if merge_denovo_asm_out == None:    # 组装失败
        print("Merge Assembly failed, Done")
        exit(0)
    fa_dic = fasta_parser.read_sequence_dict(merge_denovo_asm_out)
    ## mapping
    logger.info("****************** map denovo asm to reference start ******************")
    bam_out = out_dir + "/denovoasm_to_ref.sort.bam"
    minimap2_cmd_ls = ["minimap2", "-ax", "asm20", "-t", str(threads), reference, merge_denovo_asm_out, \
        "|", "samtools", "sort", "-O", "BAM", "-@", str(threads), "-o", bam_out, \
        "&&", "samtools index -@", str(threads), bam_out] 
    run_cmd_ls(minimap2_cmd_ls)
    logger.info("******************map denovo asm to reference Done!!!******************")
    ## 
    return fa_dic, bam_out

def test_denovo():
    import yaml
    config_file = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/Configs/Config.yaml"
    with open(config_file, "r") as f: # config参数获取
        config = yaml.safe_load(f.read())     # 获取部分参数
    threads = 10
    data_dir = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/yeast_ont/my_pipe3/step3_SV_consensus/denovo_asm1/asm"
    out_dir = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/yeast_ont/my_pipe3/step3_SV_consensus/denovo_asm1/test"
    data_type = "ont" 
    reg_ls = [['chrXIII', 0, 30000], ['chrXIII', 133500, 206000], ['chrXIII', 876000, 903361], ['chrXI', 0, 13500], ['chrXI', 649500, 692235], ['chrVIII', 0, 25100], ['chrVIII', 150000, 262000], ['chrVIII', 518000, 547529], ['chrIX', 0, 41500], ['chrIX', 410000, 431184], ['chrXV', 0, 29500], ['chrXV', 918500, 983000], ['chrXV', 993000, 1048295], ['chrV', 0, 27000], ['chrV', 404000, 523500], ['chrV', 557500, 575802], ['chrX', 0, 29000], ['chrX', 327500, 398500], ['chrX', 679000, 725652], ['chrVI', 0, 62500], ['chrVI', 266500, 285938], ['chrXII', 0, 62500], ['chrXII', 439500, 478500], ['chrXII', 687000, 717500], ['chrXII', 879700, 970000], ['chrXII', 1016000, 1043178], ['chrXVI', 0, 71000], ['chrXVI', 875500, 901780], ['chrXIV', 0, 58300], ['chrXIV', 694500, 800685], ['chrVII', 0, 43500], ['chrVII', 675000, 856500], ['chrVII', 1008000, 1107723], ['chrIV', 0, 69500], ['chrIV', 468000, 553600], ['chrIV', 939500, 976000], ['chrIV', 1442000, 1497473], ['chrII', 0, 44500], ['chrII', 749000, 800619], ['chrI', 0, 28000], ['chrI', 126000, 197190]]
    threads = 40
    
    ###
    read_num_ls= []
    task_ls = []
    for reg in reg_ls:
        task_threads = 40    # 40: 34s, 8: 46s
        reg_id = reg_to_id(reg)
        reg_out_dir = out_dir + "/" + reg_id     # 现在的输出路径
        reg_date_dir = data_dir + "/" + reg_id
        make_dir(reg_date_dir)
        reg_fq_in = reg_date_dir + "/" + "reg_denovo.fastq"
        genome_size = min((reg[2] - reg[1]) * 2, reg[2] - reg[1] + 40000)
        reads_num = len(open(reg_date_dir + "/" + "reg_denovo_ids.bed",'r').readlines())
        task_ls.append([reg_fq_in, reg_out_dir, task_threads, genome_size, data_type, config["step3"], reads_num, reg])
        read_num_ls.append(reads_num)
    
    t0 = time.time()
    results =[]
    pool = Pool(processes=threads)
    for task in task_ls:
        print("Add {} to Pool".format(task[-1]))
        res = pool.apply_async(Run_for_denovo, args=(task[0], task[1], task[2], task[3], task[4], task[5]))
        results.append(res)
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    for idx, res in enumerate(results):
        print(task_ls[idx][1],res.get())
    t1 = time.time()
    for task in task_ls:    # 串行执行：147.79468488693237
        print("start {}".format(task[-1]))
        Run_for_denovo(task[0], task[1], task[2], task[3], task[4], task[5])
    # Run_for_denovo(fq_in, work_dir, threads, genome_size, data_type, config)
    t2 = time.time()
    print("Time: \n", t1 - t0, "\n", t2 - t1)
    pass

if __name__ == "__main__":
    test_denovo()
    # import sys 
    # sys.path.append("/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe")
    # dpinf_file = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/yeast_ont/my_pipe3/step2_candidate_regions/dp_info/dpinfo.bed"
    # from find_candidate_regions import find_reg_by_depth
    # dpinfo_dic = find_reg_by_depth.Depth_info.read_dp_info(dpinf_file)
    # params = {"min_cov_ratio": 0.2, "dp_upper_bound": 2, "dp_lower_bound": 0.2}
    # print(dpinfo_dic["chrMT"], params)
    # print(is_abnormal_ctg(dpinfo_dic["chrIII"], params))
    pass