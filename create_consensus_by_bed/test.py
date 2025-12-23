import pysam
from Bio import Seq
def fun1():
    fasta_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/test_flye_hapdup/NC_060930.1/assembly.fasta"
    ctg_ls = "contig_1,contig_11,contig_12,contig_14,contig_15,contig_16,contig_17,contig_19,contig_2,contig_20,contig_21,contig_23,contig_27,contig_28,contig_29,contig_3,contig_32,contig_35,contig_41,contig_43,contig_45,contig_52,contig_53,contig_54,contig_55,contig_58,contig_60,contig_63,contig_64,contig_66,contig_7,contig_8,contig_9".split(",")
    fasta_reader = pysam.FastaFile(fasta_in)
    for ctg in ctg_ls:
        # ref_seq = pysam.FastaFile(fasta_in).fetch(ctg)
        print(fasta_reader.get_reference_length(ctg))
        
    print("SV consensus on: {}".format(ctg_ls))

    # ctg_ls = "contig_1,contig_11,contig_12,contig_14,contig_15,contig_16,contig_17,contig_19,contig_2,contig_20,contig_21,contig_23,contig_27,contig_28,contig_29,contig_3,contig_32,contig_35,contig_41,contig_43,contig_45,contig_52,contig_53,contig_54,contig_55,contig_58,contig_60,contig_63,contig_64,contig_66,contig_7,contig_8,contig_9,"
    # print(ctg_ls.strip(",").split(","))

    len_dic = dict()
    for ctg in ctg_ls:
        len_dic[ctg] = fasta_reader.get_reference_length(ctg)
    print(len_dic)

def fun2():
    bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/flye/flye_to_ref.sort.bam"
    # bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/test_local_asm/wtdbg2/asm_to_ref.sort.bam"
    bam_reader = pysam.AlignmentFile(bam_in, "rb")
    reg = "NC_001139.9:730000-737000"
    read_spans = []
    for read in bam_reader.fetch(region=reg):
        print("reference_span:{}-{}".format(read.reference_start, read.reference_end))
        if read.is_supplementary:
            print("read_span:{}-{}".format(read.query_alignment_start + read.cigartuples[0][1], read.cigartuples[0][1]+read.query_alignment_end))
            read_spans.append([read.query_alignment_start + read.cigartuples[0][1], read.cigartuples[0][1]+read.query_alignment_end])
        else:
            print("read_span:{}-{}".format(read.query_alignment_start, read.query_alignment_end))
            read_spans.append([read.query_alignment_start, read.query_alignment_end])
    print("read_gap:{}".format(read_spans[1][0]-read_spans[0][1]))
    # pysam.depth("-a", "-J", "-Q", str(min_MQ), "-r", ctg, bam_in)
    res = pysam.depth("-a", "-J", "-Q", "20", "--excl-flags", "256", "--excl-flags", "2048", "-r", "NC_001139.9", bam_in)
    return
# fun2()
def test_cigar_reverse():
    from Bio import Seq
    reg = "NC_001145.3"
    bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/flye/flye_to_ref.sort.bam"
    bam_reader = pysam.AlignmentFile(bam_in, "rb")
    check_fa = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/flye/flye/assembly.fasta"   # 验证用的fasta文件
    fasta_reader = pysam.FastaFile(check_fa)
    # for read in bam_reader.fetch(region="NC_001139.9:6165-6239"):
    for read in bam_reader.fetch():
        if read.is_secondary or read.is_supplementary:
            continue
        # print(read.flag & 0x1)
        rec_seq = read.query_sequence
        check_seq = fasta_reader.fetch(read.query_name)
        
        if read.is_reverse:
            print(read.query_name, "is_reverse")
            # print(read.flag)
            if rec_seq == check_seq:
                print("按照原始的来")
            elif rec_seq == Seq.reverse_complement(check_seq):
                print("反向互补")
            else:
                print("不知道")
# test_cigar_reverse()

import re
def solve_by_reads(ctg, start, end, support_ls):
    '''选取最优的一条序列，我们这里根据比对长度选择比对最长的一条序列'''
    if end <= start:
        print("{}:{}-{} ERROR REGION !!!".format(ctg, start, end))
        return ""
    solve_seq = ""
    ## 选择最佳读数
    max_len = -1
    best_read = None
    for read in support_ls:
        if read.query_alignment_length > max_len:
            max_len = read.query_alignment_length
            best_read = read
    ## 根据最佳读数选择   //读数reverse  不用作反转
    # target_start = 0
    # target_end = 0
    read_pos = 0
    ref_pos = read.reference_start  # 0-based
    cigar = best_read.cigarstring
    candidate_pos = []
    for token in re.findall("[\d]{0,}[A-Z]{1}", cigar):
        op = token[-1]
        op_len = int(token[:-1])
        if op in "S":   # H (hard clip未记录在sequence中，不用管  H continue)
            read_pos += op_len
        if op in "M=X": # 同时加
            if len(candidate_pos) == 0:
                if ref_pos + op_len >= start:    # 找到左端点
                    candidate_pos.append(read_pos + (start - ref_pos))
            if len(candidate_pos) == 1:
                if ref_pos + op_len >= end:     # 找到右端点
                    candidate_pos.append(read_pos + (end - ref_pos))
                    break
            read_pos += op_len
            ref_pos += op_len
        if op == "D":
            if len(candidate_pos) == 0:
                if ref_pos + op_len >= start:    # 找到左端点
                    candidate_pos.append(read_pos)
            if len(candidate_pos) == 1:
                if ref_pos + op_len >= end:     # 找到右端点
                    candidate_pos.append(read_pos)
                    break
            ref_pos += op_len
        if op == "I":
            read_pos += op_len
        if len(candidate_pos) >= 2:
            break
    if len(candidate_pos) < 2:
        print("{}:{}-{} ERROR!!!".format(ctg, start, end))
        return ""
    target_start, target_end = candidate_pos[0], candidate_pos[1]
    solve_seq = best_read.query_sequence[target_start:target_end]
    return solve_seq


# print("other_clip_reg"[:-9])

def judge(target_interval, interval_Ls):
    l = target_interval[0]
    r = target_interval[1]
    for ival in interval_Ls:
        if l >= ival[0] and r <= ival[1]:
            return True
    return False
def check_region():
    raw_bed_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/candidate_regions2/candidate.bed"
    op_bed = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/test_pipe3/candidate_op.bed"
    raw_reg_ls = []
    with open(raw_bed_in, "r") as f:
        for line in f:
            fields = line.strip().split()
            raw_reg_ls.append((fields[0], int(fields[1]), int(fields[2])))
    op_reg_dict = dict()
    with open(op_bed, "r") as f:
        for line in f:
            fields = line.strip().split()
            # op_reg_dict[fields[0]] = op_reg_dict.get(fields[0], list()).append((int(fields[1]), int(fields[2])))
            if fields[0] in op_reg_dict:
                op_reg_dict[fields[0]].append((int(fields[1]), int(fields[2])))
            else:
                op_reg_dict[fields[0]] = [(int(fields[1]), int(fields[2]))]
    print(op_reg_dict)
    for reg in raw_reg_ls:
        if not judge((reg[1], reg[2]), op_reg_dict.get(reg[0], [])):
            print("{}:{}-{} ERROR".format(reg[0], reg[1], reg[2]))
# check_region()        

class Record():
    def __init__(self, chr_id, start, end) -> None:
        self.chr_id = chr_id
        self.start = start
        self.end = end
        self.operation = None
        self.info = None
    def add_info(self, info):
        self.info = info
    def add_operation(self, operation):
        self.operation = operation

def select_reads_from_names(fastq_in, fastq_out, readset:set):
    # reads_names = get_reads_names(bam)
    # readset = set(reads_names) # 去重一下

    fout = open(fastq_out,'w')
    with open(fastq_in,'r') as fin:
        for line in fin:
            if line.startswith('@'):
                if line[1:].strip().split()[0] in readset:
                    fout.write(line)
                    fout.write(fin.readline())
                    fout.write(fin.readline())
                    fout.write(fin.readline())
                    readset.remove(line[1:].strip().split()[0])     # 从中去除
                else:
                    fin.readline()
                    fin.readline()
                    fin.readline()
    fout.close()

def func():
    candidate_op_ls = []
    with open("/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/test_pipe3/candidate_op.bed", "r") as f:
        for line in f:
            # if not line:
            #     continue
            fields = line.strip().split()
            rec = Record(fields[0], int(fields[1]), int(fields[2]))
            rec.add_operation(fields[3])
            candidate_op_ls.append(rec)
    bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/aln.sort.bam"
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    unmapped_ids = set()
    for read in bam_reader.fetch(until_eof=True):
        if read.is_unmapped:
            unmapped_ids.add(read.query_name)
    print("unmapped reads num: {}".format(len(unmapped_ids)))
 
    ## get reads of asm_reagion
    asm_regions = [(rec.chr_id, rec.start, rec.end) for rec in candidate_op_ls if rec.operation.endswith("asm")]
    print("Candidate local denovo regions:", asm_regions)
    asm_regions_ids = set()
    for reg in asm_regions:
        for read in bam_reader.fetch(reg[0], reg[1], reg[2]):
            asm_regions_ids.add(read.query_name)
    print("asm_regions reads num:".format(len(asm_regions_ids)))    
    denovo_ids = unmapped_ids.union(asm_regions_ids)
    print("Total unmapped and asm_regions reads num:", len(denovo_ids))

    fq_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/fastq/CRR198362/CRR198362_1.fastq"
    fq_out = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/test_pipe3/denovo.fastq"
    select_reads_from_names(fq_in, fq_out, denovo_ids)

# func()

def test_denovo_asm():
    import os
    out_dir = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/test_pipe3"
    denovo_asm_dir = out_dir +"/denovo_asm"
    if not os.path.exists(denovo_asm_dir):
        os.mkdir(denovo_asm_dir)
    threads = 24
    reads_fq = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/test_pipe3/denovo.fastq"
    out_prefix = denovo_asm_dir+"/out"
    asm_cmd1 = "time wtdbg2 -x preset2 -t " + str(threads) + " -i " + reads_fq +" -fo " + out_prefix
    asm_cmd2 = "time wtpoa-cns -t " + str(threads) + " -i " + out_prefix+".ctg.lay.gz" + " -fo " + out_prefix+".ctg.fa"
    for cmd in [asm_cmd1, asm_cmd2]:

        print("\n{}".format(cmd))
    os.system(asm_cmd1)
    os.system(asm_cmd2)
    return
# test_denovo_asm()
# def mapping_asm_to_ref():
##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################



## 弃置代码
    '''
def get_best_ctg(candidate_ls):   #注意输入的ctg AlignedSegment类型数组
    best_ctg = None
    for ctg in candidate_ls:
        if best_ctg:
            if ctg.query_alignment_length > best_ctg.query_alignment_length:
                best_ctg = ctg
        else:
            best_ctg = ctg
    return best_ctg
    '''

    '''
def solve_t0(candidate_ls, rec, asm_fa_dic):    # 左端点的类型
    # 提供pass右边界的候选ctg, 初始的rec, 以及local_asm的fasta文件

    new_rec = Record(rec.chr_id, rec.start, rec.end)
    best_ctg = None
    max_len = -1
    for ctg in candidate_ls:    # get best contig
        if ctg.query_alignment_length > max_len:
            max_len = ctg.query_alignment_length
            best_ctg = ctg 
    if best_ctg is None:
        new_rec.add_info(".")
        new_rec.add_operation("N_fill")
    else:
        ctg_seq = asm_fa_dic[best_ctg.query_name]
        if best_ctg.is_reverse: # 判断是否反向链
            ctg_seq = Seq.reverse_complement(ctg_seq)

        query_pos = convert_reference_pos_to_raw_pos(best_ctg, rec.end, True)
        info = ctg_seq[:query_pos]    ## 注意是否需要reverse
        new_rec.add_info(info)
        new_rec.add_operation("asm_patch")
    
    info, op = get_info_op(candidate_ls, asm_fa_dic, "right", rec.end)
    new_rec.add_info(info)
    new_rec.add_operation(op)
    return new_rec
    '''

from collections import namedtuple
Region = namedtuple('Region', ["chr_id", "start", "end"])
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
    min_MQ = 20
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    for read in bam_reader.fetch(reg.chr_id, reg.start, reg.end):
        if read.is_secondary or read.mapping_quality < min_MQ:
            continue
        if read.reference_start < reg.start:
            left_candidate_ls.append(read)
            l_ids.add(read.query_name)
        if read.reference_end > reg.end:
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
    direction表示与参考良好match的那部分
    '''
    ## get_best_ctg
    best_ctg = None
    for ctg in candidate_ls:
        if best_ctg != None:
            if ctg.query_alignment_length > best_ctg.query_alignment_length:
                best_ctg = ctg
        else:
            best_ctg = ctg
    ## get info op
    new_info = "."
    op = "N_fill"
    if best_ctg != None:
        op = "asm_patch"
        ctg_seq = asm_fa_dic[best_ctg.query_name]
        if best_ctg.is_reverse: # 判断是否反向链
            ctg_seq = Seq.reverse_complement(ctg_seq)
        if direction == "left":
            query_pos = convert_reference_pos_to_raw_pos(best_ctg, target_pos, True)
            info = ctg_seq[query_pos:]
        elif direction == "right":
            query_pos = convert_reference_pos_to_raw_pos(best_ctg, target_pos, True)
            info = ctg_seq[:query_pos]
        else:
            # raise("ERROR direction of {}:{}-{}".format(rec.chr_id, rec.start, rec.end))
            raise("ERROR direction")
    return new_info, op

def solve_t0(candidate_ls, rec, asm_fa_dic):    # 左端点的类型
    '''
    提供pass右边界的候选ctg, 初始的rec, 以及local_asm的fasta文件
    '''
    new_rec = Record(rec.chr_id, rec.start, rec.end)
    info, op = get_info_op(candidate_ls, asm_fa_dic, "right", rec.end)  # 注意是pass右边界的
    new_rec.add_info(info)
    new_rec.add_operation(op)
    return new_rec


def solve_t1(candidate_ls, rec, asm_fa_dic):     # 右侧的类型
    '''
    提供pass左边界的候选ctg, 初始的rec, 以及local_asm的fasta文件
    '''
    new_rec = Record(rec.chr_id, rec.start, rec.end)
    info, op = get_info_op(candidate_ls, asm_fa_dic, "left", rec.end)  # 注意是pass右边界的
    new_rec.add_info(info)
    new_rec.add_operation(op)
    return new_rec


def solve_t2(left_candidate_ls, right_candidate_ls, rec, asm_fa_dic): # 非端点的类型，中间的区域
    new_rec = Record(rec.chr_id, rec.start, rec.end)

    left_ids = set()
    right_ids= set()
    for read in left_candidate_ls:
        left_ids.add(read.query_name)
    for read in right_candidate_ls:
        right_ids.add(read.query_name)
    common = left_ids.intersection(right_ids)   # 若某条读数同时跨过左右两侧，认为是最佳序列
    
    ## 
    if len(common) >= 1:
        print("{}:{}-{} can span by whole ctg".format(rec.chr_id, rec.start, rec.end))
        left_best_alignment_ls = []
        right_best_alignment_ls = []
        read_to_align_length = dict.fromkeys(common, 0)  # read_id: align_length
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
        best_ctg_id = max(read_to_align_length, key=read_to_align_length.get)  # 获取比对最长的ctg
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
            new_rec.add_info(".")
            new_rec.add_operation("N_fill")
        else:
            ctg_seq = asm_fa_dic[best_ctg_id]
            if l_strand: # 判断是否反向链
                ctg_seq = Seq.reverse_complement(ctg_seq)
            if left_query_pos < right_query_pos:
                info = ctg_seq[left_query_pos:right_query_pos]
                new_rec.add_info(info)
                new_rec.add_operation("asm_patch")
            else:
                print("{}:{}-{} ERROR query position!!!".format(rec.chr_id, rec.start, rec.end))
                new_rec.add_info(".")
                new_rec.add_operation("N_fill")
        return new_rec
    else:   # common == 0
        print("{}:{}-{} can not span by whole ctg".format(rec.chr_id, rec.start, rec.end))
        ## 
        left_info, left_op = get_info_op(left_candidate_ls, asm_fa_dic, "left", rec.start)
        right_info, right_op = get_info_op(right_candidate_ls, asm_fa_dic, "right", rec.end)
        if left_op == "N_fill" and right_op == "N_fill":
            new_rec.add_info(left_info + "," + "N_fill" + "," + right_info)
            new_rec.add_operation(left_op + "," + "." + "," + right_op)
        else:
            new_rec.add_info(left_info + "," + right_info)
            new_rec.add_operation(left_op + "," + right_op)
        return new_rec

        




    # pass

def get_rec_by_asm(bam_in, candidate_op_ls, asm_fa_dic): # 传入asm_to_ref.bam, op_Ls
    '''
    process one ctg, 
    '''
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    final_op_ls = []
    asm_candidate_ls = []
    for rec in candidate_op_ls:
        if rec.operation.endswith("reads"):
            rec.add_operation("reads_patch")
            final_op_ls.append(rec)
        else:
            asm_candidate_ls.append(rec)
    print(len(asm_candidate_ls))
    min_MQ = 20
    ##  for asm_candidate_ls
    for rec in asm_candidate_ls:
        ## 判断是否几乎whole contig 被认为是difficult region
        # ctg_len = ctg_len_dic[rec.chr_id]
        ctg_len = bam_reader.get_reference_length(rec.chr_id)
        if (rec.end - rec.start) / ctg_len > 0.8:    ## too difficult contig???
            info_ls = []
            info_ids = []
            for read in bam_reader.fetch(rec.chr_id):
                if read.is_secondary or read.is_supplementary and read.mapping_quality < min_MQ:
                    continue
                info_ls.append(read.query_sequence)     # 收集所有的primary，将序列记录下来
                info_ids.append(read.query_name)
            new_rec = Record(rec.chr_id, 0, ctg_len)
            new_info = ""
            for info in info_ls:
                new_info += str(info) + ","     # 将局部组装的记录下来
            ## 两类处理方式
            if len(new_info)/ctg_len > 0.5:
                new_rec.add_info(new_info)
                new_rec.add_operation("replace_with_denovo")
                print("{} replace_with_denovo: ".format(rec.chr_id), info_ids)
            else:
                new_rec.add_info(".")
                new_rec.add_operation("keep_ref")
                print("{} denovo failed!!!".format(rec.chr_id))
            final_op_ls = [new_rec]
            return final_op_ls
        
        ## 找到左边和右边的最佳序列
        candidate_dic = get_candidate_ctgs(bam_in, Region(rec.chr_id, rec.start, rec.end))
        left_candidate_ls = candidate_dic["left"]
        right_candidate_ls = candidate_dic["right"]     # 
        ## 
        for ctg in left_candidate_ls:
            pass
        if rec.start < 1000: # 左端  is_telomere
            print("{}:{}-{} left telomere".format(rec.chr_id, rec.start, rec.end))
            new_rec = solve_t0(left_candidate_ls, rec, asm_fa_dic)
        elif rec.end > ctg_len - 1000: # 右端
            print("{}:{}-{} right telomere".format(rec.chr_id, rec.start, rec.end))
            new_rec = solve_t1(right_candidate_ls, rec, asm_fa_dic)
        else:
            new_rec = solve_t2(left_candidate_ls, right_candidate_ls, rec, asm_fa_dic)
        final_op_ls.append(new_rec)

    return final_op_ls

def apply_rec_on_ref(rec_ls, fasta_in, N_fill_size):
    '''per ctg'''
    if len(rec_ls) == 0:
        return [fasta_in]
    if len(rec_ls) == 1:
        if rec.operation == "keep_ref":
            return [fasta_in]
        elif rec.operation == "replace_with_denovo":
            fasta_out = rec.info.split(",")
            return fasta_out

    ## 其他的情形了
    rec_ls.sort(key=lambda rec:int(rec.start))  # 防止有乱序的
    fasta_out = []    # 
    offset = 0
    new_seq = fasta_in
    for rec in rec_ls:
        rec=Record(1,1,1)  #
        start = rec.start + offset
        end = rec.end + offset
        op_ls = rec.operation.split(",")
        info_ls = rec.info.split(",")
        alter_seq = ""
        for i, op in enumerate(op_ls):
            if op == "N_fill":
                alter_seq += ("N" * N_fill_size)
            elif op.endswith("patch"):
                alter_seq += info_ls[i]
            else:
                raise(ValueError)   # 待补充？？？
        offset += len(alter_seq) - (rec.end - rec.start)
        new_seq = new_seq[:start] + alter_seq + new_seq[end:]   # seq modify
    fasta_out.append(new_seq)
    return fasta_out

def parse_fasta(fa):
    with open(fa, "r") as f:
        seqs = {}  # 用一个字典去存储解析出来的序列
        seq_name = None
        for line in f:  # 这里没有直接用read()方法，是为了避免出现文件过大，难以一次性读入内存
            if not line:break
            if line.startswith('>'):
                # 以>开头，说明是含有序列名的行
                seq_name = line.strip().split()[0][1:]
                if seq_name not in seqs:
                    seqs[seq_name] = ''
            else:
                seqs[seq_name] += line.strip('\n')
    return seqs       

def fasta_write(seq_id, seq, out_stream):
    out_stream.write(">"+seq_id+"\n")
    chars_per_line = 60
    for i in range(0, len(seq), chars_per_line):
        out_stream.write(f'{seq[i:i + chars_per_line]}\n')


def test_consensus_on_ref():
    ## 
    bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/test_pipe3/denovo_asm/asm_to_ref.sort.bam"
    op_file_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/test_pipe3/candidate_op.bed"
    denovo_fa = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/test_pipe3/denovo_asm/out.ctg.fa"
    ref_in = ""
    asm_fa_dic = parse_fasta(denovo_fa)
    candidate_op_ls = []
    with open(op_file_in, "r") as f:
        for line in f:
            if line:
                fields = line.strip().split()
                rec = Record(fields[0], int(fields[1]), int(fields[2]))
                rec.add_operation(fields[3])
                rec.add_info(fields[4])
                candidate_op_ls.append(rec)
    new_rec_ls = get_rec_by_asm(bam_in, candidate_op_ls, asm_fa_dic)

    for rec in new_rec_ls:
        print(rec)
    return

# test_consensus_on_ref()

# seq = "TATTCAGCAATAAGAAAATGTGAGCATACTATATATTAATATAAATGATATATCATAAACATAAGCGTATCCAATTTTGACATATCCTTCACGAATATTGTTAGATAACTCTGAACTGTG"
# seq = seq[::-1]
# for i in range(0, len(seq), 60):
#     print(seq[i:i+60])
# print("\n".join([seq[i:i+60]] for i in range(len(seq)//60)))


# seq1 = "AGCTAAA"
# print(Seq.reverse_complement(seq1))

# ls = ["a", "b"]
# dic1 = dict.fromkeys(ls, 0)
# print(dic1)

# dic = {"a":1, "b":3, "c": 0}
# print(max(dic, key=dic.get))
# with open("/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/create_consensus_by_bed/test.fasta", "w") as f:
#     fasta_write("ctg1", seq, f)
#     fasta_write("ctg2", seq1, f)
#     fasta_write("ctg2", seq1, f)

import logging
logger = logging.getLogger()
def _enable_logging(log_file, debug, overwrite):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    if not debug:
        console_log.setLevel(logging.INFO)

    if overwrite:
        open(log_file, "w").close()
    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)

my_log = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/create_consensus_by_bed/test.log"
_enable_logging(my_log, debug=False, overwrite=True)
ls = [0,1]
print("ANNISANO:", 3)
# logger.info("{} replace_with_denovo: {}".format("chr1" ,ls))
# logger.info("Skipped filtering phase")

# preset_type = "preset2"
# threads=24
# out_prefix="aws"
# reads_fq="reads1"
# denovo_asm_dir="denovo"
# asm_cmd1_ls = ["time wtdbg2", "-x", preset_type, "-t", str(threads), "-i ", reads_fq, " -fo ",  out_prefix, ">", denovo_asm_dir+"wtdbg2.log"]
# asm_cmd2_ls = ["time wtpoa-cns", "-t", str(threads), "-i", out_prefix+".ctg.lay.gz", "-fo", out_prefix+".ctg.fa", ">", denovo_asm_dir+"wtpoa-cns.log"]
# asm_cmd1 = " ".join(asm_cmd1_ls)
# asm_cmd2 = " ".join(asm_cmd2_ls)
# # asm_cmd1 = "time wtdbg2 -x preset2 -t " + str(threads) + " -i " + reads_fq +" -fo " + out_prefix + " > wtdbg2.log"
# # asm_cmd2 = "time wtpoa-cns -t " + str(threads) + " -i " + out_prefix+".ctg.lay.gz" + " -fo " + out_prefix+".ctg.fa" + " > wtpoa-cns.log"
# for cmd in [asm_cmd1, asm_cmd2]:
#     logger.info("\n{}".format(cmd))


# import multiprocessing
# import time

# def func(msg):
#     return multiprocessing.current_process().name + '-' + msg

# if __name__ == "__main__":
#     pool = multiprocessing.Pool(processes=4) # 创建4个进程
#     results = []
#     for i in range(10):
#         msg = "hello %d" %(i)
#         results.append(pool.apply_async(func, (msg, )))
#     pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
#     pool.join() # 等待进程池中的所有进程执行完毕
#     print ("Sub-process(es) done.")

#     for res in results:
#         print (res.get())      


# dic = {"a":[], "b":[]}
# dic = dict().fromkeys(["a", "b"], [])
# print(dic)
# dic["a"].append(1)
# print(dic)

# dic1 = {}
# dic1.setdefault("a", []).append(1)
# b = dic1.get("b", [])
# b.append(1)
# print(b)
# from collections import defaultdict
# dic = defaultdict(list)
# print(dic)
# dic["a"].append(1)
# print(dic)
# for key in dic.keys():
#     print(key)


## 


def test():
    import gzip
    import io
    fastq_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/yeast/test_snakemake/ont/CRR198362/CRR198362_1.fastq.gz"
    if fastq_in.endswith("gz"):
        fin = gzip.open(fastq_in, "rt")
    else:
        fin = open(fastq_in, "r")
    for line in fin:
        print(line)
        break
    fin.close()
# test()


def test1():
    # bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/yeast/test_snakemake/snakemake_res/CRR198362_40X/eval/polished_to_ref.sort.bam"
    bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/yeast/test_snakemake/snakemake_res/CRR198362_40X/eval/Y12_to_polished.sort.bam"
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    reg = "chrXIV:1-600000"
    for read in bam_reader.fetch(region=reg):
        print(read.query_name, ":", read.reference_start, "-", read.reference_end)
        # print(convert_reference_pos_to_raw_pos(read, 583273, True))
        print(convert_reference_pos_to_raw_pos(read, 447259, True))
        print()
# test1()


def parse_fasta2(fasta_path):
    f = open(fasta_path, 'r')
    seq_name = None
    seq = None
    for line in f:  # 这里没有直接用read()方法，是为了避免出现文件过大，难以一次性读入内存
        if not line:break
        if line.startswith('>'):
            # 以>开头，说明是含有序列名的行
            seq_name = line[1:].strip('\n')
            if seq:
                yield seq_name,seq
            seq = ''
        else:
            seq += line.strip('\n')
    yield seq_name,seq     
def test2():
    import fasta_parser as fp
    denovo_fa = "/public/home/hpc214712170/test.fasta"    # denovo asm out
    reference_fn = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/melanogaster/REF/ref1/GCA_003401745.1_ASM340174v1_genomic.fna"
    # ref_dic = parse_fasta(reference_fn)
    res = fp.read_sequence_dict(reference_fn)
    # res = fp.read_sequence_dict(denovo_fa)
    # print(res)
    logger.info("Parse fasta Done !!!")
# test2()

