
import re

import pysam
from collections import Counter #导入Counter

match_detail="agctaa,+2ccgg"
st = ''.join(re.split('[\+|\-][0-9]+[ATCGatcg]+', match_detail))
print(st)

dic = Counter(st)
print(dic)
print(dic["a"], dic["bn"])
dic = Counter(match_detail)
print(dic)

match_detail = "agc-3taa,+2ccgg"
nums = re.findall(r'[+-]?\d+', match_detail)
print(nums)
# exit()
# chm13
fastq="/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/chm13/ont/rel6_30X.fastq"
asm="data/flye_ont.fa"
REF="data/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
bam="data/aln2asm.sort.bam"
# thaliana
fastq="/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/thaliana/ont/CRR302667_40X.fastq"
asm="data/flye_ont.fa"
REF="data/GWHBDNP00000000.genome.fasta"
bam="data/aln2asm.sort.bam"

pile_file = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/chm13/my_pipe/step2_candidate_regions/pileup/NC_060925.1:3000-250208550.pileup.txt"
with open(pile_file, "r") as f:
    for line in f:
        record = line.strip().split("\t")
        # print(fields)
        match_detail=record[4]
        st = ''.join(re.split('[\+|\-][0-9]+[ATCGatcg]+', match_detail))    # 用于计算单碱基不一致
        st_counter = Counter(st)
        snp_num = st_counter['a'] + st_counter['A'] + st_counter['g'] + st_counter['G'] + st_counter['c'] + st_counter['C'] + st_counter['t'] + st_counter['T']
        # print(st_counter)
        # if '+' in match_detail:
        #     print(match_detail)
        nums = re.findall(r'[+-]?\d+', match_detail)
        for num in nums:
            if snp_num < 4: break
            if int(num) > 5 or int(num) < -5 and snp_num > 5:
                print(match_detail, st)
                continue
        
        # print(match_detail, " | st:",st)
exit()
'''import re

# 堆积信息字符串
pileup_info = ".+3GAC,.+2GG,.,*..,,.+1G.,.,,.-2GA,.,,-1g,,A,G"

# 提取插入和删除操作的信息
insertions = re.findall(r'\+(\d+)[A-Za-z]+', pileup_info)
deletions = re.findall(r'\-(\d+)[A-Za-z]+', pileup_info)
snp = 1
# 打印结果
ins_num = sum([int(i) for i in insertions])
del_num = sum()
print("插入操作:", insertions, ins_num)
print("删除操作:", deletions)'''


def cal_identity(read:pysam.AlignedRead, reg_start, reg_end):
    '''计算ref_start-ref_end之间的identity'''
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
        # 
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
    return ins_num, del_num # , sum_ins_num, sum_del_num
# bed = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/thaliana/my_pipe/step2_candidate_regions/candidate/false.bed"
bed = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/thaliana/my_pipe/step2_candidate_regions/candidate/merge.clu1.bed"
bam = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/thaliana/simu/aln2simu.sort.bam"
bam_reader = pysam.AlignmentFile(bam)
MIN_ALIGN_LENGTH = 10000
MIN_ALIGN_RATE = 0.95
MIN_MAPPING_QUALITY = 20
with open(bed, "r") as f:
    for line in f:
        fields = line.strip().split("\t")
        ctg, start, end = fields[:3]
        start, end = int(start), int(end)
        span_ls = []
        ctg_len = bam_reader.get_reference_length(ctg)
        
        for read in bam_reader.fetch(ctg, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < MIN_MAPPING_QUALITY or read.query_alignment_length < MIN_ALIGN_LENGTH or (read.query_alignment_length / read.query_length) < MIN_ALIGN_RATE:
                continue
            if read.reference_start < max(0, start - 200) and read.reference_end > min(ctg_len, end + 200):
                span_ls.append([read, 1-read.get_tag("NM")/read.query_length])
        span_ls.sort(key=lambda x:x[1], reverse=True)
        iden_ls = []
        flag = False
        for i in range(len(span_ls)):
            read = span_ls[i][0]
            iden = list(cal_identity(read, start, end))
            # iden_ls.append(iden)
            iden.append(read.query_name)
            iden_ls.append(iden)
            if iden[0] != 0 or iden[1] != 0:flag = True
            # print("{}:{}-{}:".format(ctg, start, end), iden_ls)
        if flag: print("{}:{}-{}:".format(ctg, start, end), iden_ls)
'''for read in bam_reader.fetch("GWHBDNP00000002", 1052886, 10533168):
    if read.query_name == "S2_65347":
        print(read.cigartuples)
        ins_op = 0
        ins_num = 0
        del_op = 0
        del_num = 0
        for op, op_len in read.cigartuples:
            if op == 1:
                if op_len > 30:
                    ins_op += op_len
                    ins_num += 1
            elif op == 2:
                if op_len > 30:
                    del_op += op_len
                    del_num += 1
        print(ins_num, ins_op, del_num, del_op)
        res = cal_identity(read, 10529000, 10533200)
        print(res)'''