import pysam
import re
from collections import namedtuple, defaultdict
Region = namedtuple('Region', ["chr_id", "start", "end"])


########
def cal_avg_depth2(bam_in, REG):  # cal depth of a contig, return list of depth, length == contig length
    reg = REG.chr_id + ":" + str(REG.start) + "-" + str(REG.end)
    '''time samtools depth -a -Q $min_MQ -r $REG ${bam_file} > $out_prefix.depth    
    # 不加-J del处depth为0'''
    min_MQ = 60
    min_len = 5000  # 长度设置的很大，对hifi读数不一定友好
    # reg_len = int(reg.split(":")[2]) - int(reg.split(":")[1])
    depth_ls = []
    depth_stream = pysam.depth("-a", "-Q", str(min_MQ), "-l", str(min_len), "-g", "0x800", "-r", reg, bam_in)    # depth_stream接收stdout
    for line in depth_stream.split("\n"):
        line_ls = line.split()
        if not line_ls:
            continue
        depth_ls.append(int(line_ls[2]))
    return sum(depth_ls) / len(depth_ls)

def is_del_pair(reg1:Region, reg2:Region, bam_in):    # 一前一后两个区间，判断是否为一对pair。根据bam来判断。main for dels
    if reg1.chr_id!= reg2.chr_id:
        return False
    max_connect_len = 20000      # 20000，clip之间的删除
    if reg2.start - reg1.end > max_connect_len: # 能够跨越的删除一般不会超过这个距离
        return False
    ## 根据区间内的信息判断: 1、读数支持：左右两边有相同的读数支持 2、depth支持：将删除不算作depth，计算区间depth
    MIN_DP = 5   # 删除的深度
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in + ".bai")
    ## 计算两个区间之间的avg_dep来判断是否删除？
    '''计算区间之间的区间，为了解决一种特殊的删除情形
    这种删除，在删除区域也能有部分干扰的读数匹配上，两种方式解决：1、提高删除的阈值？2、增加对两侧深度的计算验证
    '''
    ## 第一种计算方式
    # if cal_avg_depth2(bam_in, Region(reg1.chr_id, reg1.end, reg2.start)) < MIN_DP: # reg1.chr_id+":"+str(reg1.end)+"-"+str(reg2.start)
    #     return True
    # else:
    #     return False
    ## 第二种计算方式
    flank_len = 1000
    # min_dep_change = 5
    dep1 = cal_avg_depth2(bam_in, Region(reg1.chr_id, reg1.end, reg2.start)) 
    if dep1 < MIN_DP:
        return True
    dep0 = cal_avg_depth2(bam_in, Region(reg1.chr_id, reg1.start-flank_len, reg1.start)) if reg1.start-flank_len > 0 else cal_avg_depth2(bam_in, Region(reg1.chr_id, 0, reg1.start))
    dep2 = cal_avg_depth2(bam_in, Region(reg1.chr_id, reg2.end, reg2.end+flank_len)) if reg2.end+flank_len < bam_reader.get_reference_length(reg1.chr_id) else cal_avg_depth2(bam_in, Region(reg1.chr_id, reg2.end, bam_reader.get_reference_length(reg1.chr_id)))
    if dep0 > MIN_DP and dep0 < 100 and dep1 / dep0 < 0.5:  # 初衷是想永dep1、dep2代替全局的dep，但现在已经过时了
        return True
    elif dep2 > MIN_DP and dep1 < 100 and dep1 / dep2 < 0.5:
        return True
    else:
        return False

###########
def check_left_clip_read(reg1:Region, reg2:Region, read:pysam.AlignedSegment):
    if read.reference_start > reg1.start and read.reference_start < reg1.end:
        if read.reference_end > reg2.end: return True
    return False
def check_right_clip_read(reg1:Region, reg2:Region, read:pysam.AlignedSegment):
    if read.reference_end > reg2.start and read.reference_end < reg2.end:
        if read.reference_start < reg2.start: return True
    return False

def is_clip_pair(reg1:Region, reg2:Region, bam_in):
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    if reg1.chr_id != reg2.chr_id:
        return False
    max_connect_len = 10000      # 20000，clip之间的删除
    if reg2.start - reg1.end > max_connect_len: # 能够跨越的删除一般不会超过这个距离
        return False
    MIN_CLIP_LEN = 200
    left_support = 0
    right_support = 0
    for read in bam_reader.fetch(reg1.chr_id, reg1.start, reg2.end):
        if read.is_unmapped or read.is_secondary:continue
        cigar = read.cigarstring
        tokens = re.findall("[\d]{0,}[A-Z]{1}", cigar)
        left, right = tokens[0], tokens[-1]
        if left[-1] in "HS" and int(left[:-1]) > MIN_CLIP_LEN:
            if check_left_clip_read(reg1, reg2, read):
                left_support += 1
        if right[-1] in "HS" and int(right[:-1]) > MIN_CLIP_LEN:
            if (reg1, reg2, read):
                right_sucheck_right_clip_readpport += 1
        if left_support > 5 and right_support > 5: return True  # 满足合并条件，返回
    return False

def merge_by_clip_dp(clipinfo_ls, bam_in):
    ''''''
    print("Start cluster by clip or depth")
    pre_rec = None
    merged_clipinfo_ls = []
    clipinfo_ls = sorted(clipinfo_ls, key=lambda x: x[1])
    need_to_merge = []
    for rec in clipinfo_ls:
        if pre_rec:
            if is_clip_pair(Region(pre_rec[0], pre_rec[1], pre_rec[2]), Region(rec[0], rec[1], rec[2]), bam_in): # 引入新的clip region合并函数
                print("{}:{}-{} is_clip_pair".format(pre_rec[0], pre_rec[1], rec[2]))
                need_to_merge.append(rec)
            elif is_del_pair(Region(pre_rec[0], pre_rec[1], pre_rec[2]), Region(rec[0], rec[1], rec[2]), bam_in):
                print("{}:{}-{} is_del_pair".format(pre_rec[0], pre_rec[1], rec[2]))
                need_to_merge.append(rec)
            else:
                if len(need_to_merge) > 1:  # 大于1说明是合并
                    print("merge: ", need_to_merge, "->", [need_to_merge[0][0], need_to_merge[0][1], need_to_merge[-1][2]])
                    # merged_clipinfo_ls.append([need_to_merge[0][0], need_to_merge[0][1], need_to_merge[-1][2], "merged_clip_reg"])
                    merged_clipinfo_ls.append([need_to_merge[0][0], need_to_merge[0][1], need_to_merge[-1][2], "merged_reg"])
                elif len(need_to_merge) == 1:
                    # merged_clipinfo_ls.append([need_to_merge[0][0], need_to_merge[0][1], need_to_merge[0][2], "other_clip_reg"])
                    merged_clipinfo_ls.append([need_to_merge[0][0], need_to_merge[0][1], need_to_merge[-1][2], "other_reg"])
                else:
                    raise(ValueError)
                need_to_merge = [rec]
            pre_rec = rec
        else:
            pre_rec = rec
            need_to_merge.append(rec)
    if pre_rec:
        # if not merge_flag:  # 否则说明上一个是进行过合并的
        #     merged_clipinfo_ls.append([pre_rec[0], pre_rec[1], pre_rec[2], "other_clip_reg"])
        if len(need_to_merge) > 1:  # 大于1说明是合并
            print("merge: ", need_to_merge, "->", [need_to_merge[0][0], need_to_merge[0][1], need_to_merge[-1][2]])
            # merged_clipinfo_ls.append([need_to_merge[0][0], need_to_merge[0][1], need_to_merge[-1][2], "merged_clip_reg"])
            merged_clipinfo_ls.append([need_to_merge[0][0], need_to_merge[0][1], need_to_merge[-1][2], "merged_reg"])
        elif len(need_to_merge) == 1:
            # merged_clipinfo_ls.append([need_to_merge[0][0], need_to_merge[0][1], need_to_merge[0][2], "other_clip_reg"])
            merged_clipinfo_ls.append([need_to_merge[0][0], need_to_merge[0][1], need_to_merge[-1][2], "other_reg"])
        else:
            raise(ValueError)

#########
def contain_in_ls(reg, ctg_reg_ls):
    for REG in ctg_reg_ls:
        if REG[1] > reg[2]: break
        if reg[1] >= REG[1] and reg[2] <= REG[2]: return True
    return False
def main():
    import time
    bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi/my_pipe/step1_mapping/aln.sorted.bam"
    # bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi_ctgs/NC_19/my_pipe/step1_mapping/aln.sorted.bam"
    reg1 = Region("NC_000019.10", 38768500, 38770000)
    reg2 = Region("NC_000019.10", 38774001, 38775500)
    # bed_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi_ctgs/NC_19/my_pipe/step2_candidate_regions/candidate.bed"
    bed_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi/my_pipe/step2_candidate_regions/candidate.bed"
    check_bed = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi/my_pipe/step3_SV_consensus/candidate_op/candidate_op.bed"
    reg_ls = []
    pair_ls = []
    with open(bed_in, "r") as f:
        for line in f:
            fields = line.split()
            if fields[3].startswith("low_dep"):continue
            ctg, start, end  = fields[:3]
            start, end = int(start), int(end)
            reg_ls.append(Region(ctg, start, end))
    # print(is_clip_pair(reg1, reg2, bam_in))
    t1 = time.time()
    dic = defaultdict(list)
    with open(check_bed, "r") as f:
        for line in f:
            fields = line.split()
            if fields[3].startswith("low_dep"):continue
            ctg, start, end  = fields[:3]
            start, end = int(start), int(end)
            dic[ctg].append([ctg, start, end])
    pre_reg = None
    for reg in reg_ls:
        if pre_reg: 
            if is_clip_pair(pre_reg, reg, bam_in):
                print("{}:{}-{} is_clip_pair".format(pre_reg.chr_id, pre_reg.start, reg.end))
                # print(pre_reg, reg, "is_clip_pair")
                pair_ls.append([pre_reg.chr_id, pre_reg.start, reg.end])
                if not contain_in_ls([pre_reg.chr_id, pre_reg.start, reg.end], dic[reg.chr_id]):
                    print("{}:{}-{} not contain in ls".format(pre_reg.chr_id, pre_reg.start, reg.end))
            pre_reg = reg
        else:
            pre_reg = reg
    
    # if is_clip_pair(reg1, reg2, bam_in):
    #     print(reg1, reg2, "is_clip_pair")
    # for reg in pair_ls:
        
    print(time.time() -t1)
if __name__ == "__main__":
    main()