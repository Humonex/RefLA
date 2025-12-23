import pysam
import logging
from collections import namedtuple
import re
# from create_consensus_by_bed import merge_clip_regions
logger = logging.getLogger()
class Record():
    def __init__(self, chr_id, start, end) -> None:
        self.chr_id = chr_id
        self.start = start
        self.end = end
        self.operation = None
        self.info = None
        self.patch_id = None
    def add_info(self, info):
        self.info = info
    def add_operation(self, operation):
        self.operation = operation
    def add_patch_id(self, patch_id):   # 用于填充的asm_id read_id
        self.patch_id = patch_id
    
# Record = namedtuple('Record', ["chr_id", "start", "end", "operation", "info"])
Region = namedtuple('Region', ["chr_id", "start", "end"])
Chr_info = namedtuple('Chr_info', ["chr_id", "chr_len"])


################################### merge_clip_regions #############################################
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
    if dep1 < MIN_DP:
        return True
    else:
        if dep0 > MIN_DP and dep0 < 100 and dep1 / dep0 < 0.5:
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
            if check_right_clip_read(reg1, reg2, read):
                right_support += 1
        if left_support > 5 and right_support > 5: return True  # 满足合并条件，返回
    return False
#######################################################################################






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
def rec_cluster_by_dis(reg_ls_in, dis): # 根据距离进行聚类，注意前面区间的end可能比后面区间要大
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
                reg_ls_out.append([chr_id, reg_start, reg_end, "clip_reg"])
                reg_start = reg[1]
                reg_end = max(reg[2], reg_end)
                need_to_cluster = [reg]
        else:
            reg_start = reg[1]
            reg_end = reg[2]
            need_to_cluster = [reg]
    if reg_start > -1:
        reg_ls_out.append([chr_id, reg_start, reg_end, "clip_reg"])   # 
    return reg_ls_out



def get_candidate_op(bedinfo_ls_in, bam_in):  ## 返回Record ls
    '''处理通过特征得到的候选区间，进行裁剪合并，
    并筛选出哪些可以通过读数解决，哪些需要局部组装来解决
    '''
    if len(bedinfo_ls_in) < 1:
        return []
    ##
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    chr_id = bedinfo_ls_in[0][0]    #
    chr_len = bam_reader.get_reference_length(chr_id)   #

    ##
    bedinfo_ls_in = clip_cov_reg_merge(bedinfo_ls_in)   # 
    candidate_op_ls_out = []
    asm_candidate_ls = []
    ## collect clip reg and low depth reg
    clipinfo_ls = []
    for rec in bedinfo_ls_in:
        # if rec[3] == 'clip_reg':
        #     clipinfo_ls.append(rec)
        # else:   # 
        #     asm_candidate_ls.append([rec[0], rec[1], rec[2], "lowdep_asm"])
        clipinfo_ls.append(rec)    # 全部都丢到一块进行聚类等处理吧，懒得搞了
    # print("clipinfo_ls:", clipinfo_ls)

    ##################################################################################
    ## 对clip_reg的作merge，1、just for del   ###容易出错的一步，将del的两侧合并起来
    ## 2、新增另一种判断模式，只要有一种判真，则合并区间
    ## 但是有很大可能是三个要合并，见鬼了。先不管了
    logger.info("Step: clip merge !!!")
    clipinfo_ls = rec_cluster_by_dis(clipinfo_ls, 500)  # 对合并的数据进行cluster
    # print("clustered_clipinfo_ls:", clipinfo_ls)
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
    # print("\nStep1_out:",merged_clipinfo_ls)



    ### 判断一下 哪些clip_reg 是否可以使用读数进行填充
    reads_span_candidate_op_ls = []    # type: Record
    for rec in merged_clipinfo_ls:
        # res = judge_span_by_reads(rec[0], rec[1], rec[2], bam_in)     # 方法一
        can_solve, solve_seq, read_id = solve_region_by_reads.solve_by_reads(bam_in, rec[0], rec[1], rec[2])
        if can_solve: # 需要判断是否返回的False
            record = Record(rec[0], rec[1], rec[2])
            record.add_info(solve_seq)
            record.add_operation("span_reads")  # merged_clip_reg / other_clip_reg
            record.add_patch_id(read_id)
            reads_span_candidate_op_ls.append(record)
        else:   # 不能横跨
            asm_candidate_ls.append([rec[0], rec[1], rec[2], "clip_asm"])   # 先不管了
    # print(merged_clipinfo_ls)
    # print("\nStep2:",asm_candidate_ls)



    ## 对asm_candidate 进行聚类 merge    //考虑将low_asm也一起聚类
    '''有一个问题：聚类后，原来位于二者之间能够span的区域仍会保留，所以在之后还要筛选一遍
    '''
    start = -1
    end = -1
    cluster_dis = 20000     # asm_candidates 聚类尺度 # 待调整
    asm_candidate_ls = sorted(asm_candidate_ls, key=lambda x: x[1]) # 凡聚类，先排序
    for rec in asm_candidate_ls:
        if start > -1:
            if rec[1] - end < cluster_dis:
                end = max(end, rec[2])   # 注意是选最大end; 原因在于clip 和lowdep聚类分别找到的区间并没有真正的一前一后，后面区间的end可能比前面区间的end 要小
                type_set.add(rec[3])
            else:
                op_record = Record(rec[0], start, end)
                # op_record.add_info(".")
                if "lowdep_asm" in type_set:
                    op_record.add_operation("lowdep_asm")
                else:
                    op_record.add_operation(type_set.pop())
                op_record.add_info(".")
                op_record.add_patch_id(".")
                candidate_op_ls_out.append(op_record)
                start = rec[1]
                end = rec[2]
                type_set = set()
                type_set.add(rec[3])
        else:
            start = rec[1]
            end = rec[2]
            type_set = set()
            type_set.add(rec[3])
    if start > -1:
        op_record = Record(rec[0], start, end)
        if "lowdep_asm" in type_set:
            op_record.add_operation("lowdep_asm")
        else:
            op_record.add_operation(type_set.pop())
        op_record.add_info(".")
        op_record.add_patch_id(".")
        candidate_op_ls_out.append(op_record)

        
    candidate_op_ls_out.extend(reads_span_candidate_op_ls)  ## 
    candidate_op_ls_out.sort(key=lambda rec:rec.start)

    ############################ 根据已有的 asm 区间对聚类 ##########################????
    ''' 根据asm 区域的大小等信息，使用asm区域合并周围小的asm区域，并将asm区域进行往两侧进行一定程度的延伸'''
    ls_in = []
    ls_in.sort(key=lambda rec:rec.start)
    ls_out = []
    radius = 20000      # 搜索周围20000bp的区域 for hifi data
    target_ls = []    # 找到目标区域（asm 区域，合并侧翼的）    [[rec_idx, rec]....]
    for idx, rec in enumerate(ls_in):
        op = rec[3]
        if op.endswith("asm"):
            target_ls.append([idx, rec])
            
    

    for idx, target_rec in target_ls:
        start, end = target_rec.start, target_rec.end
        merge_ls = [target_rec]
        left = target_rec.start - radius if target_rec.start - radius > 0 else 0
        right = target_rec.end + radius if target_rec.end + radius < chr_len else chr_len
        for i in range(idx - 1, -1, -1):    # for j in range(i-1, -1, -1):
            if ls_in[i].end >= left:   # 合并
                merge_ls.append(ls_in[i])
                left = ls_in[i].start - radius if ls_in[i].start - radius > 0 else 0
            else:break

        for i in range(idx + 1, len(ls_in)):   # for j in range(i+1, len(bedinfo_ls)):
            if ls_in[i].start <= right:
                merge_ls.append(ls_in[i])
                right = ls_in[i].end + radius if ls_in[i].end + radius < chr_len else chr_len
            else:break

        ##  对merge_ls
        merge_ls.sort(key=lambda rec:rec.start)

    for rec in ls_in:   # 从前往后找asm
        op = rec[3]     # 
        if op.endswith("asm"):




























    #############################
    '''消除重叠的一些区间，主要是由于上面的合并会导致合并得到一些新的区间可能会将旧区间掩盖（特别是一些reads_span区间）'''
    final_candidate_op_ls = []
    start = -1
    end = -1
    need_to_merge = []
    print("candidate_op_ls_out:")
    for rec in candidate_op_ls_out: print(rec.chr_id, rec.start, rec.end, rec.operation)

    for rec in candidate_op_ls_out:     # 消除交叠的区间
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

    


    print("final_candidate_op_ls:")
    for rec in final_candidate_op_ls: print(rec.chr_id, rec.start, rec.end, rec.operation)
    return final_candidate_op_ls

if __name__ == "__main__":
    ## Debug
    pass