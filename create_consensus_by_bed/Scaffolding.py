

from collections import defaultdict, namedtuple, Counter
from Bio import Seq
import pysam
import logging
from create_consensus_by_bed.Utils import Record, Connect_info
from create_consensus_by_bed import fasta_parser
# from Utils import Record, Connect_info
# import fasta_parser
logger = logging.getLogger()
Region = namedtuple('Region', ["chr_id", "start", "end"])
Chr_info = namedtuple('Chr_info', ["chr_id", "chr_len"])

def get_rec_by_asm(ctg, bam_in, candidate_op_ls, asm_fa_dic, config): # 传入asm_to_ref.bam, op_Ls
    '''
    process one ctg, 
    '''
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    final_op_ls = []
    asm_candidate_ls = []
    min_MQ = 20
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
            # asm_candidate_ls.append(rec)
            # logger.info("{}:{}-{} asm".format(rec.chr_id, rec.start, rec.end))

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
                # logger.info("{} replace_with_denovo: {}".format(rec.chr_id, patch_ls))
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

    # min_MQ = 20
    # ##  for asm_candidate_ls
    # for rec in asm_candidate_ls:
    #     ## 判断是否几乎whole contig 被认为是difficult region
    #     # ctg_len = ctg_len_dic[rec.chr_id]
    #     ctg_len = bam_reader.get_reference_length(rec.chr_id)
    #     if ctg_len == 0:
    #         print(rec.chr_id)
    #         exit(0)
    #     '''solve difficult contig'''
    #     if (rec.end - rec.start) / ctg_len > 0.8:    ## too difficult contig???一条染色体有一半是复杂区域，我们全部进行从头组装。使用从头组装结果来替代。
    #         info_ls = []
    #         patch_ls = []
    #         for read in bam_reader.fetch(rec.chr_id):
    #             if read.is_secondary or read.is_supplementary or read.mapping_quality < min_MQ:
    #                 continue
    #             info_ls.append(read.query_sequence)     # 收集所有的primary，将序列记录下来
    #             patch_ls.append(read.query_name)
    #         new_rec = Record(rec.chr_id, 0, ctg_len)
    #         new_info = ",".join(info_ls)    # 将局部组装的记录下来
    #         patch_ids = ",".join(patch_ls)
    #         for info in info_ls:
    #             new_info += str(info) + ","     # 将局部组装的记录下来
    #         ## 两类处理方式     由于部分染色体质量过差，选择全部从头组装。阈值设为0.7
    #         # new_rec.add_info(new_info)
    #         # new_rec.add_patch_id(patch_ids)     #
    #         new_rec.add_operation("replace_with_denovo")
    #         new_rec.add_info(".")
    #         new_rec.add_patch_id(".")   # 全部丢给从头组装
    #         print("{} replace_with_denovo: {}".format(rec.chr_id, patch_ls))
    #         # logger.info("{} replace_with_denovo: {}".format(rec.chr_id, patch_ls))
    #         final_op_ls = [new_rec]
    #         return ctg, final_op_ls
        
    #     '''普通情形'''
    #     ## 找到pass左边和pass右边的最佳序列
    #     candidate_dic = get_candidate_ctgs(bam_in, Region(rec.chr_id, rec.start, rec.end), config)
    #     left_candidate_ls = candidate_dic["left"]
    #     right_candidate_ls = candidate_dic["right"]     # 
    #     ## 
        
    #     if rec.start < 1000: # 左端  取pass右端
    #         # print("{}:{}-{} left telomere".format(rec.chr_id, rec.start, rec.end))
    #         new_rec = solve_t0(right_candidate_ls, rec, asm_fa_dic)
    #     elif rec.end > ctg_len - 1000: # 右端，取pass左端
    #         # print("{}:{}-{} right telomere".format(rec.chr_id, rec.start, rec.end))
    #         new_rec = solve_t1(left_candidate_ls, rec, asm_fa_dic)
    #     else:
    #         # print("{}:{}-{} middle type".format(rec.chr_id, rec.start, rec.end))
    #         new_rec = solve_t2(left_candidate_ls, right_candidate_ls, rec, asm_fa_dic)
    #     final_op_ls.append(new_rec)

    return ctg, final_op_ls

def convert_reference_pos_to_raw_pos(read, candidate_pos, include_hardclip:bool):  # 输入read信息和要求的参考上的位置，将参考的位置转换为read上的位置
    '''将参考上某个坐标转换到原始读数上的坐标，所以需要将hardclip计算进去'''
    # candidate_pos = set(candidate_pos)
    # raw_ref_pos_map={}
    if candidate_pos < read.reference_start or candidate_pos > read.reference_end:
        print("Error:", candidate_pos)
        raise ValueError
        # pass
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

def solve_t0(candidate_ls, rec, asm_fa_dic):    # 染色体左端点的类型
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

class UnionFind:
    def __init__(self):
        self.root = {}

    def find(self, x):
        # 如果 x 不在映射中，那么它是最新的ID，直接返回 x
        if x not in self.root:
            return x
        # 路径压缩：更新沿路径上的每个节点，使其指向最终的根节点
        if self.root[x] != x:
            self.root[x] = self.find(self.root[x])
        return self.root[x]

    def union(self, x, y, new_id):
        '''
        x:ctg_id1
        y:ctg_id2
        new_id:new_scaff_id
        '''
        rootX = self.find(x)
        rootY = self.find(y)
        # 将 rootY 和 rootX 的合并后的新ID存储起来
        self.root[rootX] = self.root[rootY] = new_id
        # 返回新生成的ID

    def connected(self, x, y):
        return self.find(x) == self.find(y)


class Bridge:
    def __init__(self, bridge_id, bridge_seq) -> None:
        self.bridge_id = bridge_id
        self.bridge_seq = bridge_seq
        self.bridge_length = len(bridge_seq)
        self.connect_info_ls = []
    def get_reverse_pos(self, pos):
        return self.bridge_length - pos
    def get_reverse_seq(self):
        return Seq.reverse_complement(self.bridge_seq)
    def add_connect(self, contig_id, bound, bound_pos, bridge_start, bridge_end, strand, align_info:pysam.AlignedSegment):
        '''注意要加上clip长度
        依据参考来，如果是反向比对，则是反向互补链的位置
        bridge_start: bridge 开始比对的位置
        bridge_end: 
        '''
        try:
            bridge_bound_pos = convert_reference_pos_to_raw_pos(align_info, bound_pos, True)    # 反向的比对则是反向互补链中的位置
        except:
            print(contig_id, bound, bound_pos, bridge_start, bridge_end, strand, align_info.reference_start, align_info.reference_end)
            raise ValueError
        if strand == "-": 
            # print("Reverse:", bridge_start, bridge_end, bridge_bound_pos, ", length:", self.bridge_length)
            # bridge_start, bridge_end = self.bridge_length - bridge_end, self.bridge_length - bridge_start
            # bridge_bound_pos = self.bridge_length - bridge_bound_pos
            # print("to:", bridge_start, bridge_end, bridge_bound_pos)
            raw_bridge_start = self.bridge_length - bridge_end
            raw_bridge_end = self.bridge_length - bridge_start
            raw_bridge_bound_pos = self.bridge_length - bridge_bound_pos
        else:
            raw_bridge_start, raw_bridge_end, raw_bridge_bound_pos = bridge_start, bridge_end, bridge_bound_pos
        self.connect_info_ls.append({
        "contig_id":contig_id, \
        "bound":bound, \
        "bound_pos":bound_pos, \
        "ref_id":align_info.reference_name, \
        "ref_start": align_info.reference_start, \
        "ref_end": align_info.reference_end, \
        "bridge_start":bridge_start, \
        "bridge_end":bridge_end, \
        "bridge_bound_pos": bridge_bound_pos, \
        "raw_bridge_start": raw_bridge_start, \
        "raw_bridge_end": raw_bridge_end, \
        "raw_bridge_bound_pos": raw_bridge_bound_pos, \
        "overlap":bridge_end - bridge_start, \
        "strand":strand, \
        "valid":True, \
        # "patch_seq":patch_seq, \
        "align_info":align_info})
        
        '''align_info = Align_info(ref_id, ref_start, ref_end, self.bridge_id, bridge_start, bridge_end, strand)
        contig_id = ref_id + ":" + str(ref_start) + "-" + str(ref_end)
        self.connect_info_ls.append([contig_id, bound, align_info])'''
    def __str__(self):
        connect_info_strs = []
        for info in self.connect_info_ls:
            align_info_repr = repr(info["align_info"])  # You might need to adjust this depending on the representation of align_info
            connect_info_strs.append(
                f"{{'contig_id': {info['contig_id']}, 'bound': {info['bound']}, 'bound_pos':{info['bound_pos']}, 'ref_start': {info['ref_start']}, 'ref_end': {info['ref_end']}, "
                f"'bridge_start': {info['bridge_start']}, 'bridge_end': {info['bridge_end']}, 'bridge_bound_pos': {info['bridge_bound_pos']}, "
                f"'raw_bridge_start': {info['raw_bridge_start']}, 'raw_bridge_end': {info['raw_bridge_end']}, 'raw_bridge_bound_pos': {info['raw_bridge_bound_pos']}, "
                f"'overlap': {info['overlap']}, 'strand': {info['strand']}, "
                f"'align_info': {align_info_repr}}}"
            )
        connect_infos = "[" + ", ".join(connect_info_strs) + "]"
        return (f"Bridge(bridge_id={self.bridge_id}, bridge_seq='.', bridge_length={self.bridge_length}, "
                f"connect_info_ls={connect_infos})")

    def __repr__(self):
        return (f"Bridge(bridge_id={self.bridge_id}, bridge_seq='.', "
                f"connect_info_ls={len(self.connect_info_ls)} items)")

reverse_strand_dic = {"+":"-", "-":"+"}
class Contig:
    def __init__(self, contig_id, seq) -> None:
        self.id = contig_id
        self.left_bound_pos = None
        self.right_bound_pos = None
        self.seq = seq
        self.gap_ls = []    # 使用组装结果进行填充的gap //目前不加上读数填充
        self.strand = "+"   # 记录当前contig的方向，初始为 "+"。不用于记录scaff的方向，scaff的方向没有意义，因为由多个组成
        self.length = len(seq)
        self.span = None
        # 
        self.left_bridge_id = None
        self.left_align_pos = None
        self.left_strand = None
        self.left_align_len = None
        self.right_bridge_id = None
        self.right_align_pos = None
        self.right_strand = None
        self.right_align_len = None
        # 
        self.died = False   # 表示当前contig是否已经死掉，丢弃
        self.l_solved = False
        self.r_solved = False
        self.l_patch_len = 0
        self.r_patch_len = 0
        # 
        self.is_scaffold = False
        # 用于存储桥接的contig_id信息，在对scaffolder进行操作的时候，方便记录原始片段的相关变化
        self.contig_id_ls = [contig_id]
        self.bridge_id_ls = []
        # 用于
        self.repre_loca = None
    def get_reverse_seq(self):
        return Seq.reverse_complement(self.seq)
    def add_bound_pos(self, start, end):
        if start >= end:
            raise ValueError
        self.span = end - start
        self.left_bound_pos, self.right_bound_pos = start, end
    def add_gap_ls(self, gap_ls):
        self.gap_ls = gap_ls
    def add_l_connect(self, bridge_id, strand, pos, align_len):
        self.left_bridge_id = bridge_id
        self.left_align_pos = pos
        self.left_strand = strand
        self.left_align_len = align_len

    def add_r_connect(self, bridge_id, strand, pos, align_len):
        self.right_bridge_id = bridge_id
        self.right_align_pos = pos
        self.right_strand = strand
        self.right_align_len = align_len
    def get_bound_pos(self, bound):
        if bound == "l":
            return self.left_bound_pos
        elif bound == "r":
            return self.right_bound_pos
        else:
            raise ValueError
    def solve_bound(self, bound, bridge_id, seq):
        # if self.id == 'CM008973.1:716793-794042':
        #     exit(0)
        if bound == "l" and not self.l_solved:
            self.seq = seq + self.seq
            self.l_solved = True
            self.l_patch_len = len(seq)
            self.bridge_id_ls = [bridge_id] + self.bridge_id_ls
        elif bound == "r" and not self.r_solved:
            self.seq = self.seq + seq
            self.r_solved = True
            self.r_patch_len = len(seq)
            self.bridge_id_ls = self.bridge_id_ls + [bridge_id]
        else:
            print(bound, self.l_solved,self.r_solved, bridge_id)
            raise ValueError
    def kill(self):
        self.died = True
    def cal_loca_bound(self):
        ctg_id_ls = self.contig_id_ls
        loca = None
        for ctg_id in ctg_id_ls:
            start, end = ctg_id.split(":")[1].split("-")
            start, end = int(start), int(end)
            if loca is None:
                loca = start
            else:
                loca = min(start, loca)
            if self.left_bound_pos is None:
                self.left_bound_pos = start
                self.right_bound_pos = end
            else:
                self.left_bound_pos = min(self.left_bound_pos, start)
                self.right_bound_pos = max(self.right_bound_pos, end)
        self.repre_loca = loca
    @staticmethod
    def cal_reversed_portion(contig, ctg_dic):
        whole_span = 0
        revered_span = 0
        whole_num = 0
        revered_num = 0
        for ctg_id in contig.contig_id_ls:
            ctg = ctg_dic[ctg_id]
            whole_span += ctg.span
            whole_num += 1
            if ctg.strand == "-":
                revered_span += ctg.span
                revered_num += 1
        rev_len_por = revered_span / whole_span
        rev_num_por = revered_num / whole_num
        return rev_len_por, rev_num_por
    @staticmethod
    def sort_by_ref(contig_ls):
        
        pass
    @staticmethod
    def generate_id(chr_id, start, end):
        return chr_id + ":" + str(start) + "-" + str(end)
    @staticmethod
    def write_contig(file_out, contig_dic:dict):
        new_dic = {}
        for contig in contig_dic.values():
            if not contig.died:
                print(contig.id, contig.l_solved, contig.r_solved)
                new_dic[contig.id] = contig.seq
        fasta_parser.write_fasta_dict(new_dic, file_out)
        return new_dic
    @staticmethod
    def scaffolding(contig1, contig2, bridge_id, new_seq):
        ''''''
        if contig1.died or contig2.died: raise ValueError
        contig1.kill()
        contig2.kill()
        contig_id = contig1.id + "_" + contig2.id
        contig = Contig(contig_id, new_seq)
        contig.contig_id_ls = contig1.contig_id_ls + contig2.contig_id_ls
        contig.bridge_id_ls = contig1.bridge_id_ls + [bridge_id] +  contig2.bridge_id_ls
        contig.is_scaffold = True
        contig.l_solved = contig1.l_solved
        contig.r_solved = contig2.r_solved
        contig.l_patch_len, contig.r_patch_len = contig1.l_patch_len, contig2.r_patch_len
        return contig
    @staticmethod
    def reverse_scaffold(contig, contig_dic):
        # contig = Contig()
        contig.contig_id_ls.reverse()
        contig.bridge_id_ls.reverse()
        contig.seq = Seq.reverse_complement(contig.seq)
        for contig_id in contig.contig_id_ls:    # reverse block_contig strand
            contig_dic[contig_id].strand = reverse_strand_dic[contig_dic[contig_id].strand]
        contig.l_solved, contig.r_solved = contig.r_solved, contig.l_solved
        contig.l_patch_len, contig.r_patch_len = contig.r_patch_len, contig.l_patch_len
        contig.id = "_".join(contig.contig_id_ls)
    # @staticmethod
    # def rever_id(id):
    #     contig_id_ls = id.split("_")
    #     new_id = ""
    #     for contig_id in contig_id_ls[::-1]:
    #         chr_id, reg = contig_id.split(":")
    #         start, end = reg.split("-")
    #         new_contig_id = chr_id + ":" + end + "-" + start
    #         if new_id == "":
    #             new_id = new_contig_id
    #         else:
    #             new_id += "_" + new_contig_id
    #     return id
class Scaffold:
    def __init__(self, ):
        pass


def get_best_bridge(candidate_ls):
    if len(candidate_ls) == 0: return None
    best_ctg = None
    for ctg in candidate_ls:    # 对候选ctg作过滤，加上对clip的限制，和长度的限制
        if best_ctg != None:
            if ctg.query_alignment_length > best_ctg.query_alignment_length:
                best_ctg = ctg
        else:
            best_ctg = ctg
    return best_ctg

def get_patch_seq(align_info:pysam.AlignedSegment, bridge:Bridge, bound, target_pos): # 
    '''
    target_pos为参考上跨过边界的点
    direction表示与参考良好match的那部分是左侧还是右侧
    '''
    bridge_seq = bridge.bridge_seq
    ##
    if align_info.is_reverse:   # 判断是否反向互补链
        bridge_seq = Seq.reverse_complement(bridge_seq)
    ##
    if bound == "r":    # 
        query_pos = convert_reference_pos_to_raw_pos(align_info, target_pos, True)
        # if align_info.is_reverse:   # 判断是否反向互补链
        #     query_pos = bridge.get_reverse_pos(query_pos)
        patch_seq = bridge_seq[query_pos:]
    elif bound == "l":  # contig left bound
        query_pos = convert_reference_pos_to_raw_pos(align_info, target_pos, True)
        # if align_info.is_reverse:   # 判断是否反向互补链
        #     query_pos = bridge.get_reverse_pos(query_pos)
        patch_seq = bridge_seq[:query_pos]
    else:
        # raise("ERROR bound of {}:{}-{}".format(rec.chr_id, rec.start, rec.end))
        raise("ERROR bound")
    return patch_seq

def get_raw_pos(read:pysam.AlignedSegment):
    '''(operation, length)'''
    '''对于反向比对序列，位置是反向互补之后的序列中的位置'''
    cigar = read.cigartuples
    left_clip, right_clip = 0, 0
    if cigar[0][0] == 5:
        left_clip = cigar[0][1]
    return read.query_alignment_start + left_clip, read.query_alignment_end + left_clip
    # if cigar[-1][0] == 5:
    #     right_clip = cigar[-1][1]
    # if read.is_reverse: # 反向比对
    #     return read.query_alignment_start + right_clip, read.query_alignment_end + right_clip
    # else:
    #     return read.query_alignment_start + left_clip, read.query_alignment_end + left_clip

def check_bound_strand(bound1, strand1, bound2, strand2):
    if (bound1 != bound2 and strand1 == strand2) or (bound1 == bound2 and strand1 != strand2):
        return True
    return False

def scaffold_2_contig(new_contig_dic, Old2new, bridge, connect_info1, connect_info2, failed_scaff_set, ):
    '''Process the two connect types'''
    '''
    ++, rl
    --, rl
    +-, rr
    +-, ll
    '''
    bridge_id = bridge.bridge_id
    # connect_info_ls = bridge.connect_info_ls
    # connect_info1 = connect_info_ls[0]
    # connect_info2 = connect_info_ls[1]
    bound1 = connect_info1["bound"]
    strand1 = connect_info1["strand"]
    bound2 = connect_info2["bound"]
    strand2 = connect_info2["strand"]
    if connect_info1['contig_id'] == connect_info2['contig_id']:
        failed_scaff_set.add(bridge_id)
        print("Failed scaff for {}, may become a loop:{}".format(bridge_id, bridge))
        return
    if bound1 != bound2 and strand1 == strand2:
        ''' same strand, do not reverse the contig strand
        1、keep/reverse the bridge seq
        2、put the right bound on the first
        '''
        # 1、
        if strand1 == "-":
            bridge_seq = bridge.get_reverse_seq()
        else:
            bridge_seq = bridge.bridge_seq
        # 2、swap the block, 确保第一个是right bound用于连接
        if bound1 == "l":
            print("Reverse the block:{}".format(connect_info1, connect_info2))
            connect_info1, connect_info2 = connect_info2, connect_info1 # 注意只是交换引用
        # 3、get bridge_patch_seq
        bridge_bound_pos1, bridge_bound_pos2 = connect_info1['bridge_bound_pos'], connect_info2['bridge_bound_pos']
        if bridge_bound_pos1 > bridge_bound_pos2:
            failed_scaff_set.add(bridge_id)
            print("Failed scaff for {}, error bound_pos, type1:{}".format(bridge_id, bridge))
            return
        else:
            bridge_patch_seq = bridge_seq[bridge_bound_pos1:bridge_bound_pos2]
            # 4、判断当前实际桥接的scaff中，处于桥接位置的contig的strand
            target_strand1 = "+"    # 确保桥接contig为+
            target_strand2 = "+"
            contig_id1 = connect_info1['contig_id'] # 桥接contig
            contig_id2 = connect_info2['contig_id']
            scaff_id1 = Old2new.find(contig_id1)
            scaff_id2 = Old2new.find(contig_id2)
            if scaff_id1 == scaff_id2:
                failed_scaff_set.add(bridge_id)
                print("Failed scaff for {}, may become a loop:{}".format(bridge_id, bridge))
                return
            if new_contig_dic[contig_id1].strand != target_strand1:
                Contig.reverse_scaffold(new_contig_dic[scaff_id1], new_contig_dic)
            if new_contig_dic[contig_id2].strand != target_strand2:
                Contig.reverse_scaffold(new_contig_dic[scaff_id2], new_contig_dic)
            # check bridge contig id
            if new_contig_dic[scaff_id1].contig_id_ls[-1] != contig_id1 or new_contig_dic[scaff_id2].contig_id_ls[0] != contig_id2:
                print("{},{};Error bridge contig id:{}, {}".format(contig_id1, scaff_id2, new_contig_dic[scaff_id1].contig_id_ls, new_contig_dic[scaff_id2].contig_id_ls))
                raise ValueError
            seq1 = new_contig_dic[scaff_id1].seq
            seq2 = new_contig_dic[scaff_id2].seq
            new_seq = seq1 + bridge_patch_seq + seq2
            new_contig = Contig.scaffolding(new_contig_dic[scaff_id1], new_contig_dic[scaff_id2], bridge_id, new_seq)
            new_contig_dic[new_contig.id] = new_contig
            # 5、更新并查集及其他数据
            connect_info1['valid'] = False
            connect_info2['valid'] = False
            Old2new.union(scaff_id1, scaff_id2, new_contig.id)
            print("Scaff {},{}->{} done!!!".format(scaff_id1, scaff_id2, new_contig.id))
    elif bound1 == bound2 and strand1 != strand2:
        '''由于bridge一部分正向比对到一条seq，一部分反向比对到一条seq，bridge不作修改，对bridge反向比对的桥接contig进行reverse。'''
        '''保证原来contig为+的不动进行修正
        +- ll -> ++ lr -> ++ rl
        +- rr -> ++ rl
        '''
        '''reverse the revered contig'''
        '''
        +-ll
        1、ll:R(ctg1)+B_seq+ctg2
        2、rr:ctg1+B_seq+R(ctg2)
        '''
        # 1、
        bridge_seq = bridge.bridge_seq
        # 2、-> ++
        '''rr +-    ll -+'''
        if (bound2 == 'l' and strand2 == '-') or (bound2 == 'l' and strand2 == '-'):    # block是
            connect_info1, connect_info2 = connect_info2, connect_info1 # 注意只是交换引用
            bound1, bound2 = bound2, bound1
            strand1, strand2 = strand2, strand1
        if bound1 == 'l':   # 
            target_strand1 = "-"    # 确保桥接contig1为+
            target_strand2 = "+"    # 桥接contig2为-，与bridge-比对的target_strand为-。->确保ctg2为l+或者
        else:
            target_strand1 = "+"    # 确保桥接contig1为+
            target_strand2 = "-"    # 桥接contig2为-，与bridge-比对的target_strand为-。->确保ctg2为l+或者
        # 3、get bridge_patch_seq
        bridge_bound_pos1, bridge_bound_pos2 = connect_info1['raw_bridge_bound_pos'], connect_info2['raw_bridge_bound_pos']
        if bridge_bound_pos1 > bridge_bound_pos2:
            failed_scaff_set.add(bridge_id)
            print("Failed scaff for {}, error bound_pos, type2:{}".format(bridge_id, bridge))
        else:
            bridge_patch_seq = bridge_seq[bridge_bound_pos1:bridge_bound_pos2]
            # 4、判断当前实际桥接的scaff中，处于桥接位置的contig的strand
            contig_id1 = connect_info1['contig_id'] # 桥接contig
            contig_id2 = connect_info2['contig_id']
            scaff_id1 = Old2new.find(contig_id1)
            scaff_id2 = Old2new.find(contig_id2)
            if scaff_id1 == scaff_id2:
                failed_scaff_set.add(bridge_id)
                print("Failed scaff for {}, may become a loop:{}".format(bridge_id, bridge))
                return
            if new_contig_dic[contig_id1].strand != target_strand1:
                Contig.reverse_scaffold(new_contig_dic[scaff_id1], new_contig_dic)
            if new_contig_dic[contig_id2].strand != target_strand2:
                Contig.reverse_scaffold(new_contig_dic[scaff_id2], new_contig_dic)
            # check bridge contig id
            if new_contig_dic[scaff_id1].contig_id_ls[-1] != contig_id1 or new_contig_dic[scaff_id2].contig_id_ls[0] != contig_id2:
                print("{},{};Error bridge contig id:{}, {}".format(contig_id1, scaff_id2, new_contig_dic[scaff_id1].contig_id_ls, new_contig_dic[scaff_id2].contig_id_ls))
                raise ValueError
            seq1 = new_contig_dic[scaff_id1].seq
            seq2 = new_contig_dic[scaff_id2].seq
            new_seq = seq1 + bridge_patch_seq + seq2
            new_contig = Contig.scaffolding(new_contig_dic[scaff_id1], new_contig_dic[scaff_id2], bridge_id, new_seq)
            new_contig_dic[new_contig.id] = new_contig
            # 5、更新并查集及其他数据
            connect_info1['valid'] = False
            connect_info2['valid'] = False
            Old2new.union(scaff_id1, scaff_id2, new_contig.id)
            print("Scaff {},{}->{} done!!!".format(scaff_id1, scaff_id2, new_contig.id))
    else:
        print("Failed scaff for {}, error: bound/strand conflict. {}".format(bridge_id, bridge))
        failed_scaff_set.add(bridge_id)

def judge_free(contig, ctg_dic):
    chrs = set()
    for ctg_id in contig.contig_id_ls:
        chr = ctg_id.split(":")[0]  # chr1:0-1000
        chrs.add(chr)
    if len(chrs) > 1:
        print("Dirrerent chr free:", contig.id)
        return True
    rev_len_por, rev_num_por = Contig.cal_reversed_portion(contig, ctg_dic)
    if rev_len_por > 0.3 or rev_num_por > 0.3:
        print(rev_len_por, rev_num_por)
        print("Reversed free:", contig.id)
        return True
    return False

def scaffold(candidate_scaff_dic, new_contig_dic):
    failed_scaff_set = set()
    # 使用并查集
    Old2new = UnionFind()
    ## round1: process the seq bridge
    for bridge_id, bridge in candidate_scaff_dic.items():
        # bridge = Bridge()
        connect_info_ls = bridge.connect_info_ls
        if len(connect_info_ls) == 2:   # 处理两个的类型
            # continue
            scaffold_2_contig(new_contig_dic, Old2new, bridge, connect_info_ls[0], connect_info_ls[1], failed_scaff_set) # 
        elif len(connect_info_ls) == 3:
            print("3 connect type:", bridge_id, connect_info_ls)
            failed_scaff_set.add(bridge_id)   # 3个的类型，必有一个处理不了
            '''可能其中两个是同一个contig_id'''
            contig_id_counter = Counter([connect_info['contig_id'] for connect_info in connect_info_ls]).most_common()
            print(contig_id_counter[0])
            if contig_id_counter[0][1] > 1: # dup contig_id
                '''dup contig_id
                组装结果横跨某条contig左右两端，判断使用dup_contig的哪端与single连接
                根据桥接的方向进行: (bound1 != bound2 and strand1 == strand2) or (bound1 == bound2 and strand1 != strand2)
                '''
                dup_contig_id = contig_id_counter[0][0]
                not_dup_contig_id = contig_id_counter[1][0]
                #
                bound1 = None
                strand1 = None
                bound2 = None
                strand2 = None
                processed = False
                for connect_info in connect_info_ls:
                    if connect_info['contig_id'] == not_dup_contig_id:
                        connect_info1 = connect_info
                        bound1 = connect_info1['bound']
                        strand1 = connect_info1['strand']
                for connect_info in connect_info_ls:
                    if connect_info['contig_id'] == dup_contig_id:
                        connect_info2 = connect_info
                        bound2 = connect_info2['bound']
                        strand2 = connect_info2['strand']
                        if check_bound_strand(bound1, strand1, bound2, strand2):
                            print("3 connect, chooese:{}, {}".format(connect_info1, connect_info2))
                            scaffold_2_contig(new_contig_dic, Old2new, bridge, connect_info1, connect_info2, failed_scaff_set)
                            processed = True
                            break
                if not processed:
                    print("Not process the 3 connect bridge")
                    failed_scaff_set.add(bridge_id)
            else:
                '''
                no dup contig_id，由于桥接只能两条，所以需要做放弃一个桥接点。
                1、优先桥接same chr
                2、根据桥接序列位置/方向筛选
                3、根据桥接序列overlap筛选
                '''
                processed = False
                chr_id_counter = Counter([connect_info['ref_id'] for connect_info in connect_info_ls]).most_common()
                if chr_id_counter[0][1] >= 2:    # dup chr
                    dup_chr_id = chr_id_counter[0][0]
                    new_conn_ls = []
                    for connect_info in connect_info_ls:
                        if connect_info['ref_id'] == dup_chr_id:
                            new_conn_ls.append(connect_info)
                    new_conn_ls.sort(key=lambda connect_info: connect_info["overlap"], reverse=True)
                    if check_bound_strand(new_conn_ls[0]['bound'], new_conn_ls[0]['strand'], new_conn_ls[1]['bound'], new_conn_ls[1]['strand']):
                        print("3 connect, chooese:{}, {}".format(new_conn_ls[0], new_conn_ls[1]))
                        scaffold_2_contig(new_contig_dic, Old2new, bridge, new_conn_ls[0], new_conn_ls[1], failed_scaff_set)
                        processed = True
                else:
                    connect_info_ls.sort(key=lambda connect_info: connect_info["overlap"], reverse=True)
                    if check_bound_strand(connect_info_ls[0]['bound'], connect_info_ls[0]['strand'], connect_info_ls[1]['bound'], connect_info_ls[1]['strand']):
                        print("3 connect, chooese:{}, {}".format(connect_info_ls[0], connect_info_ls[1]))
                        scaffold_2_contig(new_contig_dic, Old2new, bridge, connect_info_ls[0], connect_info_ls[1], failed_scaff_set)
                        processed = True
                if not processed:
                    print("Not process the 3 connect bridge")
                    failed_scaff_set.add(bridge_id)
        else:   # too more connect, do edge expand
            print("Failed scaff for {}, too more connect".format(bridge_id))
            failed_scaff_set.add(bridge_id)
            '''更复杂的情况，放弃桥接，仅做边界扩展'''
    
    ## round2: process the unused bridge, do edge expand
    '''延展对象：通过connect_info的valid判断这个桥接是否已经被使用
    1、桥接失败
    2、3connect types
    '''
    apply_all = True
    print("failed_scaff_set:", failed_scaff_set)
    for bridge_id in failed_scaff_set:
        bridge = candidate_scaff_dic[bridge_id]
        bridge.connect_info_ls.sort(key=lambda connect_info: connect_info["overlap"], reverse=True)
        connect_info_ls = [connect_info for connect_info in bridge.connect_info_ls if connect_info['valid']]    # check valid
        connect_info_ls.sort(key=lambda connect_info: connect_info["overlap"], reverse=True)
        print("{}, Old:{}, New:{}".format(bridge_id, bridge.connect_info_ls, connect_info_ls))
        apply_conn_ls = []
        if not apply_all:   # only apply on the best align pos
            apply_conn_ls.append(connect_info_ls[0])
        else:
            apply_conn_ls = connect_info_ls
        for connect_info in apply_conn_ls:
            contig_id = connect_info['contig_id']
            bound = connect_info["bound"]
            bound_pos = connect_info["bound_pos"]
            align_info = connect_info["align_info"]
            contig = new_contig_dic[contig_id]
            scaff_id = Old2new.find(contig_id)
            scaff_tig = new_contig_dic[scaff_id]
            # 
            target_strand = "+" # 延展需要先将桥接contig调整至目标strand
            if contig.strand != target_strand:
                Contig.reverse_scaffold(new_contig_dic[scaff_id], new_contig_dic)
            # print("bound:{}, Boundpos:{}, scaff:{}".format(bound, bound_pos, scaff_tig))
            scaff_tig.solve_bound(bound, bridge_id, get_patch_seq(align_info, bridge, bound, bound_pos))
            if bridge_id == '112':
                print(new_contig_dic[scaff_id].l_solved, new_contig_dic[scaff_id].r_solved)
    print("Old2new:", Old2new.root)
    ## round3: judge reverse
    '''根据读数原始片段的方向信息'''
    print("-----------------------judge reverse-----------------------")
    final_dic = defaultdict()
    reversed_portion = 0.5
    for contig in new_contig_dic.values():
        if contig.died: continue
        print(contig.id, contig.contig_id_ls, contig.bridge_id_ls)
        # 
        # print(Contig.cal_reversed_portion(contig, new_contig_dic))
        rev_len_por, rev_num_por = Contig.cal_reversed_portion(contig, new_contig_dic)
        if rev_len_por > reversed_portion or rev_num_por > reversed_portion:
            print("Reverse:{}".format(contig.id))
            Contig.reverse_scaffold(contig, new_contig_dic)
            print("To:{}".format(contig.id))
        # print(Contig.cal_reversed_portion(contig, new_contig_dic))
        final_dic[contig.id] = contig
    ## Get connect info
    chr_ctg_dic = defaultdict(list)
    chr_ctg_dic['freed'] = []
    # 1、put the contig to the right chr
    print("----------------------------Anchored contig to the ref----------------------------")
    for contig in final_dic.values():
        if judge_free(contig, new_contig_dic):
            chr_ctg_dic['freed'].append(contig)
        else:
            chr = contig.id.split(":")[0]  # chr1:0-1000
            chr_ctg_dic[chr].append(contig)
    # 2、sort
    print("----------------------------Sort contig on the ref----------------------------")
    for chr, ctg_ls in chr_ctg_dic.items(): # RuntimeError: dictionary changed size during iteration
        '''注意运行期间需保证dic不更改'''
        id_ls = [contig.id for contig in ctg_ls]
        print("Chr:{}, ctg_ls:{}".format(chr, id_ls))
        if chr == "freed" or len(ctg_ls) == 0:
            continue
        for contig in ctg_ls:
            contig.cal_loca_bound()
            # print(contig.left_bound_pos, contig.right_bound_pos)
        ctg_ls.sort(key=lambda contig: contig.repre_loca)
        new_ctg_ls = []
        pre_contig = ctg_ls[0]
        for contig in ctg_ls[1:]:
            if contig.left_bound_pos < pre_contig.right_bound_pos:
                if contig.length > pre_contig.length:   # remove precontig
                    print("Remove {} out of candidate scaffold_conn".format(pre_contig.id))
                    # new_ctg_ls.append(contig)
                    chr_ctg_dic['freed'].append(pre_contig)
                    pre_contig = contig
                else:   # skip now contig
                    print("Remove {} out of candidate scaffold_conn".format(contig.id))
                    # new_ctg_ls.append(pre_contig)
                    chr_ctg_dic['freed'].append(contig)
            else:
                new_ctg_ls.append(pre_contig)
                pre_contig = contig
        if pre_contig is not None:
            new_ctg_ls.append(pre_contig)
        chr_ctg_dic[chr] = new_ctg_ls
        # print("new_ctg_ls:", new_ctg_ls)
        # id_ls = [contig.id for contig in new_ctg_ls]
        # print("new_ctg_ls:", id_ls)
    print("Finally:")
    for chr, ctg_ls in chr_ctg_dic.items():
        id_ls = [contig.id for contig in ctg_ls]
        print("Chr:{}, new_ctg_ls:{}".format(chr, id_ls))
        # print("new_ctg_ls:", id_ls)
    # Get connectinfo ls
    print("---------------------------Collect connect info---------------------------")
    connect_info_ls = []
    for chr_id, ctg_ls in chr_ctg_dic.items():
        if chr_id == "freed":
            id_ls = [contig.id for contig in ctg_ls]
            print("Chr:{}, freed_ls:{}".format(chr_id, id_ls))
            cnt = 0
            for contig in ctg_ls:
                connect_info = Connect_info("freed_" + str(cnt), [contig.id], []) # 
                connect_info_ls.append(connect_info)
                cnt += 1
        else:
            if len(ctg_ls) == 0:
                raise ValueError
            connect_ls = [contig.id for contig in ctg_ls]
            connect_info = Connect_info(chr_id, connect_ls, []) # 
            pre_contig = ctg_ls[0]
            id_ls = [contig.id for contig in ctg_ls]
            print("Chr:{}, scaff_ls:{}".format(chr_id, id_ls))
            for contig in ctg_ls[1:]:
                # connect_info.gap_ls.append(str(end - start - patch_len) + "_GAP")   # 记录下difference or 直接记录end - start距离   ？？
                # contig = Contig()
                # print(contig.id, contig.left_bound_pos, pre_contig.id, pre_contig.right_bound_pos)
                Gap_diff = contig.left_bound_pos - pre_contig.right_bound_pos - (contig.l_patch_len + pre_contig.r_patch_len)
                print("{} - {}, Gap:{}, len1:{}, len2:{}".format(pre_contig.id, contig.id, Gap_diff, pre_contig.length, contig.length))
                connect_info.gap_ls.append(str(Gap_diff) + "_GAP")   # 
                pre_contig = contig
            connect_info_ls.append(connect_info)
    print("---------------------------Collect connect info done---------------------------")
    return connect_info_ls
    # return patch_ids

def generate_rec_id(rec):
    return "{}:{}-{}".format(rec.chr_id, rec.start, rec.end)

def consensus(bam, ref_dic, rec_dic, asm_fa_dic, fa_out, work_dir, config, exe_scaffold):
    '''
    将ref拆碎成片段，然后建立桥接信息，进行桥接。
    1、asm_patch的位置，如果是中间的单个asm_patch，填充。
    2、否则该位置是接头位置，从该位置生成新的contig。
    '''
    print("--------------------------------------Start scaffolding--------------------------------------")
    # config
    min_MQ = 40 # add to config dic
    min_clip = 1000
    min_align_length = 5000 # config["min_align_length"]
    min_dis =  1000# config["min_dis"]    # 只保留patch超过一定长度的
    patch_ids = set()
    # 
    contig_ls = list(ref_dic.keys())
    # new_contig_ls = []
    bridge_dic = defaultdict(Bridge)    # bridge_id -> Bridge
    bam_reader = pysam.AlignmentFile(bam, "rb", index_filename=bam+".bai")
    new_contig_dic = defaultdict()  # contig_id -> Contig
    # generate raw bridge_dic
    for bridge_id, bridge_seq in asm_fa_dic.items():
        bridge = Bridge(bridge_id, bridge_seq)
        bridge_dic[bridge_id] = bridge
    # 
    for ctg in contig_ls:
        rec_ls = rec_dic[ctg]
        # rec_ls.sort(key=lambda rec:rec.start)
        if len(rec_ls) == 0:    # 说明当前contig没有修正记录
            continue
        elif len(rec_ls) == 1:
            if rec_ls[0].operation == "keep_ref":   # 摆烂式，弃置
                continue
            elif rec_ls[0].operation == "replace_with_denovo":  # 不记录吧，后面通过碎片给记录进去
                continue
        ctg_seq = ref_dic[ctg]
        ctg_len = len(ctg_seq)
        ## -----------------第一轮遍历获取全部切割的染色体-----------------
        connect_rec_ls = []
        now_ctg_ls = []
        pre_end = 0
        new_seq = ""
        new_contig_start, new_contig_end = 0, 0
        rec2ctg = defaultdict()
        pre_recid = None
        gap_ls = []
        for rec in rec_ls:
            start = rec.start
            end = rec.end
            op_ls = rec.operation.split(",")
            info_ls = rec.info.split(",")
            patch_ls = rec.patch_id.split(",")
            new_seq += ctg_seq[pre_end:rec.start]    # 加上上一段rec的末尾到现在rec前面的中间这段序列
            new_contig_end = start
            # print(ctg, start, end, rec.operation)
            if len(op_ls) == 1: # 
                if start == 0:
                    print("left start rec")
                    new_contig_start, new_contig_end = end, end
                    connect_rec_ls.append(rec)  # 接头位置
                    rec_id = generate_rec_id(rec)
                    rec2ctg[rec_id] = [None, None]
                    pre_recid = rec_id
                elif end >= ctg_len - 1:
                    print("Right end rec")
                    if len(new_seq) > 0:
                        new_contig = Contig(Contig.generate_id(ctg, new_contig_start, new_contig_end), new_seq)
                        new_contig.add_bound_pos(new_contig_start, new_contig_end)
                        new_contig.add_gap_ls(gap_ls)
                        new_contig_dic[new_contig.id] = new_contig
                        now_ctg_ls.append(new_contig)
                        new_contig_start, new_contig_end = end, end
                        gap_ls = []
                    if pre_recid is not None:
                        rec2ctg[pre_recid][1] = new_contig.id
                    rec_id = generate_rec_id(rec)
                    rec2ctg[rec_id] = [None, None]
                    if len(new_seq) > 0:
                        rec2ctg[rec_id][0] = new_contig.id
                    pre_recid = rec_id

                    new_seq = ""    # 更新new_seq
                    connect_rec_ls.append(rec)  # 接头位置
                    
                else:   # 中间的区域    no bridge info
                    ## 
                    if rec.operation.endswith("skip"):# 原数据填充
                        new_seq += ctg_seq[rec.start:rec.end]
                    elif rec.operation.endswith("patch"):# 区域填充
                        new_seq += info_ls[0]   # 接上这段
                        if rec.operation == "asm_patch":
                            gap_ls.append([start, end])
                            patch_ids.add(patch_ls[0])
                    elif rec.operation == "N_fill":
                        '''Cut seq, and add no bridge info'''
                        if len(new_seq) > 0:
                            new_contig = Contig(Contig.generate_id(ctg, new_contig_start, new_contig_end), new_seq)
                            new_contig.add_bound_pos(new_contig_start, new_contig_end)
                            new_contig.add_gap_ls(gap_ls)
                            new_contig_dic[new_contig.id] = new_contig
                            now_ctg_ls.append(new_contig)
                            new_contig_start, new_contig_end = end, end
                            gap_ls = []
                        if pre_recid is not None:
                            rec2ctg[pre_recid][1] = new_contig.id
                        rec_id = generate_rec_id(rec)
                        rec2ctg[rec_id] = [None, None]
                        if len(new_seq) > 0:
                            rec2ctg[rec_id][0] = new_contig.id
                        pre_recid = rec_id

                        new_seq = ""    # 更新new_seq
                        connect_rec_ls.append(rec)  # 接头位置
                    else:
                        raise ValueError
            elif len(op_ls) == 3:   # type: asm/N_fill + N_fill + asm/N_fill
                # 截断区域，两端为可能的接头位置，不需要填充
                '''-------------Cut seq, get new contig-------------'''
                if len(new_seq) > 0:
                    new_contig = Contig(Contig.generate_id(ctg, new_contig_start, new_contig_end), new_seq)
                    new_contig.add_bound_pos(new_contig_start, new_contig_end)
                    new_contig.add_gap_ls(gap_ls)
                    new_contig_dic[new_contig.id] = new_contig
                    now_ctg_ls.append(new_contig)
                    new_contig_start, new_contig_end = end, end
                    gap_ls = []
                    # if new_contig.id == "CM008962.1:50500-350855":
                    #     print(len(new_contig.seq))
                if pre_recid is not None:
                    rec2ctg[pre_recid][1] = new_contig.id
                rec_id = generate_rec_id(rec)
                rec2ctg[rec_id] = [None, None]
                if len(new_seq) > 0:
                    rec2ctg[rec_id][0] = new_contig.id
                pre_recid = rec_id
                new_seq = ""    # 更新new_seq
                connect_rec_ls.append(rec)  # 接头位置
            else:
                print("op value Error")
                raise ValueError
            pre_end = end
        new_seq += ctg_seq[pre_end:ctg_len]  # 加上最后一个rec后面的序列
        if len(new_seq) > 0: 
            new_contig = Contig(Contig.generate_id(ctg, new_contig_start, new_contig_end), new_seq)
            new_contig.add_bound_pos(new_contig_start, new_contig_end)
            new_contig.add_gap_ls(gap_ls)
            new_contig_dic[new_contig.id] = new_contig
            now_ctg_ls.append(new_contig)
            
            if pre_recid is not None:
                rec2ctg[pre_recid][1] = new_contig.id
            rec_id = generate_rec_id(rec)
            rec2ctg[rec_id] = [None, None]
            rec2ctg[rec_id][0] = new_contig.id
            pre_recid = rec_id
        # exit(0)
        ## 第二轮遍历获取全部的桥接信息     
        '''
        桥接位置
        1、左右边界。len(op_ls) =0 and op == asm_patch
        2、中间存在断点的区域 len(op) = 3
        '''
        # print("Size:, ctg_num:{}, connect_num:{}".format(len(now_ctg_ls), len(connect_rec_ls)))
        # print(rec2ctg)
        for contig in now_ctg_ls:
            print(contig.id)
        for idx, rec in enumerate(connect_rec_ls):
            start = rec.start
            end = rec.end
            op_ls = rec.operation.split(",")
            info_ls = rec.info.split(",")
            # 
            left_candidate_ls = []
            right_candidate_ls = []
            left_bound = max(start - min_dis, 1)
            right_bound = min(end + min_dis, ctg_len - 1)
            # 
            left_gap_ls, right_gap_ls = [], []
            rec_id = generate_rec_id(rec)
            if rec2ctg[rec_id][0] is not None:
                left_gap_ls = new_contig_dic[rec2ctg[rec_id][0]].gap_ls
            if rec2ctg[rec_id][1] is not None:
                right_gap_ls = new_contig_dic[rec2ctg[rec_id][1]].gap_ls
            l_contig_id = rec2ctg[rec_id][0]
            r_contig_id = rec2ctg[rec_id][1]
            # if idx > 0 and idx < len(connect_rec_ls) - 1:
            #     left_gap_ls = now_ctg_ls[idx - 1].gap_ls
            #     right_gap_ls = now_ctg_ls[idx].gap_ls
            # elif idx == 0:
            #     right_gap_ls = now_ctg_ls[idx].gap_ls
            # elif idx == len(connect_rec_ls) - 1:
            #     left_gap_ls = now_ctg_ls[idx - 1].gap_ls
            # print("{}:{}-{}, l-Gap_ls:{}, r-Gap_ls:{}".format(rec.chr_id, rec.start, rec.end, left_gap_ls, right_gap_ls))
            for read in bam_reader.fetch(rec.chr_id, rec.start, rec.end):
                if read.is_secondary or read.mapping_quality < min_MQ or read.query_alignment_length < min_align_length:
                    continue
                cigar = read.cigartuples
                left_clip, right_clip = 0, 0
                if cigar[0][1] == 5 or cigar[0][1] == 4:
                    left_clip = cigar[0][1]
                if cigar[-1][1] == 5 or cigar[-1][1] == 4:
                    right_clip = cigar[-1][1]
                if len(left_gap_ls) != 0:
                    left_gap_pos = max(1, left_gap_ls[-1][1])   # 取左边块最后一个asm reg的右端点
                else:   # 左边块的左bound
                    # left_gap_pos = 1 if idx == 0 else new_contig_dic[l_contig_id].left_bound_pos
                    left_gap_pos = 1 if l_contig_id is None else new_contig_dic[l_contig_id].left_bound_pos
                if len(right_gap_ls) != 0:
                    right_gap_pos = min(ctg_len - 1, right_gap_ls[0][0])   # 取右边块第一个asm reg的左端点
                else:   # 右边块的右bound
                    # right_gap_pos = ctg_len - 1 if idx == len(connect_rec_ls) - 1 else new_contig_dic[l_contig_id].right_bound_pos
                    right_gap_pos = ctg_len - 1 if r_contig_id is None else new_contig_dic[r_contig_id].right_bound_pos
                ## 
                if read.reference_start <= left_bound:    # patch左边界的点，如果比对起始在左边块内部，则需要没有大型left_clip
                    if read.reference_start > left_gap_pos and left_clip > min_clip: continue
                    left_candidate_ls.append(read)
                if read.reference_end >= right_bound:        # patch右边界的点，如果比对结束在右边块内部，则需要没有大型right_clip
                    if read.reference_end < right_gap_pos and right_clip > min_clip: continue
                    right_candidate_ls.append(read)
            if len(op_ls) == 1:
                if op_ls[0] == "asm_patch":
                    if start == 0:  # contig start, left bound
                        best_bridge = get_best_bridge(right_candidate_ls)
                        if best_bridge is not None:
                            # best_bridge = pysam.AlignedSegment()
                            bridge_id = best_bridge.query_name
                            # contig_id = now_ctg_ls[0].id
                            contig_id = r_contig_id
                            bridge = bridge_dic[bridge_id]
                            strand = "-" if best_bridge.is_reverse else "+"
                            bridge_start, bridge_end = get_raw_pos(best_bridge)
                            # bridge = Bridge()
                            # bridge.add_connect(contig_id, "l", now_ctg_ls[0].left_bound_pos, bridge_start, bridge_end, strand, info_ls[0], best_bridge)
                            bridge.add_connect(contig_id, "l", new_contig_dic[contig_id].left_bound_pos, bridge_start, bridge_end, strand, best_bridge)
                            new_contig.add_l_connect(bridge_id, strand, best_bridge.query_alignment_start, best_bridge.query_alignment_length)
                    elif end >= ctg_len - 1:    # contig end, right bound
                        best_bridge = get_best_bridge(left_candidate_ls)
                        if best_bridge is not None:
                            bridge_id = best_bridge.query_name
                            contig_id = l_contig_id
                            bridge = bridge_dic[bridge_id]
                            strand = "-" if best_bridge.is_reverse else "+"
                            bridge_start, bridge_end = get_raw_pos(best_bridge)
                            bridge.add_connect(contig_id, "r", new_contig_dic[contig_id].right_bound_pos, bridge_start, bridge_end, strand, best_bridge)
                            new_contig.add_r_connect(bridge_id, strand, best_bridge.query_alignment_start, best_bridge.query_alignment_length)
            elif len(op_ls) == 3:
                best_bridge = get_best_bridge(right_candidate_ls)
                # l_contig_id = now_ctg_ls[idx - 1].id
                # r_contig_id = now_ctg_ls[idx].id
                l_best_bridge = get_best_bridge(left_candidate_ls)
                r_best_bridge = get_best_bridge(right_candidate_ls)
                if l_best_bridge is not None:   # connect/bridge to l_contig right bound    "r"
                    best_bridge = l_best_bridge
                    bridge_id = best_bridge.query_name
                    bridge = bridge_dic[bridge_id]
                    strand = "-" if best_bridge.is_reverse else "+"
                    bridge_start, bridge_end = get_raw_pos(best_bridge)
                    # bridge.add_connect(l_contig_id, "r", now_ctg_ls[idx - 1].right_bound_pos, bridge_start, bridge_end, strand, info_ls[0], best_bridge)
                    bridge.add_connect(l_contig_id, "r", new_contig_dic[l_contig_id].right_bound_pos, bridge_start, bridge_end, strand, best_bridge)
                    new_contig.add_r_connect(bridge_id, strand, best_bridge.query_alignment_start, best_bridge.query_alignment_length)
                if r_best_bridge is not None:   # connect/bridge to r_contig left bound     "l"
                    best_bridge = r_best_bridge
                    bridge_id = best_bridge.query_name
                    bridge = bridge_dic[bridge_id]
                    strand = "-" if best_bridge.is_reverse else "+"
                    bridge_start, bridge_end = get_raw_pos(best_bridge)
                    # bridge.add_connect(r_contig_id, "l", now_ctg_ls[idx].left_bound_pos, bridge_start, bridge_end, strand, info_ls[2], best_bridge)
                    bridge.add_connect(r_contig_id, "l", new_contig_dic[r_contig_id].left_bound_pos, bridge_start, bridge_end, strand, best_bridge)
                    new_contig.add_r_connect(bridge_id, strand, best_bridge.query_alignment_start, best_bridge.query_alignment_length)
            else:
                # print(rec)
                # print(op_ls)
                raise ValueError

    candidate_scaff_dic = {}
    cnt = 0
    candidate_scaff_id_set = set()
    print("---------------------------Start collect bridge--------------------------------")
    
    for bridge_id, bridge in bridge_dic.items():
        # bridge = Bridge()
        if len(bridge.connect_info_ls) == 0:
            continue
        elif len(bridge.connect_info_ls) == 1:
            contig_id = bridge.connect_info_ls[0]["contig_id"]
            bound = bridge.connect_info_ls[0]["bound"]
            align_info = bridge.connect_info_ls[0]["align_info"]
            contig = new_contig_dic[contig_id]
            contig.solve_bound(bound, bridge_id, get_patch_seq(align_info, bridge, bound, contig.get_bound_pos(bound)))
            patch_ids.add(bridge.bridge_id)
        else:
            bridge.connect_info_ls.sort(key=lambda connect_info: connect_info["raw_bridge_start"])
            candidate_scaff_dic[bridge_id] = bridge
            cnt += 1
            print("\n", cnt, bridge_id, bridge)
            # print("Bridge: ", ",".join(id_ls), bridge_id)
            for connect_info in bridge.connect_info_ls:
                candidate_scaff_id_set.add(connect_info['contig_id'])
            patch_ids.add(bridge.bridge_id)
    with open(work_dir + "/candidate_scaff.bed", "w") as f:
        for id in candidate_scaff_id_set:
            ctg, reg = id.split(":")
            start, end = reg.split("-")
            f.write("{}_{}_{}\n".format(ctg, start, end))
    if not exe_scaffold:
        final_dic = Contig.write_contig(fa_out, new_contig_dic)
        return [], patch_ids, final_dic
    # -------------------------------Perform scaffolding in several rounds-------------------------------
    # 多轮处理
    '''
    0、2个的类型
    1、源自同染色体的片段桥接
    2、源自不同chr
    '''
    print("---------------------------------------Perform scaffolding in several rounds---------------------------------------")
    print("Candidate_scaff:{}".format(candidate_scaff_dic.keys()))
    connect_info_ls = scaffold(candidate_scaff_dic, new_contig_dic)
    final_dic = Contig.write_contig(fa_out, new_contig_dic)
    # 
    # Connect_info.write_connect_info(connect_info_ls, file = "/public/home/hpc214712170/Test/tests/bench/reinhardtii_hifi/my_pipe2/step3_SV_consensus/connect.txt")
    return connect_info_ls, patch_ids, final_dic


if __name__ == "__main__":
    # 
    import yaml
    # import fasta_parser
    config_f = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/Configs/Config.yaml"
    # work_dir = "/public/home/hpc214712170/Test/tests/bench/reinhardtii_ont/my_pipe2"
    # asm_fa = work_dir + "/step3_SV_consensus/merge_denovo/shasta/Assembly.fasta"
    work_dir = "/public/home/hpc214712170/Test/tests/bench/reinhardtii_hifi/my_pipe2"
    asm_fa = work_dir + "/step3_SV_consensus/merge_denovo/hifiasm/out.p_ctg.fasta"
    bam = work_dir + "/step3_SV_consensus/merge_denovo/denovoasm_to_ref.sorted.bam"
    ref_fa = work_dir + "/step3_SV_consensus/semi_consensus.fa"
    rec_file = work_dir + "/step3_SV_consensus/semi_consensus.bed"
    # fa_out = "work_dir + /step3_SV_consensus/semi.fasta"
    fa_out = work_dir + "/step3_SV_consensus/scaff.fasta"
    ref_dic = fasta_parser.read_sequence_dict(ref_fa)
    asm_fa_dic = fasta_parser.read_sequence_dict(asm_fa)
    rec_ls = Record.read_record(rec_file)
    candidate_rec_dic = defaultdict(list)
    ctg_ls = list(ref_dic.keys())
    with open(config_f, "r") as f:
        config = yaml.safe_load(f.read())     # 获取部分参数
    for rec in rec_ls:
        candidate_rec_dic[rec.chr_id].append(rec)
    rec_dic = defaultdict(list)
    for ctg in ctg_ls:
        ctg, ctg_rec_ls = get_rec_by_asm(ctg, bam, candidate_rec_dic[ctg], asm_fa_dic, config["step3"]["solve_by_denovo"])
        rec_dic[ctg] = ctg_rec_ls
    consensus(bam, ref_dic, rec_dic, asm_fa_dic, fa_out, work_dir, config, True)  # bam, ref_dic, rec_dic, asm_fa_dic, fa_out, work_dir, config, exe_scaffold
    pass

# contig_ls = [["chr1", 10, 20, "a"*10], ["chr1", 30, 40, "b"*10], ["chr2", 10, 20, "c"*10], ["chr2", 30, 40, "d"*10]]
# dic = {"seq1":"sssssssssssrrrrrrrrrr", "seq2":"iiiiiiiiiiiiiiiiiiinnnnnnnnnnnnnnnnn", "seq3":"mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm"}
# bridge_seqs = Bridge_seqs(dic)
# # bridge_seqs.add_bridge("seq1", "chr1:10-20", "left")

# # contig = Contig("chr1", 10, 20, "a"*10)
# # print(contig.seq)
# print(list(dic.keys()))

# contig1 = Contig("chr1", 10, 20, "a"*10)
# contig2 = Contig("chr1", 30, 40, "b"*10)
# contig3 = Contig("chr2", 10, 20, "c"*10)
# contig4 = Contig("chr2", 30, 40, "d"*10)
# # print(contig3.seq)

# contig1.add_l_connect("seq3", "-", 10, 5)
# bridge_seqs.add_bridge("seq3", contig1.id, "left")
# contig1.add_r_connect("seq1", "+", 20, 5)
# bridge_seqs.add_bridge("seq1", contig1.id, "right")

# contig2.add_l_connect("seq2", "+", 10, 5)
# bridge_seqs.add_bridge("seq2", contig2.id, "left")

# contig3.add_r_connect("seq3", "-", 20, 5)
