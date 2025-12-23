from Bio import Seq
import pysam
from create_consensus_by_bed import fasta_parser

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
        raw_存了原始位置
        bridge_start: bridge 开始比对的位置
        bridge_end: 
        '''
        bridge_bound_pos = convert_reference_pos_to_raw_pos(align_info, bound_pos, True)    # 反向的比对则是反向互补链中的位置
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
        # 
        self.is_scaffold = False
        # 用于存储桥接的contig_id信息，在对scaffolder进行操作的时候，方便记录原始片段的相关变化
        self.contig_id_ls = [contig_id]
        self.bridge_id_ls = []
    def get_reverse_seq(self):
        return Seq.reverse_complement(self.seq)
    def add_bound_pos(self, start, end):
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
    def solve_bound(self, bound, seq):
        if bound == "l" and self.l_solved:
            self.seq = seq + self.seq
            self.l_solved = True
        elif bound == "r" and self.r_solved:
            self.seq = self.seq + seq
            self.r_solved = True
        else:
            raise ValueError
    def kill(self):
        self.died = True
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
    @staticmethod
    def scaffolding(contig1, contig2, bridge_id, new_seq):
        contig1.kill()
        contig2.kill()
        contig_id = contig1.id + "_" + contig2.id
        contig = Contig(contig_id, new_seq)
        contig.contig_id_ls = contig1.contig_id_ls + contig2.contig_id_ls
        contig.bridge_id_ls = contig1.bridge_id_ls + [bridge_id] +  contig2.bridge_id_ls
        contig.is_scaffold = True
        contig.l_solved = contig1.l_solved
        contig.r_solved = contig2.r_solved
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

class UnionFind:
    def __init__(self):
        self.root = {}
        # self.next_id = 0

    def find(self, x):
        # 如果 x 不在映射中，那么它是最新的ID，直接返回 x
        if x not in self.root:
            return x
        # 路径压缩：更新沿路径上的每个节点，使其指向最终的根节点
        if self.root[x] != x:
            self.root[x] = self.find(self.root[x])
        return self.root[x]

    def union(self, x, y, new_id):
        # new_id = 'id' + str(self.next_id)  # 生成新的ID
        # self.next_id += 1
        rootX = self.find(x)
        rootY = self.find(y)
        # 将 rootY 和 rootX 的合并后的新ID存储起来
        self.root[rootX] = self.root[rootY] = new_id
        # 返回新生成的ID

    def connected(self, x, y):
        return self.find(x) == self.find(y)

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
# # 示例使用
# uf = UnionFind()
# new_id1 = uf.union('id1', 'id2') # 假设生成 'id0'
# new_id2 = uf.union('id0', 'id3') # 假设生成 'id1'
# # 使用旧的id去检索，应该都能得到最新的id
# print(uf.find('id1'))  # 输出 'id1', 即最新生成的id
# print(uf.find('id2'))  # 输出 'id1'
# print(uf.find('id3'))  # 输出 'id1'
# print(uf.find('id0'))  # 输出 'id1'