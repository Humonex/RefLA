
# from create_consensus_by_bed import fasta_parser
# from create_consensus_by_bed import Utils
import fasta_parser, Utils

def get_semi_rec(candidate_op_ls):
    '''
    1、process one ctg
    2、仅处理read patch区域
    '''
    final_op_ls = []
    for rec in candidate_op_ls:
        if rec.operation.endswith("reads"):
            rec.add_operation("reads_patch")
            final_op_ls.append(rec)
        elif rec.operation.endswith("asm"):
            final_op_ls.append(rec)
        else: raise ValueError

    return final_op_ls

def apply_read_patch(chr_id, rec_ls, ref_dic:dict):
    ''' 提供参考序列一条contig，及对其的所有操作集，将处理后的contig写入consensus_fa_dic中 '''
    consensus_fa_dic = {}
    ref_seq_in = ref_dic[chr_id]

    new_rec_ls = []
    # 
    rec_ls.sort(key=lambda rec:int(rec.start))  # 防止有乱序的
    new_seq = ""
    pre_end = 0    # pre rec end
    for rec in rec_ls:
        new_rec = rec
        # 
        op_ls = rec.operation.split(",")
        info_ls = rec.info.split(",")
        new_seq += ref_seq_in[pre_end:rec.start]    # 加上上一段rec的末尾到现在rec前面的中间这段序列，相当于是保留的，bias不变
        new_rec.start = len(new_seq)
        ## process
        for i, op in enumerate(op_ls):
            if op == "reads_patch":
                new_seq += info_ls[i]   # 接上这段
                # new_rec.skip = True
                new_rec.operation = "skip"
            else:   # 保留，并更改坐标
                pass
        new_rec.end = len(new_seq)
        # 
        pre_end = rec.end
        new_rec_ls.append(new_rec)
    consensus_fa_dic[chr_id] = new_seq
    return new_rec_ls, consensus_fa_dic

'''file = "/public/home/hpc214712170/Test/tests/yeast_ont/test/my_pipe/step3_SV_consensus/semi_consensus.bed"
fa_in = "/public/home/hpc214712170/Test/tests/yeast_ont/test/my_pipe/step3_SV_consensus/semi_consensus.fa"
fa_dic = fasta_parser.read_sequence_dict(fa_in)

rec_ls = Utils.Record.read_record(file)
for rec in rec_ls:
    if rec.operation == "reads_patch":
        info = rec.info
        if info == fa_dic[rec.chr_id][rec.start:rec.end]:
            print("Check right")
        else:
            print(len(info), rec.end - rec.start)
            # print("Error")

fa_in2 = "/public/home/hpc214712170/Test/tests/yeast_ont/test/my_pipe/corrected_ref/reference.fasta"
file2 = "/public/home/hpc214712170/Test/tests/yeast_ont/test/my_pipe/step3_SV_consensus/candidate_op/candidate_op.bed"
rec_ls2 = Utils.Record.read_record(file2)
fa_dic2 = fasta_parser.read_sequence_dict(fa_in2)
pre_end1 = 0
pre_end2 = 0
for i in range(len(rec_ls)):
    rec1 =rec_ls[i]
    rec2 = rec_ls2[i]
    # rec2.print_info()
    # if rec1.operation == "reads_patch":
    #     print(rec1.info ==  rec2.info)
    print(fa_dic[rec1.chr_id][pre_end1:rec1.start] == fa_dic2[rec2.chr_id][pre_end2:rec2.start])
    pre_end1 = rec1.end
    pre_end2 = rec2.end'''

'''from collections import defaultdict
asm_N_fill_file = "/public/home/hpc214712170/Test/tests/bench/reinhardtii_hifi/can.bed"
with open(asm_N_fill_file, "r") as f:
    dic = defaultdict(list)
    dic2 = defaultdict(list)
    for line in f:
        fields = line.strip().split("\t")
        chr, start, end, op, patch = fields
        start, end = int(start), int(end)
        patch_ls = patch.split(",")
        dic[chr].append([chr, start, end, op, patch])
        idx = 1
        for patch_info in patch_ls:
            if patch_info != ".":
                reg_id = chr + ":" + str(start) + "-" + str(end) + "-" + str(idx)
                dic2[patch_info].append(reg_id)
            idx += 1
    for key in dic2.keys():
        # print(key)
        if len(dic2[key]) > 1:
            print(key, dic2[key])'''
'''


ptg000058l {'CM008963.1:3622156-3766256', 'CM008963.1:2863659-2991359'}
ptg000253l {'CM008972.1:2180442-2496942', 'CM008972.1:2563620-2710020', 'CM008963.1:2863659-2991359'}
ptg000251l {'CM008963.1:3622156-3766256', 'CM008972.1:1304912-1480912'}
ptg000072l {'CM008963.1:6956120-7037620', 'CM008970.1:3099466-3209966'}
ptg000077l {'CM008963.1:6956120-7037620', 'CM008970.1:4006129-4113529'}
ptg000009l {'CM008978.1:4958828-5162728', 'CM008963.1:8897087-9149587'}
ptg000010l {'CM008978.1:4697620-4836420', 'CM008969.1:183513-342013'}
ptg000221l {'CM008970.1:4006129-4113529', 'CM008970.1:2651870-2694370'}
ptg000074l {'CM008970.1:3099466-3209966', 'CM008970.1:2651870-2694370'}
ptg000248l {'CM008972.1:2180442-2496942', 'CM008972.1:1203684-1264984', 'CM008972.1:1304912-1480912'}
ptg000263l {'CM008973.1:277675-401675', 'CM008973.1:793133-1034733'}
ptg000264l {'CM008973.1:277675-401675', 'CM008973.1:643972-728272'}
ptg000375l {'CM008978.1:4958828-5162728', 'CM008978.1:4697620-4836420'}，ptg000058l比对到了后面的两个位置，说明这两个位置有关联，请写一个基因组Scaffolding算法 
'''
'''
我是做生物基因组组装的，我现在有组装结果，和一堆桥接序列，我将桥接序列比对到了组装结果中。
组装结果中存在一些我标记出来的区域，这些区域我理解为是缺口，我希望通过桥接序列来修复这些缺口。
1、目前我实现了最简单的情况，也就是区域前后的序列顺序一致，因此中间的桥接序列能够完全跨越缺口。
2、对于区域前后顺序不一致的情况，从比对上来看是缺口的左右两侧比对上了不一样的桥接序列。因此需要对这些序列进行重新排序链接等操作。例如一个例子是gap1左右比对上了seq1和seq2，gap2左右比对上了seq1和seq3，因此需要将gap1的左端连上seq1再脸上gap2的左端。
你给出一个大致的代码

'''
'''
我不觉得这样合理，我觉得不应该按gap处理。
因为我已经处理了简单的能用一条序列横跨一个gap的情况。
现在处理的是两条序列比对到一个gap两侧的情况，所以这时候gap两侧的序列不应该衔接在一起。
我的目的是通过遍历所有桥接序列，如果他们连接上了不同的gap侧，则将这gap两侧的序列通过桥接序列连接上。
事实上，应该一开始根据未解决gap的位置，对基因组序列进行切割，将每个块分配一个id，
同时记录下每个块左右两侧连接的桥接序列id，然后进行后续处理
'''
class Mis_rec:
    def __init__(self, ctg, start, end, seg_id_ls, unalign_len, mis_info_ls) -> None:
        self.ctg = ctg
        self.start = start
        self.end = end
        self.seg_id_ls = seg_id_ls
        self.unalign_len = unalign_len
        self.mis_info_ls = mis_info_ls
        ## 
        self.ref = None
        self.ref_start = None
        self.ref_end = None
    
    def Add_ref(self, ref1, start1, end1, ref2, start2, end2, is_bound):
        if is_bound:
            self.ref = ref1
            self.ref_start = min(start1, end1)
            self.ref_end = max(start1, end1)
        elif ref1 == ref2:
            self.ref = ref1
            start1, end1 = min(start1, end1), max(start1, end1)
            start2, end2 = min(start2, end2), max(start2, end2)
            if start1 > start2:
                start1, end1, start2, end2 = start2, end2, start1, end1
            self.ref_start = min(end1, start2)
            self.ref_end = max(end1, start2)
        else:
            self.ref = []
    @staticmethod
    def write_rec(mis_rec_ls, bed):
        with open(bed, "w") as fo:
            fo.write("#ctg\tstart\tend\tseg_ids\tunalign_len\tmis_info_ls\n")
            for mis_rec in mis_rec_ls:
                # mis_rec = Mis_rec()
                fo.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(mis_rec.ctg, str(mis_rec.start), str(mis_rec.end), ",".join('%s' %id for id in mis_rec.seg_id_ls), str(mis_rec.unalign_len), ",".join(mis_rec.mis_info_ls)))

class Align_Info:
    def __init__(self, contig, start, end, ref, ref_start, ref_end, IDY, seg_id) -> None:
        self.contig = contig
        self.start = start
        self.end = end
        ##
        self.ref = ref
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.IDY = IDY
        ## 
        self.seg_id = seg_id
        self.is_left = False
        self.is_right = False
        # 
        self.align_orient = "+"
        if self.start > self.end:
            self.align_orient = "-"
            self.start, self.end = self.end, self.start

tsv = "/public/home/hpc214712170/Test/tests/bench/reinhardtii_hifi/my_pipe2/step3_SV_consensus/quast_eval_large/contigs_reports/all_alignments_semi.tsv"
ctg_bed = "/public/home/hpc214712170/Test/tests/bench/reinhardtii_hifi/my_pipe2/step3_SV_consensus/candidate_ctg.bed"
fou = "/public/home/hpc214712170/Test/tests/bench/reinhardtii_hifi/my_pipe2/step3_SV_consensus/candidate_ctg_align.bed"
ctg_ls = []
print("Start")
with open(ctg_bed, "r") as f:
    for line in ctg_bed:
        ctg_ls.append(line.strip())


with open(tsv, "r") as fin:
    print("Start")
    line = fin.readline()
    # print(line)
    line = fin.readline()
    ctg_align_info_ls = []
    ctg_mis_info_ls = []
    mis_rec_ls = []
    seg_id = 0
    bias = 10
    unalign_num = 0
    correct_num = 0
    while line:
        # print(line)
        if line.startswith("CONTIG"):
            ## 处理last contig info
            fields = line.strip().split("\t")
            contig_type = fields[3]
            ctg, ctg_len = fields[1:3]
            ctg_len = int(ctg_len)
            if len(ctg_mis_info_ls) == 0 and (contig_type == "mis_unaligned" or contig_type == "unaligned" or contig_type == "correct_unaligned"):
                print("---------{} is {}, no mis ---------".format(ctg, contig_type))
                unalign_num += 1
            elif len(ctg_mis_info_ls) == 0 and contig_type == "correct":
                print("---------{} is {}---------".format(ctg, contig_type))
                correct_num += 1
            else:   # misassembled contig、correct contig(has local mis)
                if contig_type == "mis_unaligned" or contig_type == "unaligned" or contig_type == "correct_unaligned":
                    print("{}, {}:{}".format(ctg, contig_type, ctg_mis_info_ls))
                print("{} start, type:{}".format(ctg, contig_type))
                ctg_align_info_ls[0].is_left = True
                ctg_align_info_ls[-1].is_right = True
                # 
                pre_align_info = None
                mis_idx = 0
                for align_info in ctg_align_info_ls:
                    # align_info = Align_Info()
                    if align_info.is_left:  # Type1: left
                        mis_rec = Mis_rec(ctg, 1, align_info.start, [align_info.seg_id], align_info.start - 1, ['left_unalign'])
                        # mis_rec.
                        # print("left unalign:1-{}, length:{}".format(str(align_info.start), str(align_info.start - 1)))
                        # print(mis_rec.__dict__)
                        mis_rec_ls.append(mis_rec)
                    elif ctg_mis_info_ls[mis_idx][0] == "unknown":  # unknown type
                        mis_idx += 1
                        continue
                    else:
                        # Type2
                        mis_info = ctg_mis_info_ls[mis_idx]
                        # pre_align_info = Align_Info()
                        # align_info = Align_Info()
                        mis_start = max(pre_align_info.end - bias, pre_align_info.start)
                        mis_end = min(align_info.start + bias, align_info.end)
                        seg_id_ls = [pre_align_info.seg_id, align_info.seg_id]
                        unalign_len = align_info.start - pre_align_info.end
                        mis_rec = Mis_rec(ctg, mis_start, mis_end, seg_id_ls, unalign_len, mis_info)
                        if mis_rec.start > mis_rec.end:
                            print("pre:{}-{}, now:{}-{}".format(pre_align_info.start, pre_align_info.end, align_info.start, align_info.end))
                            print("Error")
                            print(pre_align_info.__dict__, align_info.__dict__)
                            raise ValueError
                        mis_idx += 1
                        mis_rec_ls.append(mis_rec)
                        # Type3
                        if align_info.is_right: # 右边界的末端-end
                            mis_rec = Mis_rec(ctg, align_info.end, ctg_len, [align_info.seg_id], ctg_len - align_info.end, ['right_unalign'])
                            # print("right unalign:{}-{}, length:{}".format(str(align_info.end), ctg_len, str(ctg_len - align_info.end)))
                            mis_rec_ls.append(mis_rec)
                    pre_align_info = align_info
            # 
            print(ctg + ":", ctg_len, "finished\n")
            # break
            seg_id = 0
            ctg_align_info_ls = []
            ctg_mis_info_ls = []
            line = fin.readline()
        else:
            if line.strip() == "unknown": 
                line = fin.readline()
                continue
            ## parse alignsegment info
            align_seg_fields = line.strip().split("\t") 
            # print(align_seg_fields)
            # 181', '12634', '140', '12550', 'chrI', 'NC_001133.9', '96.82', '', 'True
            # S1      E1      S2      E2      Reference       Contig  IDY     Ambiguous       Best_group
            # print(align_seg_fields)
            ref_start, ref_end, ctg_start, ctg_end, ref, contig, IDY = align_seg_fields[:7]
            # if int(ctg_start) > int(ctg_end):
            #     print(ref_start, ref_end, ctg_start, ctg_end, ref, contig, "is reverse")
            # print(ref_start, ref_end, ctg_start, ctg_end, ref, contig, IDY, seg_id)
            # contig, start, end, ref, ref_start, ref_end, IDY, seg_id
            align_info = Align_Info(contig, int(ctg_start), int(ctg_end), ref, int(ref_start), int(ref_end), float(IDY), seg_id) 
            # print(align_info.__dict__)
            ctg_align_info_ls.append(align_info)
            line = fin.readline()
            if line.startswith("CONTIG"):   # End
                continue
            else:   # parse mis info
                mis_info = line.strip().split(",")
                # mis_info = line.strip()
                # print(mis_info)
                ctg_mis_info_ls.append(mis_info)
                line = fin.readline()
            seg_id += 1
    new_mis_rec_ls = []
    for mis_rec in mis_rec_ls:
        if mis_rec.ctg not in ctg_ls: continue
        if mis_rec.mis_info_ls[0] == 'left_unalign' or mis_rec.mis_info_ls[0] == 'right_unalign':
            # if mis_rec.unalign_len > 10000:
            #     new_mis_rec_ls.append(mis_rec)
            continue
        elif mis_rec.mis_info_ls[0].startswith("fake") or mis_rec.mis_info_ls[0].startswith("indel"):
            continue
        else:
            new_mis_rec_ls.append(mis_rec)
    print("-----------Sum:---------\nmis_num:", len(new_mis_rec_ls))
    print("unalign contig:{}".format(unalign_num))
    print("correct contig:{}".format(correct_num))
    Mis_rec.write_rec(new_mis_rec_ls, fou)