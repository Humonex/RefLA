
from collections  import defaultdict
bed_in = "/public/home/hpc214712170/Test/tests/bench/reinhardtii_hifi/my_pipe2/step3_SV_consensus/candidate_ctg_align.bed"
bed_out = "/public/home/hpc214712170/Test/tests/bench/reinhardtii_hifi/my_pipe2/step3_SV_consensus/candidate_ctg_align_redunt.bed"
pair_bed = "/public/home/hpc214712170/Test/tests/bench/reinhardtii_hifi/my_pipe2/step3_SV_consensus/candidate_ctg_pair.bed"

pair_ls = []
# with open(pair_bed, "r") as f:
#     for line in f:
#         pair = line.strip().split("\t")
#         print(pair)
#         pair_ls.append(pair)

# with open(bed_in, "r") as f, open(bed_out, "w") as fo:
#     ctg_ls = []
#     dic = defaultdict(list)
#     dic2 = defaultdict(list)
#     for line in f:
#         # S1 E1 S2 E2 Reference Contig IDY
#         S1, E1, S2, E2, Reference, Contig = line.strip().split("\t")[:6]
#         S1, E1, S2, E2 = int(S1), int(E1), int(S2), int(E2)
#         S1, E1 = min(S1, E1), max(S1, E1)
#         dic[Contig].append([Reference, S1, E1])
#     for ctg, align_ls in dic.items():
#         # print(ctg, align_ls)
#         print("\n----------------" + ctg + "--------------------")
#         fo.write("----------{}----------\n".format(ctg))
#         align_dic = defaultdict(list)
#         for align_reg in align_ls:
#             align_dic[align_reg[0]].append(align_reg)
#         for ref, align_ls in align_dic.items():
#             # print(ctg, ref, align_ls)
#             align_ls.sort()
#             merge_start, merge_end = align_ls[0][1], align_ls[-1][2]
#             # align_dic[ref] = [ref, merge_start, merge_end]
#             print(ref, merge_start, merge_end)
#             fo.write("{}\t{}\t{}\t{}\n".format(ref, merge_start, merge_end, ctg))
#         dic2[ctg].append([ref, merge_start, merge_end])
    
#     print("Pair align info")
#     for pair in pair_ls:
#         print("\n")
#         print(pair)
#         for ctg in pair:
#             ctg = ctg[:-2]
#             print(ctg, dic2.get(ctg))
# CP097821.1:1376714-2075573


# from collections import Counter
# # 假设有一个列表
# lst = [1, 2, 2, 3, 3, 3, 4, 5, 5, 5, 5]

# # 创建Counter对象
# counter = Counter(lst)
# print(counter)
# # 获取按计数降序排序的列表
# sorted_elements_by_count = counter.most_common()

# for element, count in sorted_elements_by_count:
#     print(f'Element: {element}, Count: {count}')
# print(sorted_elements_by_count[0][1])
# for idx, rec in enumerate(lst):
#     print(idx, rec)

# import pysam
# bam = "/public/home/hpc214712170/Test/tests/bench/reinhardtii_hifi/my_pipe2/step3_SV_consensus/merge_denovo/denovoasm_to_ref.sorted.bam"
# reg = "CM008973.1:9529991-9623827"
# # reg = "CM008973.1:9474359-9523378"
# # reg = "CM008973.1:2206570-2206575"
# reg = "CM008973.1:4685631-4685635"
# # target = "ptg000286l"
# # target = "ptg000287l"
# # target = "ptg000270l"
# target  = "ptg000179l"
# bam_reader = pysam.AlignmentFile(bam, "rb", index_filename=bam+".bai")
# for read in bam_reader.fetch(region=reg):
#     if read.query_name == target:
#         if read.is_reverse:
#             print(read.query_name, read.reference_start, read.reference_end)
#             print(read.query_alignment_start, read.query_alignment_end) # 按参考的方向，读数进行反向互补，在反向互补链上的位置
#             print(read.cigarstring) # cigar按参考的方向给出，若是反向互补，也是反向互补之后的
        
'''
- primary
15422 69024
15422S7303M3I1703M1D19805M30I3171M1I680M1D133M1D3526M3I10634M2D1297M1D68M2D692M1D1854M1D1877M1I293M1I174M1D353M

- supp
0 105874
1185M1I50M.....I112M883653H

- supp
0 1636
108898H1636M135801H

19349 27571
19349S188M1D166M1D100M1I484M1I21M92D14M1I611M1D936M1I1213M744I1125M745I542M1I55M200I3M2I1068M23430S
0 8228
39942H961M91D15M1I1199M1I56M1I114M1I241M1D44M1I230M1I876M744I36M1I113M1I165M1I811M745I501M1D95M200I3M2I1068M2831H

'''
# ls1 = [1, 1, 1]
# ls2 = [2, 2, 2]
# ls3 = [ls1, ls2]
# ls1, ls2 = ls2, ls1
# print(ls1, ls2, ls3)
# print(type(ls3[0]))

s = "{}:{}-{}".format("ch", 1, 3)
print(s)