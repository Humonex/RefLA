from collections import defaultdict

species = "thaliana", "chm13"
# fai = "/public/home/hpc214712170/Test/mis_detect/asm/thaliana/data/flye_ont.fa.fai"
# file = "/public/home/hpc214712170/Test/mis_detect/asm/thaliana/quast_large/contigs_reports/mis.bed"
fai = "/public/home/hpc214712170/Test/mis_detect/asm/thaliana_hifi/data/flye_hifi.fa.fai"
file = "/public/home/hpc214712170/Test/mis_detect/asm/thaliana_hifi/my_pipe/step2_candidate_regions/filtered/merge.bed"
# file = "/public/home/hpc214712170/Test/mis_detect/asm/thaliana_hifi/my_pipe/step2_candidate_regions/filtered/false.bed"
ctg_cnt = defaultdict(int)
dic = defaultdict(int)
mis_len = defaultdict(int)
with open(fai, "r") as f:
    for line in f:
        fields = line.strip().split("\t")
        ctg, ctg_len = fields[:2]
        ctg_len = int(ctg_len)
        dic[ctg] = ctg_len
with open(file, "r") as f:
    for line in f:
        if line.startswith("#"):continue
        fields = line.strip().split("\t")
        ctg, start, end = fields[:3]
        start, end = int(start), int(end)
        ctg_cnt[ctg] += 1
        mis_len[ctg] += (end - start)
for ctg in ctg_cnt.keys():
    print("{}:{}, cnt:{}, mis_len:{}".format(ctg, dic[ctg], ctg_cnt[ctg], mis_len[ctg]))
ls = []
portion_ls = []
for ctg in ctg_cnt.keys():
    # ls.append([dic[ctg], ctg_cnt[ctg]])
    # ls.append(dic[ctg])
    ls.append([dic[ctg], ctg, ctg_cnt[ctg], mis_len[ctg] / dic[ctg]])
    portion_ls.append(mis_len[ctg] / dic[ctg])
ls.sort()
print(ls)
# portion_ls.sort()
# print(portion_ls)
# import pysam
# ref="/public/home/hpc214712170/Test/mis_detect/asm/chm13_ont/flye/data/flye_ont.fa"

# pysam.faidx(ref)

a = input()
print(a)
