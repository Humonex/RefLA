import os
import time
from collections import defaultdict
# file_in = "/public/home/hpc214712170/Test/tests/sativa_hifi/my_pipe2/step1_mapping/reads.pos"
# file_out = "/public/home/hpc214712170/Test/tests/sativa_hifi/my_pipe2/step1_mapping/reads.bed"
# file_in = "/public/home/hpc214712170/Test/tests/sativa_hifi/my_pipe2/step1_mapping/CP132241.1:18130284-18201609/pri.pos"
# file_out = "/public/home/hpc214712170/Test/tests/sativa_hifi/my_pipe2/step1_mapping/CP132241.1:18130284-18201609/pri.bed"
file_in = "/public/home/hpc214712170/Test/tests/sativa_hifi/my_pipe2/step1_mapping/CP132245.1:6957842-6985098/reads.pos"
file_out = "/public/home/hpc214712170/Test/tests/sativa_hifi/my_pipe2/step1_mapping/CP132245.1:6957842-6985098/reads.bed"
def cluster(arr, threshold=100):
    clusters = []
    # 排序数组
    arr = sorted(arr)
    # 初始化第一个簇
    current_cluster = [arr[0]]
    # 从第二个元素开始遍历
    for i in range(1, len(arr)):
        # 如果当前元素与簇中最后一个元素的差异小于阈值，则加入当前簇
        if abs(arr[i] - current_cluster[-1]) <= threshold:
            current_cluster.append(arr[i])
        else:
            # 否则，当前簇结束，将其加入到聚类结果中，并开始一个新的簇
            clusters.append(current_cluster)
            current_cluster = [arr[i]]
    # 添加最后一个簇
    clusters.append(current_cluster)
    return clusters
    
with open(file_in, "r") as fin, open(file_out, "w") as fo:
    dic = defaultdict(list)
    for line in fin:
        fields = line.strip().split("\t")
        # fo.write("{}\t{}\t{}\n".format(fields[0], fields[1], int(fields[1]) + 100))
        ctg, pos = fields[0], fields[1]
        pos = int(pos)
        dic[ctg].append(pos)
    # print(dic)
    for ctg, ls in dic.items():
        print(ctg)
        clusters = cluster(ls, 3000)
        for arr in clusters:
            # if len(arr)
            print(arr)
            if len(arr) > 5:
                fo.write("{}\t{}\t{}\t{}\n".format(ctg, arr[0], arr[-1], len(arr)))
        print(len(clusters))
# print(prefix)
# ls1 = [str(i) + prefix for i in range(10000000)]
# ls1 = [str(i) for i in range(10000000)]
# t0 = time.time()
# ls = [i + "\n" for i in ls1]
# file = "/public/home/hpc214712170/Test/tests/sativa_hifi/test/write.txt"
# t1 = time.time()
# with open(file, "w") as f:
#     f.writelines(ls)
#     # for i in ls1:
#     #     f.write("{}\n".format(i))
# print("cost1:{}, cost2:{}".format(t1 - t0, time.time() - t1))
# dic = defaultdict(set)
# keys = [str(i) for i in range(200)]
# for i in range(200):
#     s = set()
#     for j in range(100000):
#         id = i * 10000 + j
#         s.add(id)
#     dic[keys[i]] = s
# s = set()
# t0 = time.time()
# for key, val in dic.items():
#     s.update(val)
# print(time.time() - t0)

# s = set()
# for i in range(10):
#     s.add(i)
# print(" ".join(s))
# fq=../simu/ont/fq/merge.fastq
# asm=../simu/Ref/template_simu.fasta
# bam=../simu/ont/aln2simu.sort.bam
# ref=../simu/Ref/template_ref.fasta

# ../simu/ont/hifi/merge.fastq
# ../simu/hifi/aln2simu.sort.bam