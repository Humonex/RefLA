import pysam
import re
import numpy as np
from scipy.cluster.hierarchy import linkage     #导入linage函数用于层次聚类
from scipy.cluster.hierarchy import dendrogram  #dendrogram函数用于将聚类结果绘制成树状图
from scipy.cluster.hierarchy import fcluster    #fcluster函数用于提取出聚类的结果
from sklearn.datasets import make_blobs         #make_blobs用于生成聚类算法的测试数据
from sklearn.cluster import AgglomerativeClustering  #自底向上层次聚类算法
import matplotlib.pyplot as plt                 #导入matplotlib绘图工具包

def cal_ident(read):
    ciger = read.ci
def test():
    bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi_ctgs/NC_19/my_pipe/step1_mapping/aln.sorted.bam"
    ctg = "NC_000019.10"
    start = 27369501
    end = 27371000
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")

    min_support_reads = 3
    MIN_MAPPING_QUALITY = 60
    MIN_ALIGN_LENGTH = 10000
    MIN_ALIGN_RATE = 0.95
    # MIN_CLIP_LEN = 1000
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    support_reads_ls = []
    for read in bam_reader.fetch(ctg, start, end):
        ## filter1
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < MIN_MAPPING_QUALITY:
            continue    # 只保留primary
        ## filter2
        if read.query_alignment_length < MIN_ALIGN_LENGTH or (read.query_alignment_length / read.query_length) < MIN_ALIGN_RATE:
            continue
        ## 
        # cigar = read.cigarstring
        # tokens = re.findall("[\d]{0,}[A-Z]{1}", cigar)
        # left, right = tokens[0], tokens[-1]
        # if left[-1] in "HS" and int(left[:-1]) > MIN_CLIP_LEN: continue
        # if right[-1] in "HS" and int(right[:-1]) > MIN_CLIP_LEN: continue
        if read.reference_start > start-1000 or read.reference_end < end + 1000:continue
        support_reads_ls.append(read)
        # if read.reference_start < start-1000 and read.reference_end > end + 1000:
        #     support_reads_ls.append(read)
    iden_ls = []
    read_sim_ls = []
    for read in support_reads_ls:
        print("{}:{}-{}".format(read.query_name, read.reference_start, read.reference_end))
        print(read.query_alignment_length, read.query_alignment_length / read.query_length)
        print(read.query_name, read.get_tag("NM"))
        print(1-read.get_tag("NM")/read.query_length)
        iden_ls.append(1-read.get_tag("NM")/read.query_length)
        read_sim_ls.append([read, 1-read.get_tag("NM")/read.query_length])
        print("\n")
    print(iden_ls)
    print(len(iden_ls))
    if len(support_reads_ls) < min_support_reads:
        return False
    read_sim_ls.sort(key=lambda x:x[1], reverse=True)
    # print(sorted(iden_ls))
    # print(iden_ls)
    from sklearn.cluster import KMeans
    # 实例化KMeans,聚成3类 
    kmeans = KMeans(n_clusters=3, n_init='auto')
    data = np.array([x[1] for x in read_sim_ls])
    # 训练模型
    kmeans.fit(data.reshape(-1,1))
    # 获得聚类标签
    labels = kmeans.labels_ 
    print(labels)
    for i in range(len(data)):
        print("{}:{}\t->{}".format(read_sim_ls[i][0].query_name, data[i], labels[i]))
    
    high_ident_lable = labels[0]
    high_ident_nums = np.bincount(labels)[high_ident_lable]
    print(high_ident_nums)
    cnt  = 0
    for i in labels:
        if i == high_ident_lable:
            cnt += 1
    print(cnt)

test()

def test_cluster():
    ##
    from sklearn.cluster import KMeans
    import time
    t1 = time.time()
    # Reads相似度列表
    # sim = np.array([[0.5807367449783519, 0.5871770042301616, 0.5893947806774014, 0.5939368942049907, 0.5971657346212507, 0.5980131999727836, 0.62181240063593, 0.6238909069255048, 0.6379091294089176, 0.6441538277797683, 0.6503691645834544, 0.6536404322672198, 0.654084743980654, 0.6566518966821571, 0.6607922214242308, 0.6682704744166796, 0.6683333333333333, 0.6684922932477897, 0.6734932307057935, 0.6772517939380958, 0.6796005564387917, 0.6812327737409171, 0.6819558203028047, 0.6836982968369829, 0.6842779694430754, 0.6880889925191878, 0.6932061514502628, 0.69775390625, 0.7065608068236682, 0.7078007558707375, 0.7098902157164869, 0.7117698534237584, 0.7124120461581762, 0.7217032094057834, 0.7269936938445225, 0.7283240319889535, 0.7291675438965851, 0.7317216981132075, 0.7403766691313993, 0.7417929292929293, 0.7459517812162648, 0.7462729221243153, 0.7517407854423916, 0.7713457899716178, 0.7732012267044114, 0.7908961593172119, 0.7941254651780968, 0.8168890559418656, 0.8388119057975517, 0.8514050987807263, 0.9934805237495207, 0.9939124326855537, 0.9946998123827392, 0.9948219688154526, 0.9952253528516711, 0.9957546479285609, 0.9963770502601936, 0.996535648603567, 0.9968624497991968, 0.9976722532588455, 0.9977456223130964, 0.9980146537461593, 0.9980969372584002, 0.9983292891523131, 0.9984218832193582, 0.9984827137747914, 0.9985391223046807, 0.9986757565353868, 0.9990820238673794, 0.9990853658536586, 0.9991051721052318, 0.9991682015328858, 0.9993195735994557, 0.9993461422392114, 0.9993986770895971, 0.9994019964119785, 0.9995701037892281, 0.9996543858436442, 0.9996681048788583, 0.9997902537929105, 0.9998707425838558]])
    sim = np.array([0.7908961593172119, 0.9991682015328858, 0.9986757565353868, 0.9957546479285609, 0.9968624497991968, 0.9980146537461593, 0.7517407854423916, 0.9997902537929105, 0.7462729221243153, 0.9990820238673794, 0.9984827137747914, 0.69775390625, 0.7117698534237584, 0.9939124326855537, 0.7065608068236682, 0.9980969372584002, 0.6796005564387917, 0.6880889925191878, 0.7459517812162648, 0.9996543858436442, 0.7317216981132075, 0.9952253528516711, 0.9946998123827392, 0.9984218832193582, 0.9948219688154526, 0.7124120461581762, 0.9991051721052318, 0.9983292891523131, 0.7078007558707375, 0.6842779694430754, 0.7269936938445225, 0.6932061514502628, 0.9996681048788583, 0.9998707425838558, 0.9995701037892281, 0.6683333333333333, 0.9990853658536586, 0.6819558203028047, 0.7291675438965851, 0.9993461422392114, 0.6734932307057935, 0.6682704744166796, 0.6536404322672198, 0.996535648603567, 0.9934805237495207, 0.6836982968369829, 0.6812327737409171, 0.5971657346212507, 0.6607922214242308, 0.6503691645834544, 0.6379091294089176, 0.6772517939380958, 0.7098902157164869, 0.9993986770895971, 0.9993195735994557, 0.6566518966821571, 0.8388119057975517, 0.5980131999727836, 0.9985391223046807, 0.5939368942049907, 0.62181240063593, 0.6238909069255048, 0.9963770502601936, 0.6684922932477897, 0.7417929292929293, 0.7941254651780968, 0.7732012267044114, 0.8168890559418656, 0.7713457899716178, 0.8514050987807263, 0.5807367449783519, 0.5893947806774014, 0.654084743980654, 0.7217032094057834, 0.9976722532588455, 0.9994019964119785, 0.6441538277797683, 0.7403766691313993, 0.7283240319889535, 0.9977456223130964, 0.5871770042301616])
    # 实例化KMeans,聚成3类 
    kmeans = KMeans(n_clusters=3, n_init='auto')

    # 训练模型
    kmeans.fit(sim.reshape(-1,1))
    # 获得聚类标签
    labels = kmeans.labels_ 
    print(labels)
    for i in range(len(sim)):
        print("{}:{}".format(sim[i], labels[i]))
    print(time.time() - t1)
# test_cluster()
