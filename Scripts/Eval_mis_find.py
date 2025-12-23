
import os
import sys
from collections import defaultdict

'''
precision=TP/(TP+FP)
recall=TP/(TP+FN) //评估所有实际正例是否被预测出来的覆盖率占比多少
F1=2*precision*recall/(precision + recall)
'''


    
# For test
# bed1 = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast/simu/hack2/h1_new.bed"
# bed2 = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast/inspector/structural_error.bed"
# bed2 = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast/CRAQ/craq_out/runAQI_out/locER_out/out_final.CRE.bed"
# bed2 = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast/CRAQ/craq_out/runAQI_out/locER_out.bed"
# bed2 = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast/CRAQ/craq_out/runAQI_out/strER_out.bed"
# bed2 = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast/CRAQ/craq_out/runAQI_out/test.bed"
# bed2 = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast/my_pipe/step2_candidate_regions/candidate.bed"
# bed2 = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast/quast/mis.bed"
# bed2 = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast/quast/1.bed"
# bed2 = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast/quast/icarus_viewers/mis.bed" 
# failed_bed = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast/ins_failed.bed"

class Bed:
    def __init__(self) -> None:
        pass
    @staticmethod
    def read_bed(bed):
        with open(bed, "r") as fin:
            dic = defaultdict(str)
            for line in fin:
                if line.startswith("#"):continue
                fileds = line.strip().split("\t")
                chr, start, end = fileds[:3]
                if ";" in start:    # inspector HaplotypeSwitch
                    start_ls = start.split(";")
                    end_ls = end.split(";")
                    start = min(int(start_ls[0]), int(start_ls[1]))
                    end = max(int(end_ls[0]), int(end_ls[1]))
                else:
                    start, end = int(start), int(end)
                if chr in dic.keys():
                    dic[chr].append((start, end))
                else:
                    dic[chr] = [(start, end)]
            # print("chrs in bed: ", [chr for chr in dic.keys()], os.path.abspath(bed))
            sorted_dic = defaultdict(str)
            for chr, reg_ls in dic.items():
                reg_ls.sort()
                sorted_dic[chr] = reg_ls
        return sorted_dic
        # pass
    @staticmethod
    def write_bed(bed, reg_ls):
        with open(bed, "w") as fo:
            for i in reg_ls:
                fo.write("\t".join(i) + "\n")

def call(truth_reg, call_reg, threshold=0.4, bias=0):
    # print("threshold:", threshold, "bias:", bias)
    # threshold = 0.4
    # bias = 0
    call_reg = [max(call_reg[0]-bias, 0), call_reg[1]+bias]
    # if 
    inter = min(truth_reg[1], call_reg[1]) - max(truth_reg[0], call_reg[0])
    # threshold = 0.01
    # print(truth_reg)
    # print("inter:", inter)
    # print("len(truth_reg): ", truth_reg[1] - truth_reg[0])
    # if inter / len(truth_reg) > 0.5:return True
    # if truth_reg[1] == truth_reg[0]:
    #     print(truth_reg)
    if truth_reg[1] - truth_reg[0] <= 0:
        print("Error:", truth_reg)
        raise ValueError
    if inter / (truth_reg[1] - truth_reg[0]) > threshold:return True
    return False
    # pass

def stats_no_zero(ls):
    num = 0
    for i in ls:
        if i > 0:
            num += 1
    return num
def stats_zero(ls):
    num = 0
    for i in ls:
        if i == 0:
            num += 1
    return num
def eval(true_reg_ls, call_reg_ls, chr, threshold=0, bias=0):
    # print("threshold:", threshold, "bias:", bias)
    '''
    eval single chr
    可能还需要改进
    '''
    recall_ls = [0] * len(true_reg_ls)
    precise_ls = [0] * len(call_reg_ls)
    succes_ls = []
    '''i, j = 0, 0
    while i < len(true_reg_ls) and j < len(call_reg_ls):
        while i < len(true_reg_ls) and j < len(call_reg_ls) and true_reg_ls[i][1] < call_reg_ls[j][0]: i+=1
        while i < len(true_reg_ls) and j < len(call_reg_ls) and call_reg_ls[j][1] < true_reg_ls[i][0]: j+=1
        if i >= len(true_reg_ls) or j >= len(call_reg_ls): break
        ## call
        if call(true_reg_ls[i], call_reg_ls[j]):
            recall_ls[i] += 1
            precise_ls[j] += 1
        ## 
        if true_reg_ls[i][1] >= call_reg_ls[j][1]:
            j += 1
        else:
            i += 1'''
            
    # for recall ls cal
    j = 0
    for i in range(len(true_reg_ls)):
        while j < len(call_reg_ls) and call_reg_ls[j][1] < true_reg_ls[i][0] - bias: j+=1
        if i >= len(true_reg_ls) or j >= len(call_reg_ls): break
        if call(true_reg_ls[i], call_reg_ls[j], threshold, bias):
            recall_ls[i] += 1
            succes_ls.append([chr, true_reg_ls[i][0], true_reg_ls[i][1]])
    # for precise ls cal
    i = 0
    for j in range(len(call_reg_ls)):
        while i < len(true_reg_ls) and true_reg_ls[i][1] < call_reg_ls[j][0] - bias: i+=1
        if i >= len(true_reg_ls) or j >= len(call_reg_ls): break
        if call(true_reg_ls[i], call_reg_ls[j], threshold, bias):
            precise_ls[j] += 1
    ## stats
    # print("Recall:{}".format(recall_ls))
    # print("Precise:{}".format(precise_ls))
    recall_num = 0
    for i in recall_ls:
        if i > 0: recall_num += 1
    # print("Recall rate:{}".format(recall_num / len(recall_ls)))
    # print("Precise_rate:{}".format(sum(precise_ls) / len(precise_ls)))
    # if len(recall_ls) > 0:
    #     print("Recall rate:{}".format(stats_no_zero(recall_ls) / len(recall_ls)))
    # else:
    #     print("Recall rate:None")
    # if len(precise_ls) > 0:
    #     print("Precise rate:{}".format(stats_no_zero(precise_ls) / len(precise_ls)))
    # else:
    #     print("Precise rate:None")
    return recall_ls, precise_ls, succes_ls
    # pass

def cal_f1(precision, recall):
    # if precision
    if precision + recall != 0:
        return 2*precision*recall/(precision + recall)
    else:
        return 0
def cal_jaccarrd2(comm1, comm2, a, b):
    '''
    将comm1和a放缩为comm2和a_
    '''
    if comm1 == 0 or comm2 == 0: return 0
    a_ = (comm2 / comm1) * a
    return comm2/(comm2 + a_ + b)

def evalAll(bed1, bed2, failed_bed, false_bed, ex_ls, threshold=0, bias=0):
    fo = open(bed2+".eval", "w")
    dic1 = Bed.read_bed(bed1)   # truth bed
    dic2 = Bed.read_bed(bed2)   # calling bed for eval
    chr_ls = set()
    for chr in dic1.keys():
        chr_ls.add(chr)
    for chr in dic2.keys():
        chr_ls.add(chr)
    # print(chr_ls)
    recall_ls = []
    precise_ls = []
    failed_recall_ls = []
    False_pre_ls = []
    succes_ls = []
    fo.write("Truth bed:{}\n".format(os.path.abspath(bed1)))
    fo.write("threshold:{}, bias:{}\n".format(threshold, bias))
    fo.write("Ex_ls:{}\n".format(ex_ls))
    for chr in chr_ls:
        if chr in ex_ls: 
            # print("Skip:{}".format(chr))
            continue
        # print("\n" + chr + ":")
        # print("Truth:", dic1[chr], ", Call:", dic2[chr])
        print(dic1[chr], dic2[chr])
        chr_recall_ls, chr_precise_ls, chr_succes_ls = eval(dic1[chr], dic2[chr], chr, threshold, bias)
        recall_ls.extend(chr_recall_ls)
        precise_ls.extend(chr_precise_ls)
        succes_ls.extend(chr_succes_ls)
        for i in range(len(chr_recall_ls)):
            if chr_recall_ls[i] == 0:
                failed_recall_ls.append([chr, str(dic1[chr][i][0]), str(dic1[chr][i][1])])
        for i in range(len(chr_precise_ls)):
            if chr_precise_ls[i] == 0:
                False_pre_ls.append([chr, str(dic2[chr][i][0]), str(dic2[chr][i][1])])
        # break
        # pass
    fo.write("Failed recall:{}\n".format(failed_recall_ls))
    fo.write("False_pre_ls:{}\n".format(False_pre_ls))
    fo.write("Succes ls:{}\n".format(succes_ls))
    # print("Failed recall:{}".format(failed_recall_ls))
    # print("False_pre_ls:{}".format(False_pre_ls))
    # print("Succes ls:{}".format(succes_ls))
    # print("-----------------------------Sum stats-----------------------------")
    # print("Truth set:{}".format(len(recall_ls)))
    # print("Predict set:{}".format(len(precise_ls)))
    TP1 = stats_no_zero(recall_ls)
    TP2 = stats_no_zero(precise_ls)
    FP = stats_zero(precise_ls)
    FN = stats_zero(recall_ls)
    # print("TP1:{}, TP2:{}".format(TP1, TP2))
    # print("FP:{}".format(FP))
    # print("FN:{}".format(FN))
    fo.write("Truth:{}, Call:{}\n".format(len(recall_ls), len(precise_ls)))
    fo.write("TP1:{}, TP2:{}\n".format(TP1, TP2))
    fo.write("FP:{}, FN:{}\n".format(FP, FN))
    # fo.write("FN:{}\n".format(FN))q
    Precise_val = stats_no_zero(precise_ls) / len(precise_ls)   # predict 
    Recall_val = stats_no_zero(recall_ls) / len(recall_ls)
    print(Precise_val, Recall_val)
    print(precise_ls, recall_ls)
    F1_val = cal_f1(Precise_val, Recall_val)
    # print("Recall rate:{}".format(Recall_val))
    # print("Precise rate:{}".format(Precise_val))
    # print("F1 score:{}".format(F1_val))
    fo.write("Recall rate:{}\n".format(Recall_val))
    fo.write("Precise rate:{}\n".format(Precise_val))
    fo.write("F1 score:{}\n".format(F1_val))
    # print("Recall rate:{}".format(sum(recall_ls) / len(recall_ls)))
    # print("Precise_rate:{}".format(sum(precise_ls) / len(precise_ls)))
    # failed_bed = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast/ins_failed.bed"
    Bed.write_bed(failed_bed, failed_recall_ls)
    Bed.write_bed(false_bed, False_pre_ls)
    mis_num = len(recall_ls)
    call_num = len(precise_ls)
    # 
    jacc = cal_jaccarrd2(TP2, TP1, FP, FN)  # TP1: QUAST TP, TP2: LRMD TP, FP: LRMD only, FN: QUAST only
    # jacc = cal_jaccarrd2(TP1, TP2, FN, FP)    # 
    print("jacc: ", jacc)
    fo.write("jacc:{}\n".format(jacc))

    fo.close()
    return Recall_val, Precise_val, F1_val, mis_num, call_num

if __name__ == "__main__":
    ## 
    # python /public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/Eval__mis_SV/Eval_mis_find.py ../../../../../quast_eval_large/contigs_reports/mis_cl.bed merge.bed 0 0 ../../../../../quast_eval_large/contigs_reports/ex.txt
    bed1 = sys.argv[1]  # truth bed
    bed2 = sys.argv[2]   # calling bed for eval
    # failed_bed = sys.argv[3] # failed bed
    # false_bed = sys.argv[4]
    false_bed = bed2 + ".false.bed"
    failed_bed = bed2 + ".failed.bed"
    threshold = float(sys.argv[3])
    bias = int(sys.argv[4])
    ex_f = sys.argv[5]
    ex_ls = []
    with open(ex_f, "r") as f:
        for line in f:
            if line.startswith("#"):continue
            if line:
                ex_ls.append(line.strip())
    # print("ex_ls:", ex_ls)
    print("threshold:", threshold, "bias:", bias)
    Recall_val, Precise_val, F1_val, mis_num, call_num = evalAll(bed1, bed2, failed_bed, false_bed, ex_ls, threshold, bias)
    print("{}|{}|{}|{}".format(Recall_val, Precise_val, F1_val, call_num))
    exit()
    #################################################   For single simu data
    ex_ls = []
    bed1 = sys.argv[1]  # truth bed
    bed2 = sys.argv[2]   # calling bed for eval
    # failed_bed = sys.argv[3] # failed bed
    # false_bed = sys.argv[4]
    false_bed = bed2 + ".false.bed"
    failed_bed = bed2 + ".failed.bed"
    threshold = float(sys.argv[3])
    bias = int(sys.argv[4])
    print("threshold:", threshold, "bias:", bias)
    Recall_val, Precise_val, F1_val, mis_num, call_num = evalAll(bed1, bed2, failed_bed, false_bed, ex_ls, threshold, bias)
    print("{}|{}|{}|{}".format(Recall_val, Precise_val, F1_val, call_num))
    exit()
    #################################################   For simu data parallel
    ex_ls = []
    # sample = "sativa"
    # sample = "chm13"
    sample = "GRCH38"
    # data_type = "ont"
    # data_type = "hifi"
    data_type = "clr"
    res1 = "filtered/merge/merge.bed"
    res2 = "filtered2/merge/merge.bed"
    quast1 = "quast"
    quast2 = "quast_large"
    threshold = 0
    bias = 2000
    bed1 = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/" + sample + "/simu/Ref/mis_simu.bed"
    work_dir = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/" + sample + "/" + data_type
    # dic = {"gaep":"gaep/mis_breakpoints.txt", "craq":"craq_out/runAQI_out/mis.bed", "inspector":"inspector/structural_error.bed", "quast":"quast_large/contigs_reports/mis.bed", "my_pipe":"my_pipe/step2_candidate_regions/filtered2/merge.bed"}   # filtered | filtered2
    dic = {"gaep":"gaep/mis_breakpoints.txt", "craq":"craq_out/runAQI_out/mis.bed", "inspector":"inspector/structural_error.bed", "quast":"quast/contigs_reports/mis.bed", "my_pipe":"my_pipe/step2_candidate_regions/" + res2}   # filtered | filtered2
    tool_ls = ["gaep", "craq", "inspector", "quast", "my_pipe"]
    print("threshold:{}, bias:{}".format(threshold, bias))
    print("Truth:{}".format(bed1))
    print("Sample:{}, type:{}".format(sample, data_type))
    with open("/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/Eval__mis_SV/eval/" + sample + "_" + data_type + "_stats", "w") as f:
        for tool in tool_ls:
            print("---------------------{}---------------------".format(tool))
            # f.write("---------------------{}---------------------\n".format(tool))
            bed2 = work_dir + "/" + dic[tool]
            failed_bed = work_dir + "/" + dic[tool] + ".failed.bed"
            false_bed = work_dir + "/" + dic[tool] + ".false.bed"
            if os.path.isfile(bed2):
                Recall_val, Precise_val, F1_val, mis_num, call_num = evalAll(bed1, bed2, failed_bed, false_bed, ex_ls, threshold, bias)
                print("{}|{}|{}|{}".format(Recall_val, Precise_val, F1_val, call_num))
                # f.write("{}|{}|{}\n".format(Recall_val, Precise_val, F1_val))
            else:
                print("{} is none".format(bed2))
    pass
