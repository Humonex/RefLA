import numpy as np
import sys
import os
import pandas as pd


def stats(bed_file):
    with open(bed_file, "r")as f:
        ls = []
        for line in f:
            if not line or line.startswith("#"): continue
            ctg, start, end = line.strip().split("\t")[:3]
            start, end = int(start), int(end)
            ls.append(end - start)
        
        ## stats
        print(sorted(ls))
        print("reg num: ", len(ls))
        print("mean: ", np.mean(ls))
        print("median: ", np.median(ls))
        step = 0.01
        quantile_ls = [i * step for i in range(int(1 / step))]
        quantile_num_ls = np.quantile(ls, quantile_ls)
        # if reverse_flag:
        #     quantile_num_ls = np.flip(quantile_num_ls)
        for i in range(len(quantile_ls)):
            print("{:.3f}:{:.0f}".format(quantile_ls[i], quantile_num_ls[i]))

def stats_dic(dic, reverse_flag_dic, out_file):
    fo = open(out_file, "w")
    step = 0.001
    quantile_ls = [i * step for i in range(int(1 / step))]
    final_out = [[] for i in range(len(quantile_ls))]
    key_ls= [] 
    for key, ls in dic.items():
        key_ls.append(key)
        new_ls = [a for a in ls if a is not None]
        # print("# Reg num: {}, Error_portion{:.4f}".format(len(new_ls), 1 - len(new_ls)))
        # print("# mean: {}, median: {}".format(np.mean(new_ls), np.median(new_ls)))
        fo.write("# --------{}--------\n".format(key))
        fo.write("# Normal reg num: {}, Error reg portion: {:.4f}\n".format(len(new_ls), 1 - len(new_ls)/len(ls)))
        fo.write("# mean: {}, median: {}\n".format(np.mean(new_ls), np.median(new_ls)))
        ### 
        # print(new_ls)
        quantile_num_ls = np.quantile(new_ls, quantile_ls)
        if reverse_flag_dic[key]:
            quantile_num_ls = np.flip(quantile_num_ls)
        # ressults = []
        for i in range(len(quantile_ls)):
            # print("{:.3f}:{:.4f}".format(quantile_ls[i], quantile_num_ls[i]))
            # ressults.append(quantile_num_ls[i])
            final_out[i].append(quantile_num_ls[i])
    # print(final_out)
    ## 
    fo.write("#threashold\t" + "\t".join(key_ls) + "\n")
    for idx, ls in enumerate(final_out):
        fo.write("{:.3f}\t".format(quantile_ls[idx]))
        for a in ls:
            fo.write("{:.4f}\t".format(a))
        fo.write("\n")
    fo.close()
def add_ls(ls, a):
    if a:
        ls.append(float(a))
    else:
        ls.append(None)


def merge_pileup_feature(data_dir):
    print("--------------------Start stats for {}--------------------".format(os.path.abspath(data_dir)))
    head = ["contig", "start_pos", "end_pos", "correct_portion", "disagree_portion", "differ_portion"]
    part_dir = data_dir + "/parts"
    all_files = os.listdir(part_dir)
    # print(all_files)
    ls1 = []
    ls2 = []
    ls3 = []
    for file in all_files:
        if not file.endswith("feature.txt"): continue
        file_path = os.path.join(part_dir, file)
        print("Collect", file_path)
        with open(file_path, "r") as f:
            for line in f:
                if line.startswith("#") or line.startswith("contig\t"): continue
                # print(line)
                fields = line.strip("\n").split("\t")
                # break
                ## 
                add_ls(ls1, fields[3])
                add_ls(ls2, fields[4])
                add_ls(ls3, fields[5])
    ## correct_portion", "disagree_portion", "differ_portion
    dic = {"correct_portion":ls1, "disagree_portion":ls2, "differ_portion":ls3}
    reverse_flag_dic = {"correct_portion":False, "disagree_portion":True, "differ_portion":True}
    out_file = data_dir + "/pileup.stats"
    stats_dic(dic, reverse_flag_dic, out_file)

    # print("\ncorrect_portion")
    # stats_ls(ls1, False)
    # print("\ndisagree_portion")
    # stats_ls(ls2, True)
    # print("\ndiffer_portion")
    # stats_ls(ls3, True)

if __name__ == "__main__":

    ## run to get
    '''
    sys.path.append("/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err")
    from find_candidate_regions import find_from_pileup
    threads = 40
    ctg_ls = ['chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII', 'chrIX', 'chrX', 'chrXI', 'chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI', 'chrMT']
    ref = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/yeast_ont/my_pipe3/corrected_ref/reference.fasta"
    bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/yeast_ont/my_pipe3/step1_mapping/aln.sorted.bam"
    out_dir = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/yeast_ont/my_pipe3/step2_candidate_regions/pileup2"
    # params = {'reg_size': 1000000, 'win_size': 400, 'step_size': 200, \
    #     'min_correct_portion': 0.9, 'max_differ_portion': 0.1, \
    #     'max_disagree_portion': 0.02, 'cluster_dis': 1000}
    params = {'reg_size': 1000000, 'win_size': 400, 'step_size': 200, \
        'min_correct_portion': 0.9, 'max_differ_portion': 0.1, \
        'max_disagree_portion': 0.03, 'cluster_dis': 1000}
    find_from_pileup.parse_pileup_parallel(ctg_ls, ref, bam_in, threads, params, out_dir)'''
    
    ## merge and stats
    out_dir = sys.argv[1]
    merge_pileup_feature(out_dir)

    ## 
    '''
    yeast: params = {'reg_size': 1000000, 'win_size': 400, 'step_size': 200, \
        'min_correct_portion': 0.9, 'max_differ_portion': 0.15, \
        'max_disagree_portion': 0.03, 'cluster_dis': 1000}
    '''
    # bed_file = sys.argv[1]
    # stats(bed_file)
