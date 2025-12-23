import yaml
from create_consensus_by_bed import fasta_parser
from create_consensus_by_bed.Utils import Connect_info
# import fasta_parser
# from Utils import Connect_info


def get_scaffold(new_seq, scaffold_id, fa_out_dic, scaffold_ls):
    if len(new_seq) > 0:
        fa_out_dic[scaffold_id] = new_seq  # 记录scaffold序列
        if len(scaffold_ls) > 1:
            print("Scaffolding from {} -> {}".format(",".join(scaffold_ls), scaffold_id))
        else:
            print("Keep from {} -> {}".format(scaffold_ls[0], scaffold_id))
    else:
        print("no seq")

def scaffold(fa_in, connect_file, fa_out, config):
    fa_in_dic = fasta_parser.read_sequence_dict(fa_in)
    print("Read fasta in done!!!")
    fa_out_dic = {}
    connect_info_ls = Connect_info.read_connect_info(connect_file)
    print("Read connect file done!!!")
    N_fill_size = config["Scaffold_Params"]["N_fill_size"]
    N_fill_seq = N_fill_size * "N"
    ## 性能参数
    MAX_GAP = config["Scaffold_Params"]["MAX_GAP"]  # 1000000
    MIN_CONNECT_SEQ = config["Scaffold_Params"]["MIN_CONNECT_SEQ"]  # 20000000 超过这个距离就不再scaffold
    for info in connect_info_ls:
        if len(info.connect_ls) == 0:
            print("{} has been discarded".format(info.chr_id))
            continue
        if len(info.connect_ls) == 1:   # 一条完整的seq，不用scaffolding
            ## 有种可能是原序列不存在（抛光过程中丢弃了那条序列）
            if fa_in_dic.get(info.connect_ls[0], ".") != ".":            
                fa_out_dic[info.chr_id] = fa_in_dic[info.connect_ls[0]]    # 保留原序列
                print("{} keep whole contig".format(info.chr_id))
            else:print("Missed {} in fasta".format(info.chr_id))
            continue
        ## 
        new_seq = ""
        try:
            new_seq += fa_in_dic[info.connect_ls[0]]    # 加上首段seq
        except:
            raise Exception("Scaffold序列不存在")  ## 有种可能是原序列不存在（抛光过程中丢弃了那条序列）
        scaffold_ls = [info.connect_ls[0]]
        cnt = 0
        for idx, seq_id in enumerate(info.connect_ls[1:]):  # 从第二个seq开始
            # 查询第i-1段与第i段之间的gap情况。注意上面是
            gap_len, gap_type = info.gap_ls[idx].split("_")
            gap_len = abs(int(gap_len))
            now_seq = fa_in_dic[seq_id] # 当前seq的序列
            if gap_type == "INV": # INV：不进行scaffolding， cut off
                scaffold_id = info.chr_id + "_" + str(cnt)
                get_scaffold(new_seq, scaffold_id, fa_out_dic, scaffold_ls)
                ## 截断后的一些重置
                cnt += 1    # 
                scaffold_ls = [seq_id]  # 重置
                new_seq = now_seq
            elif gap_type == "GAP":
                ## 
                if gap_len < MAX_GAP and len(new_seq) > MIN_CONNECT_SEQ and len(now_seq) > MIN_CONNECT_SEQ: # scaffolding
                    new_seq += N_fill_seq + now_seq
                    scaffold_ls.append(seq_id)
                    # new_seq += now_seq
                else:   # cut off
                    scaffold_id = info.chr_id + "_" + str(cnt)
                    get_scaffold(new_seq, scaffold_id, fa_out_dic, scaffold_ls)
                    ## 截断后的一些重置
                    cnt += 1    # 
                    scaffold_ls = [seq_id]  # 重置
                    new_seq = now_seq
            else:
                raise ValueError
        ## 
        scaffold_id = info.chr_id + "_" + str(cnt)
        get_scaffold(new_seq, scaffold_id, fa_out_dic, scaffold_ls)
    print("write fasta_dict")
    fasta_parser.write_fasta_dict(fa_out_dic, fa_out)
      
if __name__ == "__main__":
    # fa = "/public/home/hpc214712170/Test/tests/chm13_hifi_ctgs/NC_18/my_pipe/step3_SV_consensus/fa/merge.fasta"
    # fa_out = "/public/home/hpc214712170/Test/tests/chm13_hifi_ctgs/NC_18/my_pipe/step3_SV_consensus/merge.fasta"
    # fa_dic = fasta_parser.read_sequence_dict(fa)
    # new_fa_dic = {}
    # for key, val in fa_dic.items():
    #     N_pos = [i for i,c in enumerate(val) if c == "N"]
    #     # N_pos = re.findall("N", val)
    #     # print(N_pos)
    #     start, end = -1, -1
    #     ls = []
    #     for pos in N_pos:
    #         if start >= 0:
    #             if pos - end < 2:
    #                 end = pos
    #             else:
    #                 ls.append([start, end])
    #                 start, end = pos, pos + 1
    #         else:
    #             start, end = pos, pos + 1
    #     if start > -1:
    #         ls.append([start, end])
    #         print(ls)
    #         seq_ls = []
    #         pre_end = 0
    #         for start, end in ls:
    #             seq_ls.append(val[pre_end:start])
    #             print("seq from {}-{}".format(pre_end, start))
    #             pre_end = end
    #         seq_ls.append(val[pre_end:len(val)])
    #         print("seq from {}-{}".format(pre_end, len(val)))
    #         # print(len(seq_ls))
    #         for i, seq in enumerate(seq_ls):
    #             new_chr = key + "_" + str(i)
    #             new_fa_dic[new_chr] = seq
    #     else:
    #         new_fa_dic[key] = val
    #         pass
    
    # fasta_parser.write_fasta_dict(new_fa_dic, fa_out)
    
    ## statstic of GAP len
    # import numpy as np
    # import math
    # file = "/public/home/hpc214712170/Test/tests/chm13_hifi/my_pipe/step3_SV_consensus/connect.txt"
    # ls = apply_rec.Connect_info.read_connect_info(file)
    # len_ls = []
    # for gap in ls.gap_ls:
    #     gap_len = int(gap.split(",")[0])
    #     len_ls.append(gap_len)
    # len_ls.sort()
    # print(len_ls)
    # print(np.average(len_ls))

    ## 
    # import time
    # t1 = time.time()
    # fa_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi_ctgs/NC_18/my_pipe/step4_polish/racon.fasta"
    # connect_file = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi_ctgs/NC_18/my_pipe/step3_SV_consensus/connect.txt"
    # fa_out = "/public/home/hpc214712170/Test/tests/chm13_hifi_ctgs/NC_18/my_pipe/step5_scaffold/final.fasta"

    # # fa_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi/my_pipe/step4_polish/racon.fasta"
    # # connect_file = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi/my_pipe/step3_SV_consensus/connect.txt"
    # # fa_out = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi/my_pipe/step5_scaffold/final.fasta"

    # N_fill_size = 20
    # scaffold(fa_in, connect_file, fa_out, N_fill_size)
    # t2 = time.time()

    work_dir = "/public/home/hpc214712170/Test/tests/bench/reinhardtii_hifi/my_pipe2"
    yaml_file = "/public/home/hpc214712170/Test/tests/bench/reinhardtii_hifi/Config.yaml"
    with open(yaml_file, "r") as f:
        config = yaml.safe_load(f.read())
        print(type(config), config)
    fa_in = work_dir + "/step4_polish/racon.fasta"
    connect_file =  work_dir + "/step3_SV_consensus/connect.txt"
    fa_out = work_dir + "/step5_scaffolding/final.fasta"
    scaffold(fa_in, connect_file, fa_out, config["step5"])
