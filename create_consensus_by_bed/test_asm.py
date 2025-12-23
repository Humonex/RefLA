from collections import defaultdict
def check_asm_reg():
    bed_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/melanogaster_ont/my_pipe/step3_SV_consensus/candidate_op/candidate_op.bed"
    asm_len = 0
    big_len = 0
    cnt = 0
    big_cnt = 0
    with open(bed_in, "r") as f:
        for line in f:
            fields = line.strip().split()
            ctg, start, end, op = fields[0], fields[1], fields[2], fields[3]
            start, end = int(start), int(end)
            if op.endswith("asm"):
                cnt += 1
                if end - start > 100000:
                    big_len += end -start
                    big_cnt += 1
                    print("{}:{}-{}\t{}bp".format(ctg, start, end, end-start))
                asm_len += end - start
        print(asm_len, big_len)
        print("asm num:{}\n big num:{}".format(cnt, big_cnt))

# check_asm_reg()

def check_chm13():
    bed_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi/my_pipe/step3_SV_consensus/consensus.bed"
    pre = "NC_0000"
    for i in range(1,10):
        # print("\"" + pre + "0" + str(i) + "\"", end=",")
        print(pre + "0" + str(i), end=",")
    for i in range(10, 25):
        # print("\"" + pre + str(i) + "\"", end=",")
        print(pre + str(i), end=",")
    chr_ls = ["NC_000001", "NC_000002", "NC_000003", "NC_000004", "NC_000005", "NC_000006", "NC_000007", "NC_000008", "NC_000009", "NC_000010", "NC_000011", "NC_000012", "NC_000013", "NC_000014", "NC_000015", "NC_000016", "NC_000017", "NC_000018", "NC_000019", "NC_000020", "NC_000021", "NC_000022", "NC_000023", "NC_000024"]
    # ctg_to_reg = defaultdict(list)
    # ctg_to_len = defaultdict(int)
    # with open(bed_in, "r") as f:
    #     N_fill_cnt = 0
    #     for line in f:
    #         fields = line.strip().split()
    #         ctg, start, end, ops = fields[0], fields[1], fields[2], fields[3]
    #         start, end  = int(start), int(end)
    #         ctg = ctg.split(".")[0]
    #         op_ls = ops.split(",")
    #         # print(ctg, op_ls)
    #         if "N_fill" in op_ls:
    #             N_fill_cnt += 1
    #             if ctg in chr_ls:
    #                 ctg_to_reg[ctg].append([start, end])
    #     for key, val in ctg_to_reg.items():
    #         for reg in val:
    #             ctg_to_len[key] += reg[1] - reg[0]
    #     for key, val in ctg_to_len.items():
    #         print(key, val)
# check_chm13()

def ref_process():
    fa_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/chm13/REF/ref1/GCF_000001405.40_GRCh38.p14_genomic.fna"
    fa_out = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/chm13/REF/ref1/new_GRCH38.fna"
    import fasta_parser
    ctg_dic = fasta_parser.read_sequence_dict(fa_in)
    new_ctg_dic = {}
    for ctg, seq in ctg_dic.items():
        if ctg.startswith("NC_000"):
            new_ctg_dic[ctg] = seq
    fasta_parser.write_fasta_dict(new_ctg_dic, fa_out)
# ref_process()

def ref_stats():
    fa_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/chm13/REF/ref1/new_GRCH38.fna"
    import fasta_parser
    ctg_dic =  fasta_parser.read_sequence_dict(fa_in)
    chr_id = "NC_000019.10"
    for ctg, seq in ctg_dic.items():
        # if ctg != chr_id:continue   # 查看特定ctg
        ctg_len = len(seq)
        AGCT = 0
        N_num = 0
        for ch in seq:
            if ch in ("A", "G", "C", "T", "a", "g", "c", "t"):
                AGCT += 1
            elif ch == "N":
                N_num += 1
            else:
                # print(ch)
                pass
        print(ctg, "AGCT_num:{}\tN_num:{}, 占比：{}".format(AGCT, N_num, N_num/ctg_len))
ref_stats()


def get_ctgs(fai_in, prefix):
    with open(fai_in, "r") as f:
        ctg_ls = []
        for line in f:
            ctg = line.strip().split()[0]
            if ctg.startswith(prefix):
                ctg_ls.append(ctg)
        print(ctg_ls, "\n")
        for ctg in ctg_ls:
            print(ctg, end=",")
        print("\n")
        for ctg in ctg_ls:
            print(ctg, end="\t")
        print()
fai_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/chm13/REF/ref1/GCF_000001405.40_GRCh38.p14_genomic.fna.fai"
prefix = "NC_00"     # 用于筛选的prefix
# get_ctgs(fai_in, prefix)