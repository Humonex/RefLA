
# bed_in = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/create_consensus_by_bed/operation.bed"
bed_in = "/public/home/hpc214712170/Test/tests/chm13_hifi_ctgs/NC_18/my_pipe2/step3_SV_consensus/consensus.bed"
# bed_in = "/public/home/hpc214712170/Test/tests/chm13_hifi/my_pipe/step3_SV_consensus/consensus.bed"
with open(bed_in, "r") as fin:
    for line in fin:
        # ctg, start, end, INFO, operation = line.strip().split()[:5]
        fields = line.strip().split("\t")
        ctg, start, end, op_ls, info_ls, patch_ls = fields[:6]
        start, end = int(start), int(end)
        if "N_fill" in op_ls:
            patch_len = len(info_ls)
            # print("patch length:", patch_len)
            print(ctg, start, end, op_ls,":", end - start, "bp", "\tpatch_len_diff: ", end - start - patch_len, "bp")

