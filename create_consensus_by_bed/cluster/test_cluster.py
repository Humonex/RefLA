out = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/create_consensus_by_bed/cluster/cluster.bed2"
bed_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi/my_pipe/step3_SV_consensus/candidate_op/candidate_op.bed"
all_ls = []
ctg_asm_ls = []
pre_ctg = None
big_reg_threshold = 50000     # 100000 | 
fo = open(out, "w")
with open(bed_in, "r") as f:
    for line in f:
        fields = line.strip().split("\t")
        ctg, start, end, op = fields[0], int(fields[1]), int(fields[2]), fields[3]
        if op.endswith("asm"):
            if ctg == pre_ctg:
                if end - start > big_reg_threshold:
                    ctg_asm_ls.append([ctg, start, end])
            else:
                all_ls.append(ctg_asm_ls)
                ctg_asm_ls = [[ctg, start, end]]
                pre_ctg = ctg
all_ls.append(ctg_asm_ls)

cluster_dis = 1000000   # 1000000
ls = []
cnt = 0
for ctg_asm_ls in all_ls:
    ctg, start, end = None, -1, -1
    cnt += len(ctg_asm_ls)
    print(ctg_asm_ls)
    for reg in ctg_asm_ls:
        if start >= 0:
            if reg[1] - end < cluster_dis:
                end  = reg[2]
            else:
                ls.append([ctg, start, end])
                fo.write("{}\t{}\t{}\n".format(ctg, start, end))
                start, end = reg[1], reg[2]
        else:
            ctg, start, end = reg[0], reg[1], reg[2]

    if start >= 0:
        ls.append([ctg, start, end])
        fo.write("{}\t{}\t{}\n".format(ctg, start, end))

fo.close()
print(ls)
print(cnt, len(ls))