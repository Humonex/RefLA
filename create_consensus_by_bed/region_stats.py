import numpy as np
bed_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi/my_pipe/step3_SV_consensus/consensus.bed"
ls = []
len_ls = []
gap_ls = []
pre_ctg = None
pre_start = -1
pre_end = -1
with open(bed_in, "r") as f:
    for line in f:
        fields = line.strip().split("\t")
        ctg, start, end, op_ls  = fields[0], int(fields[1]), int(fields[2]), fields[3].split(",")
        if "N_fill" in op_ls:
            if ctg == pre_ctg:
                # print("gap: ", start - pre_end)
                gap_ls.append(start - pre_end)
                if start - pre_end > 500000:
                    print("big gap:", ctg, pre_start, pre_end, "-", ctg, start, end)    # 
            ls.append([ctg, start, end])
            len_ls.append(end - start)
            ## 
            if end - start < 60950:
                print(ctg, start, end)
        if "N_fill" in op_ls or "asm_patch" in op_ls:
            pre_ctg = ctg
            pre_start = start
            pre_end = end
len_ls.sort()
gap_ls.sort()
print(ls)
print(len_ls)
print("gap_ls:", gap_ls)
print(np.quantile(len_ls, [0.25, 0.5, 0.75]))
print(np.quantile(gap_ls, [0.25, 0.5, 0.75, 0.95]))