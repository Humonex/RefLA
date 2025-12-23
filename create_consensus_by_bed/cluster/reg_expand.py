bed_in = ""

all_ls = []
ctg_asm_ls = []
with open(bed_in, "r")as f:
    for line in f:
        fields = line.strip().split("\t")
        ctg, start, end, op = fields[0], int(fields[1]), int(fields[2]), fields[3]
        if op.endswith("asm"):
            if ctg == pre_ctg:
                # if end - start > big_reg_threshold:
                ctg_asm_ls.append([ctg, start, end])
            else:
                all_ls.append(ctg_asm_ls)
                ctg_asm_ls = [[ctg, start, end]]
                pre_ctg = ctg
all_ls.append(ctg_asm_ls)

for ctg_asm_ls in all_ls:
    for reg in ctg_asm_ls:
        if start >= 0: