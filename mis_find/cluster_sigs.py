

import pysam


def cluster_all(bam, reg, out_dir):
    print("cluster for:", reg)
    bam_reader = pysam.AlignmentFile(bam, "rb", index_filename=bam + ".bai")
    sigs_file = "%s/signatures/_%s_%d_%d.bed"%(out_dir, reg[0], reg[1], reg[2])
    cluster_file = "%s/signatures/_%s_%d_%d.bed"%(out_dir, reg[0], reg[1], reg[2])
    final_ls = []
    ins_ls = []
    del_ls = []
    inv_ls = []
    tra_ls = []
    other_ls = []
    with open(sigs_file, "r") as f:
        for line in f:
            if line.startswith("#"): continue
            if line.startswith("INS"): ins_ls.append(line.strip().split("\t"))
            if line.startswith("DEL"): del_ls.append(line.strip().split("\t"))
            if line.startswith("INV"): continue
    ## 
    lower_dp = 5
    if reg[2] - reg[1] > 20000:
        final_ls.append([reg[0], reg[1], reg[2], "Large", ])
    pass



def clustar_inv():
    pass

if __name__ == "__main__":
    bed = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/my_pipe/step2_candidate_regions/candidate/merge.clu1.bed"
    reg_ls = []
    with open(bed, "r") as f:
        for line in f:
            if line.startswith("#"): continue
            reg = line.strip().split("\t")
            reg[1] = int(reg[1])
            reg[2] = int(reg[2])
            reg_ls.append(reg)
    bam = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/aln2alt.sort.bam"
    ref = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/GCF_000146045.fna"
    out_dir = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/my_pipe/step2_candidate_regions"
    
    for reg in reg_ls:
        cluster_all(bam, reg, out_dir)
    pass