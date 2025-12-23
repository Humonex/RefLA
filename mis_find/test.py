

# bed = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/my_pipe/step2_candidate_regions/candidate/false.bed2"
bed = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/my_pipe/step2_candidate_regions/candidate/merge.clu1.bed2"
with open(bed, "r") as f:
    ls = []
    for line in f:
        ctg, start, end = line.strip().split("\t")[:3]
        start, end = int(start), int(end)
        ls.append(end - start)
    # print(sorted(ls))
bed1 = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/my_pipe/step2_candidate_regions/candidate/false.bed2"
bed2 = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/my_pipe/step2_candidate_regions/candidate/merge.clu1.bed2"
bed3 = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/my_pipe/step2_candidate_regions/candidate/false.bed2_info" 
with open(bed1, "r") as f1, open(bed2, "r") as f2, open(bed3, "w") as f3:
    ls = []
    for line in f1:
        ctg, start, end = line.strip().split("\t")[:3]
        # start, end = int(start), int(end)
        ls.append(ctg + start + end)
    print(len(ls))
    for line in f2:
        ctg, start, end = line.strip().split("\t")[:3]
        if ctg + start + end in ls:
            f3.write(line)