

from collections import defaultdict
import math
from multiprocessing import Pool
import random
import subprocess
import time
import pysam
import logging
logger = logging.getLogger()

def reg_to_id(reg):
    return reg[0] + ":" + str(reg[1]) + "-" + str(reg[2])
def write_read_ids(file, read_ids):
    ls = [read_id + "\n" for read_id in read_ids]
    with open(file, "w") as f:
            f.writelines(ls)

def get_asm_read_ids(bam_in, threads, reg, config):
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    ids_set = set()
    # for reg in regions:
    if config["remove_secondary"]:
        for read in bam_reader.fetch(reg[0], reg[1], reg[2]):   # exclude secondary
            if read.is_secondary: continue
            ids_set.add(read.query_name)
    else:
        for read in bam_reader.fetch(reg[0], reg[1], reg[2]):   # 加入了所有的序列
            ids_set.add(read.query_name)
    return "asm", reg, ids_set

def get_clip_read_ids(bam_in, threads, reg, mm_num, config):
    '''
    问题：容易收集过多不应该收集的读数，导致大量冗余的错误组装，因此引入条件，确保该区域为目标区域
    1、存在子区域1000bp: depth>2*whole_dp
    2、
    '''
    # --------------get high clip reads ids-------------- #
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    # print("find {}".format(reg))
    min_clip_len_portion = config["high_clip"]["min_clip_len_portion"]   # min_clip_len_portion
    # min_clip_len_portion = 0.2
    ids_set = set()
    # max_MAPQ = 60
    # for reg in regions:
    for read in bam_reader.fetch(reg[0], reg[1], reg[2]):
        # if read.is_secondary or read.is_supplementary or read.is_unmapped: continue
        if read.is_secondary or read.is_unmapped: continue
        # if read.mapping_quality > max_MAPQ: continue    # 该条件是否添加
        cigar = read.cigartuples
        clip_length = 0
        left = cigar[0]
        if left[0] == 4 or left[0] == 5:
            clip_length += left[1]
        right = cigar[-1]
        if right[0] == 4 or right[0] == 5:
            clip_length += right[1]
        if read.query_length > 5000 and clip_length / read.query_length > min_clip_len_portion:
            ids_set.add(read.query_name) # add read
    if len(ids_set) < mm_num:
        ids_set = set()
    return "clip", reg, ids_set

def get_region_read_ids(work_dir, bam_in, threads, candidate_op_ls, whole_dp, config):
    '''
    return: reg_read_ids_dic, 
    "all_read"
    '''
    t0 = time.time()
    asm_regions = []
    read_regions = []
    whole_ctg_asm_regions = []
    for rec in candidate_op_ls:
        if rec.operation.endswith("asm"):   # reg_asm、whole_ctg_asm
            asm_regions.append([rec.chr_id, rec.start, rec.end])
        elif rec.operation.endswith("reads"):   # span_reads
            read_regions.append([rec.chr_id, rec.start, rec.end])
        else:
            raise ValueError
    ## 
    clip_reg_ls = []
    reg_read_ids_dic = defaultdict(set)     # asm_reg_id -> reg_read_ids
    clip_read_ids_dic = defaultdict(set)     # asm_reg_id -> reg_read_ids
    all_ids = set()
    asm_ids = set()
    clip_ids = set()
    results = []
    pool = Pool(processes=threads)
    # # results = [pool.apply_async(SV_consensus_on_ref, args=(ctg, candidate_op_dic[ctg], ref_dic[ctg], asm_fa_dic, asm_to_ref_bam, Nfill_size)) for ctg in all_chrs]
    # results = [pool.apply_async(SV_consensus_on_ref, args=(ctg, semi_candidate_op_dic[ctg], ref_dic, asm_fa_dic, asm_to_ref_bam, Nfill_size)) for ctg in all_chrs]
    if config["high_clip"]["apply_on_high_clip"]:
        mm_num = max(math.ceil(config["high_clip"]["min_high_clip_num"]), config["high_clip"]["min_portion"] * whole_dp)
        print("mm_num:{}".format(mm_num))
        for reg in read_regions:
            results.append(pool.apply_async(get_clip_read_ids, args=(bam_in, threads, reg, mm_num, config)))
    for reg in asm_regions:
        results.append(pool.apply_async(get_asm_read_ids, args=(bam_in, threads, reg, config)))
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    for res in results:
        flag, reg, ids_set = res.get()
        reg_id = reg_to_id(reg)
        if flag == "asm":
            reg_id = reg_to_id(reg)
            reg_read_ids_dic[reg_id] = ids_set
            asm_ids.update(ids_set)
            # print("Asm reg:{}".format(reg))
        else:
            clip_read_ids_dic[reg_id] = ids_set
            # clip_reg_ls.append(reg)
            clip_ids.update(ids_set)
    all_ids.update(asm_ids)
    all_ids.update(clip_ids)
    reg_read_ids_dic["all_read"] = all_ids
    reg_read_ids_dic["clip_read"] = clip_ids
    reg_read_ids_dic["asm_read"] = asm_ids
    
    # print("Read Stats:")
    # stats_file = work_dir + "/" + "denovo_read.stats"
    # with open(stats_file, "w") as f:
    #     f.write("Asm read num:{}\n".format(len(asm_ids)))
    #     f.write("Clip read num:{}\n".format(len(clip_ids)))
    #     f.write("All read num:{}\n".format(len(all_ids)))
    #     f.write("Clip read add num:{}\n".format(len(all_ids)-len(asm_ids)))
    print("Get read done!!!")
    print("Get denovo read cost:{}".format(time.time() - t0))
    # print(reg_read_ids_dic.keys())
    return reg_read_ids_dic, clip_read_ids_dic

def select_reads_from_names(fastq_in, fastq_out, read_ids, threads):
    ## 使用seqkit实现 
    '''seqkit grep -f id.txt -j threads seqs.fq.gz -o result.fq.gz'''
    # seqkit_cmd = ["seqkit", "grep", "-f", read_ids, "-j", str(threads), fastq_in, "-o", fastq_out]
    # logger.info("Running: %s", " ".join(seqkit_cmd))
    cmd = ["seqtk", "subseq", fastq_in, read_ids, ">", fastq_out]   ## bbmap filterbyname.sh 命令好像更快 
    logger.info("Running: %s", " ".join(cmd))
    subprocess.check_call(" ".join(cmd), shell=True)


def get_reg_ids(bam_in, reg_ls):
    read_ids = set()
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    for reg in reg_ls:
        for read in bam_reader.fetch(reg[0], reg[1], reg[2]):
            read_ids.add(read.query_name)
    return read_ids

def get_reg_ids_parallel(bam_in, reg_ls, threads):
    random.shuffle(reg_ls)
    batch_size = len(reg_ls) // (threads + 1)   # 划分batch
    batches = []
    for i in range(threads):
        batch_reg_ls = reg_ls[i*batch_size:(i+1)*batch_size]
        batches.append(batch_reg_ls)
    results = []
    pool = Pool(processes=threads)
    results = [pool.apply_async(get_reg_ids, args=(bam_in, batch_reg_ls)) for batch_reg_ls in batches]
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    ids = set()
    for res in results:
        ids.update(res.get())

if __name__ == "__main__":
    import yaml
    high_clip_file = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/sativa_hifi/test/my_pipe2/step3_SV_consensus/high_clip.bed"
    bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/sativa_hifi/test/my_pipe2/step1_mapping/aln.sorted.bam"
    threads = 40 
    mm_num = 6
    config_f = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/Configs/Config.yaml"
    bed = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/sativa_hifi/test/my_pipe2/step3_SV_consensus/candidate_op/candidate_op.bed"
    reg_ls = []
    with open(config_f, "r") as f: # config参数获取
        config = yaml.safe_load(f.read())     # 获取部分参数
    config = config["step3"]
    # print(config)
    with open(bed, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            ctg, start, end, op = fields[:4]
            if op.endswith("reads"):
                reg_ls.append([ctg, int(start), int(end)])
                # print(ctg, start, end, op)
    print("Find")
    # with open(high_clip_file, "w") as f:
    #     for reg in reg_ls:
    #         flag, reg, ids_set = get_clip_read_ids(bam_in, threads, reg, mm_num, config)
    #         if len(ids_set) > 0:
    #             print(reg, len(ids_set), ",".join(ids_set))
    #             f.write("{}\t{}\t{}\t{}\t{}\n".format(reg[0], reg[1], reg[2], len(ids_set), ",".join(ids_set)))
    results = []
    pool = Pool(processes=threads)
    for reg in reg_ls:
        results.append(pool.apply_async(get_clip_read_ids, args=(bam_in, threads, reg, mm_num, config)))
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    with open(high_clip_file, "w") as f:
        for res in results:
            flag, reg, ids_set = res.get()
            if len(ids_set) > 0:
                print(reg, len(ids_set), ",".join(ids_set))
                f.write("{}\t{}\t{}\t{}\t{}\n".format(reg[0], reg[1], reg[2], len(ids_set), ",".join(ids_set)))
    pass
