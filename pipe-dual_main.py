import argparse
import os
import sys
import pysam
# import signal
# import multiprocessing as mp
from multiprocessing import Pool
import random
import time
import re
import subprocess
from collections import namedtuple,defaultdict
from Bio import Seq
import shutil
import logging
import gzip
import yaml
from find_candidate_regions.find_candidate_regions import run_find_candidate_parallel, run_find_candidate2_parallel
from create_consensus_by_bed.get_fasta_consensus2 import run_SVconsensus_parallel
from create_consensus_by_bed import fasta_parser, scaffold_by_rec, Utils
from diploid_asm.whatshap_pipe import run_phasing

logger = logging.getLogger()
def _enable_logging(log_file, debug, overwrite):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    if not debug:
        console_log.setLevel(logging.INFO)

    if overwrite:
        open(log_file, "w").close()
    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)

def file_check(path):
    if not os.path.isfile(path):
        logger.error("Missing output: %s", path)
        raise Exception("Missing output")
        
def convert_size(size:str):
    if size.endswith("G") or size.endswith("g"):
        return 1000000*float(size[:-1])
    elif size.endswith("M") or size.endswith("m"):
        return 1000*float(size[:-1])
    elif size.endswith("B") or size.endswith("b"):
        return float(size[:-1])
    else:
        raise ValueError("Error format of genome size!!!")


def run_minimap2(target_fa, query_fa, data_type, threads, minimap_out, params_ls):     
    ''' 运行minimap2 提供params '''
    t1 = time.time()
    params = " ".join([str(param) for param in params_ls])
    minimap_cmd = ["minimap2", "-ax", "map-"+data_type, params, target_fa, query_fa, "-t", str(threads), ">", minimap_out]
    logger.info("Running: %s", " ".join(minimap_cmd))
    subprocess.check_call(" ".join(minimap_cmd), shell=True)
    logger.info("Run minimap2 finished, cost {}s".format(time.time() - t1))

def run_samtools(sam_in, samtools_out, threads):
    t1 = time.time()
    samtools_cmd = ["samtools", "sort", "-O", "BAM", "-@", str(threads), sam_in, "-o", samtools_out, "&&", "samtools", "index", "-@", str(threads), samtools_out]
    logger.info("Running: %s", " ".join(samtools_cmd))
    subprocess.check_call(" ".join(samtools_cmd), shell=True)
    logger.info("Run samtools finished, cost {}s".format(time.time() - t1))

def print_step_time(t0, t1, step_name):
    print("####STEP {} cost: {}".format(step_name, t1 - t0))
def print_error(value):
    print("error: ", value)


def main():
    parser = argparse.ArgumentParser(description="dual reference guided asssembly")
    #####
    ### ,MIN_CLIP_NUM, MIN_DEP, MIN_DEP_REG
    # parser = argparse.ArgumentParser(description="get_candidate_regions")
    group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument("-t", "--threads", dest="threads", type=int, default=1)
    group.add_argument("--ctg-ls", dest="ctg_ls", help="ctg list to perform assembly, recommend to provide contig ls, only accept format like:chr1,chr2")   # contig list，影响结果
    group.add_argument("-a", dest="all_chrs", default=False, action='store_true', help="we will keep and process all chrs")
    ### 路径
    parser.add_argument("--work-dir", dest="work_dir", required=True, help="work directory to output results")
    # parser.add_argument("--bam", dest="bam_in", required=True, help="bam file")
    parser.add_argument("--ref", dest="ref", required=True)
    parser.add_argument("--fastq", dest="fastq_in", required=True, help="fastq reads file")
    ## 
    parser.add_argument("--data-type", dest="data_type", required=True, choices=["ont", "hifi"], help="fastq file type")
    parser.add_argument("-g", "--genome-size", dest="genome_size", default="100m", help="genome size end with: b,m,g")
    ### 性能参数
    parser.add_argument("--min_clip_num", dest="min_clip_num", default=5, help="min_clip_num of window to be selected") # important
    parser.add_argument("--min_clip_len", dest="min_clip_len", default=500, help="min_clip_len of clip to be selected as clip")
    parser.add_argument("--dp-win-size", dest="dp_win_size", default=100, help="dp win_size")
    parser.add_argument("--min-dep", dest="min_dep", default=10, help="dep threshold to be select as low dep region")
    parser.add_argument("--max-dep", dest="max_dep", default=48, help="max dep threshold")  # ->后面改为随由平均cov来计算 (1.7 * avg cov)
    parser.add_argument("--min-dep-reg", dest="min_dep_reg", default=100, help="minimum length of low dep region")
    parser.add_argument("--min-contig", dest="min_contig", help="skip SV consensus process contig shorter than this, keep with raw", default=100, type=int)   # 200000 调 1,000,000     5,000,000
    parser.add_argument("--fill-size", dest="Nfill_size", type=int, default=100, help="N_fill_size")
    # 区间聚类和延展参数
    # parser.add_argument("--cluster-params", dest="cluster_params")
    ### 选择性的性能参数，
    parser.add_argument("--ex-unmapped-denovo", dest="ex_unmapped_denovo", default=False, action="store_true", help="exclude unmapped denovo asm, may be unique sequence of sample")    # 默认保留，获得更加完整的组装结果
    '''Test'''
    # parser.add_argument("--cut-gap", dest="cut_gap", default=False, action="store_true", help="cut between big gap instead of N_fill")      
    parser.add_argument("--polish-filter", dest="polish_filter", default=False, action="store_true")
    ## 
    parser.add_argument("--overwrite", dest="overwrite", default=False, action="store_true", help="overwrites all existing results")
    parser.add_argument("-R", dest="re_run", choices=["step1", "step2", "step3", "step4", "step5"], default=False, help="re_run step")
    ## 一些参数的配置文件
    parser.add_argument("--config", default='/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/Configs/Config.yaml')
    ## mode，二倍体模式：普通的单倍体组装/二倍体组装
    parser.add_argument("--dual", dest="dual", default=False, action="store_true", help="open dual mode")
    parser.add_argument("--parallel", dest="parallel", default=False, action="store_true", help="open parallel mode")
    parser.add_argument("--diploid_config", default='/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/diploid_asm/Configs/diploid.yaml')
    parser.add_argument("--phasing_script", default="/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/diploid_asm/whatshap_pipe/run_phasing.sh")
    # 
    args = parser.parse_args()
    

    logger.info("Start RefBasm assembly pipe")
    ### 绝对路径配置
    args.work_dir = os.path.abspath(args.work_dir)
    args.ref = os.path.abspath(args.ref)
    args.fatsq = os.path.abspath(args.fastq_in)
    ### 生成文件路径配置
    work_dir = args.work_dir
    ref_dir = work_dir + "/" + "corrected_ref"
    '''
    step1_dir = work_dir + "/" + "step1_mapping"
    step2_dir = work_dir + "/" + "step2_candidate_regions"
    step3_dir = work_dir + "/" + "step3_SV_consensus"
    step4_dir = work_dir + "/" + "step4_polish"
    step5_dir = work_dir + "/" + "step5_scaffolding"
    '''
    log_dir = work_dir + "/" + "logs"
    config_dir = work_dir + "/" + "config"
    # dp_file_dir = work_dir + "/" + "depths"
    '''for Dir in [work_dir, ref_dir, dp_file_dir, step1_dir, step2_dir, step3_dir, step4_dir, step5_dir, log_dir, config_dir]:
        Utils.make_dir(Dir)'''
    Utils.make_dir(work_dir)
    Utils.make_dir(ref_dir)
    Utils.make_dir(log_dir)
    Utils.make_dir(config_dir)
    ## log和文件检查
    assembly_log = os.path.join(log_dir, "assembly.log")
    _enable_logging(assembly_log, debug=False, overwrite=True)
    if args.data_type not in ["ont", "hifi"]:
        print("Reads type could be either 'ont' or 'hifi'", file=sys.stderr)
        return 1
    file_check(args.fastq_in)
    file_check(args.ref)
    ## 一些参数的获取，和对参考序列的处理
    with open(args.config, "r") as f: # config参数获取
        config = yaml.safe_load(f.read())     # 获取部分参数
        logger.info("Config参数设置: {}".format(config))
    config_save = os.path.join(config_dir, "WorkConfig.yaml")
    shutil.copy(args.config, config_save)
    # process reference
    raw_ref_dic = fasta_parser.read_sequence_dict(args.ref)    # 参考序列字典
    corrected_ref_dic = {}
    ctg_ls = [] # 储存进行SV consensus 的contig
    if args.all_chrs:
        for chr_id,ref_seq in raw_ref_dic.items():
            corrected_ref_dic[chr_id] = ref_seq
            if len(ref_seq) > args.min_contig:
                ctg_ls.append(chr_id)
    else:   # 只保留指定的contig
        ls = args.ctg_ls.strip(",").split(",")
        for chr_id in ls:
            ref_seq = raw_ref_dic.get(chr_id)
            if not ref_seq: continue
            corrected_ref_dic[chr_id] = ref_seq
            if len(ref_seq) > args.min_contig:
                ctg_ls.append(chr_id)
    logger.info("corrected ref include:{}".format(corrected_ref_dic.keys()))
    logger.info("You will SV consensus on:{}".format(",".join(ctg_ls)))
    corrected_ref = ref_dir + "/" + "reference.fasta"
    fasta_parser.write_fasta_dict(corrected_ref_dic, corrected_ref)     # 写入新的参考基因组中
    faidx_cmd = ["samtools faidx", corrected_ref]
    subprocess.check_call(" ".join(faidx_cmd), shell=True)

    overwrite = args.overwrite
    re_run = args.re_run
    genome_size = convert_size(args.genome_size)
    if re_run: print("Will re run from {}".format(re_run))
    ### params 
    threads = args.threads
    print("Threads: {}".format(threads))
    data_type = args.data_type
    Nfill_size = args.Nfill_size
    ex_unmapped_denovo = args.ex_unmapped_denovo
######################################################### pipeline #########################################################
    
    '''## step1 mapping 
    t1 = time.time()
    minimap_out1 = os.path.join(step1_dir, "aln.sam")
    samtools_out1 = os.path.join(step1_dir, "aln.sorted.bam")
    if re_run == "step1": overwrite = True
    if os.path.isfile(samtools_out1) and not overwrite:
        logger.info("Skipped step1")
    else:
        logger.info("Run step1 mapping")
        minimap2_params_ls1 = config["step1"]["minimap2_params_ls"]
        # run_minimap(corrected_ref, args.fastq_in, args.data_type, threads, minimap_out1)
        try:
            run_minimap2(corrected_ref, args.fastq_in, args.data_type, threads, minimap_out1, minimap2_params_ls1)
        except:
            raise Exception("Run mapping minimap2 failed!!!")
        try:
            run_samtools(minimap_out1, samtools_out1, threads) # sam -> sorted.bam, and index
            os.remove(minimap_out1)
        except:
            raise Exception("Run mapping samtools failed!!!")
        overwrite = True'''
    
    ## for dual assembly, mapping、phasing
    # corrected_ref = args.ref
    if args.dual:
        # genome_size = genome_size / 2
        phasing_dir = work_dir + "/" + "Phasing"
        if not os.path.isdir(phasing_dir):os.makedirs(phasing_dir)
        '''
        dual mode file out
        '''
        hap1_dir = work_dir + "/" + "hap1"
        hap2_dir = work_dir + "/" + "hap2"
        Utils.make_dir(hap1_dir)
        Utils.make_dir(hap2_dir)
        hap1_samtools_out1 = os.path.join(phasing_dir, "whatshap_out", "Merged_out", "merged_hap_1.sorted.bam")
        hap2_samtools_out1 = os.path.join(phasing_dir, "whatshap_out", "Merged_out", "merged_hap_2.sorted.bam")
        hap1_fastq_in = os.path.join(phasing_dir, "whatshap_out", "Merged_out", "hap_1.fastq")
        hap2_fastq_in = os.path.join(phasing_dir, "whatshap_out", "Merged_out", "hap_2.fastq")
        
        ## 
        if re_run == "step1": overwrite = True
        if os.path.isfile(hap1_samtools_out1) and os.path.isfile(hap1_fastq_in) and not overwrite:
            logger.info("Skipped Phasing")
        else:
            overwrite = True
            logger.info("Run Phasing")
            t_0 = time.time()
            try:
                with open(args.diploid_config, "r") as f: # config参数获取
                    diploid_config = yaml.safe_load(f.read())
                res = run_phasing.phasing(diploid_config, args.phasing_script, phasing_dir, threads, args.ref, args.fastq_in, args.data_type)
                if res != 0:
                    raise Exception("Run Phasing Failed!!!")
            except:
                raise Exception("Run Phasing Failed!!!")
            t_1 = time.time()
            print_step_time(t_0, t_1, "Run Phasing")
    else:
        raise Exception("Pipe only support diploid mode!!!")

    if args.parallel and threads > 2:
        # threads = threads // 2
        print("Parallel run steps")
        pool = Pool(processes=2)
        for hap in [1, 2]:
            if hap == 1:
                hap_dir = hap1_dir
                samtools_out1 = hap1_samtools_out1
                fastq_in = hap1_fastq_in
            elif hap == 2:
                hap_dir = hap2_dir
                samtools_out1 = hap2_samtools_out1
                fastq_in = hap2_fastq_in
            pool.apply_async(new_func, args=(config, ctg_ls, corrected_ref, re_run, genome_size, threads, hap, hap_dir, samtools_out1, fastq_in, data_type, Nfill_size, ex_unmapped_denovo), error_callback=print_error)
        pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
        pool.join() # 等待进程池中的所有进程执行完毕
    else:
        for hap in [1, 2]:
            if hap == 1:
                hap_dir = hap1_dir
                samtools_out1 = hap1_samtools_out1
                fastq_in = hap1_fastq_in
            elif hap == 2:
                hap_dir = hap2_dir
                samtools_out1 = hap2_samtools_out1
                fastq_in = hap2_fastq_in
            new_func(config, ctg_ls, corrected_ref, re_run, genome_size, threads, hap, hap_dir, samtools_out1, fastq_in, data_type, Nfill_size, ex_unmapped_denovo)

    '''## time stastic
    time_ls = [t1, t2, t3, t4, t5, t6, t7]
    print("Time stastic of the pipe: ")
    print("\t".join(["mapping raw reads,", "Find candidate regions,", "SVconsensus,", "mapping,", "polish,", "scaffolding"]))
    for i in range(len(time_ls)-1):
        print("t{}-t{}: {}s".format(i + 1, i + 2, time_ls[i + 1] - time_ls[i]))'''


def new_func(config, ctg_ls, corrected_ref, re_run, genome_size, threads, hap, hap_dir, samtools_out1, fastq_in, data_type, Nfill_size, ex_unmapped_denovo):
    logger.info("---------------------------Run hap{}---------------------------".format(str(hap)))

    ## step2 Find candidate regions
    step2_dir = hap_dir + "/" + "step2_candidate_regions"
    Utils.make_dir(step2_dir)
    candidate_bed = os.path.join(step2_dir, "candidate.bed")
    if re_run == "step2": overwrite = True
    if os.path.isfile(candidate_bed) and not overwrite:
        logger.info("Skipped step2")
    else:
        logger.info("Run step2 Find candidate regions")
        # run_find_candidate_parallel(samtools_out1, ctg_ls, step2_dir, dp_file_dir, threads, args.min_clip_num, args.min_clip_len, args.dp_win_size, args.min_dep, args.max_dep, args.min_dep_reg)
        t_0 = time.time()
        try:
            run_find_candidate2_parallel(corrected_ref, samtools_out1, ctg_ls, step2_dir, threads, config["step2"])
        except:
            raise Exception("Run step2 Failed!!!")
        t_1 = time.time()
        print_step_time(t_0, t_1, "run_find_candidate2_parallel")
        overwrite = True

    ## step3 SVconsensus
    # t3 = time.time()
    # print("STEP2 cost: ", t3 - t2)
    step3_dir = hap_dir + "/" + "step3_SV_consensus"
    Utils.make_dir(step3_dir)
    SV_consensus_fasta = os.path.join(step3_dir, "consensus.fasta")
    SV_consensus_bed = os.path.join(step3_dir, "consensus.bed")
    merge_fasta = os.path.join(step3_dir, "merge.fasta")
    if re_run == "step3": overwrite = True
    if os.path.isfile(SV_consensus_fasta) and os.path.isfile(SV_consensus_bed) and os.path.isfile(merge_fasta) and not overwrite:
        logger.info("Skipped step3")
    else:
        bam_in = samtools_out1
        logger.info("Run step3 SV consensus")
        t_0 = time.time()
        try:
            data_type, Nfill_size, ex_unmapped_denovo
            run_SVconsensus_parallel(step3_dir, bam_in, ctg_ls, candidate_bed, threads, fastq_in, data_type, corrected_ref, Nfill_size, genome_size, ex_unmapped_denovo, config["step3"])
        except:
            raise Exception("Run step3 Failed!!!")
        t_1 = time.time()
        print_step_time(t_0, t_1, "run_SVconsensus_parallel")
        overwrite = True

    ## step4 polish
    step4_dir = hap_dir + "/" + "step4_polish"
    Utils.make_dir(step4_dir)
    polish_out = os.path.join(step4_dir, "racon.fasta")
    if re_run == "step4": overwrite = True
    if os.path.isfile(polish_out) and not overwrite:
        logger.info("Skipped step4")
    else:
        logger.info("---------------------------Run step4 polish---------------------------")
        '''minimap_cmd = ["minimap2", "-ax", "map-"+data_type, params, target_fa, query_fa, "-t", str(threads), ">", minimap_out]'''
        '''racon {input.fastq} {output.sam} {input.asm} -t {threads} > {output.fasta}'''
        minimap_out2 = os.path.join(step4_dir, "aln.sam")
        minimap2_params_ls2 = config["step4"]["minimap2_params_ls"]
        racon_params = " ".join(config["step4"]["racon_params_ls"])
                
        ## 
        logger.info("Run mapping to consensus fasta")
        t_0 = time.time()
        try:
            run_minimap2(merge_fasta, fastq_in, data_type, threads, minimap_out2, minimap2_params_ls2)
        except:
            raise Exception("Run mapping to consensus fasta failed!!!")
        t_1 = time.time()
        print_step_time(t_0, t_1, "Run mapping to consensus fasta")

        logger.info("Run racon polish")
        ''' racon: -u, --include-unpolished   output unpolished target sequences'''
        polish_cmd = ["racon", racon_params, fastq_in, minimap_out2, merge_fasta, "-t", str(threads), ">", polish_out]
        logger.info("Running: %s", " ".join(polish_cmd))
        t_0 = time.time()
        try:
            subprocess.check_call(" ".join(polish_cmd), shell=True)
        except:
            raise Exception("Run racon polish failed!!!")
        t_1 = time.time()
        print_step_time(t_0, t_1, "Run racon polish")
        overwrite = True

    ## step5 scaffolding
    step5_dir = hap_dir + "/" + "step5_scaffolding"
    Utils.make_dir(step5_dir)
    connect_file = os.path.join(step3_dir, "connect.txt")
    scaffold_out = os.path.join(step5_dir, "final.fasta")
    if re_run == "step5": overwrite = True
    if os.path.isfile(scaffold_out) and not overwrite:
        logger.info("Skipped all!!!")
    else:
        logger.info("Run step5 scaffolding")
        t_0 = time.time()
        scaffold_by_rec.scaffold(polish_out, connect_file, scaffold_out, config["step5"])
        t_1 = time.time()
        print_step_time(t_0, t_1, "Run step5 scaffolding")
        logger.info("Run scaffolding Done!!!")



if __name__ == "__main__":
    main()
    pass
