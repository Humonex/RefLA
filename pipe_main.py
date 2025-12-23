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
from mis_find.find_mis_pipe import run_find_pipe
import info_stats
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
    minimap_cmd = ["/usr/bin/time -v minimap2", "-ax", "map-"+data_type, params, target_fa, query_fa, "-t", str(threads), ">", minimap_out]
    logger.info("Running: %s", " ".join(minimap_cmd))
    subprocess.check_call(" ".join(minimap_cmd), shell=True)
    logger.info("Run minimap2 finished, cost {}s".format(time.time() - t1))

def run_samtools(sam_in, samtools_out, threads):
    t1 = time.time()
    samtools_cmd = ["/usr/bin/time -v samtools", "sort", "-O", "BAM", "-@", str(threads), sam_in, "-o", samtools_out, "&&", "samtools", "index", "-@", str(threads), samtools_out]
    logger.info("Running: %s", " ".join(samtools_cmd))
    subprocess.check_call(" ".join(samtools_cmd), shell=True)
    logger.info("Run samtools finished, cost {}s".format(time.time() - t1))

def generate_bam(target_fa, query_fa, data_type, threads, bam_out, params_ls):
    ''' 运行minimap2 提供params，输出sorted bam '''
    t1 = time.time()
    params = " ".join([str(param) for param in params_ls])
    minimap_cmd = ["/usr/bin/time -v minimap2", "-ax", "map-"+data_type, params, target_fa, query_fa, "-t", str(threads), "|", "samtools", "sort", "-O", "BAM", "-@", str(threads), "-o", bam_out, "&&", "samtools", "index", "-@", str(threads), bam_out]
    logger.info("Running: %s", " ".join(minimap_cmd))
    subprocess.check_call(" ".join(minimap_cmd), shell=True)
    logger.info("Run minimap2 and samtools finished, cost {}s".format(time.time() - t1))


def main():
    parser = argparse.ArgumentParser(description="reference guided asssembly")
    #####
    ### ,MIN_CLIP_NUM, MIN_DEP, MIN_DEP_REG
    parser = argparse.ArgumentParser(description="get_candidate_regions")
    group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument("-t", "--threads", dest="threads", type=int, default=1)
    group.add_argument("--ctg-ls", dest="ctg_ls", help="ctg list to perform assembly, recommend to provide contig ls, only accept format like:chr1,chr2")   # contig list，影响结果
    group.add_argument("-a", dest="all_chrs", default=False, action='store_true', help="we will keep and process all chrs")
    ### 路径
    parser.add_argument("--work-dir", dest="work_dir", required=True, help="work directory to output results")
    # parser.add_argument("--bam", dest="bam_in", required=True, help="bam file")
    parser.add_argument("--ref", dest="ref", required=True)
    parser.add_argument("--fastq", dest="fastq", required=True, help="fastq reads file")
    ## 
    parser.add_argument("--data-type", dest="data_type", required=True, choices=["ont", "hifi", "pb"], help="fastq file type")
    parser.add_argument("-g", "--genome-size", dest="genome_size", default="100m", help="genome size end with: b,m,g")
    ### 性能参数
    parser.add_argument("--min_clip_num", dest="min_clip_num", default=5, help="min_clip_num of window to be selected") # important
    parser.add_argument("--min_clip_len", dest="min_clip_len", default=500, help="min_clip_len of clip to be selected as clip")
    parser.add_argument("--dp-win-size", dest="dp_win_size", default=100, help="dp win_size")
    parser.add_argument("--min-dep", dest="min_dep", default=10, help="dep threshold to be select as low dep region")
    parser.add_argument("--max-dep", dest="max_dep", default=48, help="max dep threshold")  # ->后面改为随由平均cov来计算 (1.7 * avg cov)
    parser.add_argument("--min-dep-reg", dest="min_dep_reg", default=100, help="minimum length of low dep region")
    parser.add_argument("--min-contig", dest="min_contig", help="skip SV consensus process contig shorter than this, keep with raw", default=20000, type=int)   # 200000 调 1,000,000     5,000,000
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
    ## 1、mis find mode
    parser.add_argument("--find_mis", dest="find_mis", default=False, action="store_true", help="Start find mode")
    parser.add_argument("--min_span", dest="min_span", type=int, default=5)
    # parser.add_argument("--boundary", dest="boundary", default=3000, type=int)
    parser.add_argument("--bam", dest="bam")
    ## 2、correct mode。纠错模式？？？
    parser.add_argument("--correct", dest="correct", default=False, action="store_true")
    parser.add_argument("--min_corr_len", dest="min_corr_len", type=int, default=500000)
    parser.add_argument("--skip_denovo", dest="skip_denovo", default=False, action="store_true")
    parser.add_argument("--cut", dest="cut", default=False, action="store_true")
    # parser.add_argument("")
    ## 一些参数的配置文件
    parser.add_argument("--config", default='/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/Configs/Config.yaml')
    ## 
    parser.add_argument("--out_denovo", dest="out_denovo", default=False, help="Provide path to the denovo fa, for testing")

    ##
    args = parser.parse_args()
    

    logger.info("Start my pipe to assembly")
    ### 绝对路径配置
    args.work_dir = os.path.abspath(args.work_dir)
    args.ref = os.path.abspath(args.ref)
    args.fatsq = os.path.abspath(args.fastq)
    ### 生成文件路径配置
    work_dir = args.work_dir
    ref_dir = work_dir + "/" + "corrected_ref"
    step1_dir = work_dir + "/" + "step1_mapping"
    step2_dir = work_dir + "/" + "step2_candidate_regions"
    step3_dir = work_dir + "/" + "step3_SV_consensus"
    step4_dir = work_dir + "/" + "step4_polish"
    step5_dir = work_dir + "/" + "step5_scaffolding"
    log_dir = work_dir + "/" + "logs"
    config_dir = work_dir + "/" + "config"
    dp_file_dir = work_dir + "/" + "depths"
    for Dir in [work_dir, ref_dir, dp_file_dir, step1_dir, step2_dir, step3_dir, step4_dir, step5_dir, log_dir, config_dir]:
        Utils.make_dir(Dir)
    ## log和文件检查
    assembly_log = os.path.join(log_dir, "assembly.log")
    _enable_logging(assembly_log, debug=False, overwrite=True)
    if args.data_type not in ["ont", "hifi", "pb"]:
        print("Reads type could be either 'ont' or 'hifi' or 'pb'", file=sys.stderr)
        exit(1)
        # return 1
    
    if args.find_mis or args.correct:
        file_check(args.bam)
    else:
        file_check(args.fastq)
    file_check(args.ref)
    file_check(args.config)
    ## 一些参数的获取
    with open(args.config, "r") as f: # config参数获取
        config = yaml.safe_load(f.read())     # 获取部分参数
        logger.info("Config参数设置: {}".format(config))
    config_save = os.path.join(config_dir, "WorkConfig.yaml")
    shutil.copy(args.config, config_save)
    
    ## 参考序列的处理
    corrected_ref = ref_dir + "/" + "reference.fasta"
    if args.find_mis:
        print("Find mis mode, Skip reference process")
        if not os.path.isfile(args.ref + ".fai"):
            pysam.faidx(args.ref)
        # shutil.copy(args.ref, corrected_ref)
        ctg_ls = []
        with open(args.ref + ".fai", "r") as f:
            for line in f:
                chr = line.strip().split("\t")[0]
                ctg_ls.append(chr)
    else:
        ref_dic = fasta_parser.read_sequence_dict(args.ref)    # 参考序列字典
        corrected_ref_dic = {}
        ctg_ls = [] # 储存进行SV consensus 的contig
        if args.all_chrs:
            for chr_id,ref_seq in ref_dic.items():
                corrected_ref_dic[chr_id] = ref_seq
                if len(ref_seq) > args.min_contig:
                    ctg_ls.append(chr_id)
        else:   # 只保留指定的contig
            ls = args.ctg_ls.strip(",").split(",")
            for chr_id in ls:
                ref_seq = ref_dic.get(chr_id)
                if not ref_seq: continue
                corrected_ref_dic[chr_id] = ref_seq
                if len(ref_seq) > args.min_contig:
                    ctg_ls.append(chr_id)
        logger.info("corrected ref include:{}".format(corrected_ref_dic.keys()))
        logger.info("You will SV consensus on:{}".format(",".join(ctg_ls)))
        fasta_parser.write_fasta_dict(corrected_ref_dic, corrected_ref)     # 写入新的参考基因组中
        # index
        faidx_cmd = ["samtools faidx", corrected_ref]
        subprocess.check_call(" ".join(faidx_cmd), shell=True)
        if len(ctg_ls) < 1:
            exit(0)

    ## 
    overwrite = args.overwrite
    re_run = args.re_run
    genome_size = convert_size(args.genome_size)
    if re_run: print("Will re run from {}".format(re_run))

######################################################### pipeline #########################################################
    
    ## step1 mapping 
    t1 = time.time()
    if args.find_mis:   # 
        print("Find mis mode, Skip mapping")
        samtools_out1 = args.bam
        overwrite = True
    else:
        minimap_out1 = os.path.join(step1_dir, "aln.sam")
        samtools_out1 = os.path.join(step1_dir, "aln.sorted.bam")
        if re_run == "step1": overwrite = True
        if os.path.isfile(samtools_out1) and not overwrite:
            logger.info("Skipped step1")
        else:
            logger.info("Run step1 mapping")
            minimap2_params_ls1 = config["step1"]["minimap2_params_ls"]
            # run_minimap(corrected_ref, args.fastq, args.data_type, args.threads, minimap_out1)
            try:
                generate_bam(corrected_ref, args.fastq, args.data_type, args.threads, samtools_out1, minimap2_params_ls1)
            except:
                raise Exception("Run mapping minimap2 samtools failed!!!")
            # try:
            #     run_minimap2(corrected_ref, args.fastq, args.data_type, args.threads, minimap_out1, minimap2_params_ls1)
            # except:
            #     raise Exception("Run mapping minimap2 failed!!!")
            # try:
            #     run_samtools(minimap_out1, samtools_out1, args.threads) # sam -> sorted.bam, and index
            #     os.remove(minimap_out1)
            # except:
            #     raise Exception("Run mapping samtools failed!!!")
            overwrite = True
    
    ## step2 Find candidate regions
    t2 = time.time()
    print("STEP1 cost: ", t2 - t1)

    keep_ls = []
    if args.correct or args.find_mis:   # filter out small contig
        bam_reader = pysam.AlignmentFile(samtools_out1, "rb", index_filename=samtools_out1+".bai")
        print("--------correct/find mis mode-------")
        min_contig = args.min_contig
        new_ctg_ls = []
        for chr in ctg_ls:
            if bam_reader.get_reference_length(chr) > min_contig:
                new_ctg_ls.append(chr)
            else:
                keep_ls.append(chr)
        print("We will corr/find mis on:{}".format(new_ctg_ls))
        ctg_ls = new_ctg_ls
    # find_mis mode
    if args.find_mis:
        # 1、mis find1
        run_find_pipe(args.ref, samtools_out1, ctg_ls, step2_dir, args.threads, config["step2"])
        # 2、some stats
        print("Start stats")
        ctg_ls = info_stats.simple_stats(args.ref, min_contig, step2_dir)
        info_stats.get_qv(args.ref, args.bam, step2_dir, ctg_ls, args.threads)
        return
    # assembly mode
    candidate_bed = os.path.join(step2_dir, "candidate.bed")
    if re_run == "step2": overwrite = True
    if os.path.isfile(candidate_bed) and not overwrite:
        logger.info("Skipped step2")
    else:
        logger.info("Run step2 Find candidate regions")
        # run_find_candidate_parallel(samtools_out1, ctg_ls, step2_dir, dp_file_dir, args.threads, args.min_clip_num, args.min_clip_len, args.dp_win_size, args.min_dep, args.max_dep, args.min_dep_reg)
        try:
            run_find_candidate2_parallel(corrected_ref, samtools_out1, ctg_ls, step2_dir, args.threads, config["step2"], args)
        except:
            raise Exception("Run step2 Failed!!!")
        overwrite = True

    
    ## step3 SVconsensus
    t3 = time.time()
    print("STEP2 cost: ", t3 - t2)
    SV_consensus_fasta = os.path.join(step3_dir, "consensus.fasta")
    SV_consensus_bed = os.path.join(step3_dir, "consensus.bed")
    merge_fasta = os.path.join(step3_dir, "merge.fasta")
    if re_run == "step3": overwrite = True
    if os.path.isfile(SV_consensus_fasta) and os.path.isfile(SV_consensus_bed) and os.path.isfile(merge_fasta) and not overwrite:
        logger.info("Skipped step3")
    else:
        bam_in = samtools_out1
        logger.info("Run step3 SV consensus")
        try:
            run_SVconsensus_parallel(step3_dir, bam_in, ctg_ls, candidate_bed, args.threads, args.fastq, args.data_type, corrected_ref, args.Nfill_size, genome_size, args.ex_unmapped_denovo, args.out_denovo, keep_ls, args, config["step3"])
        except:
            raise Exception("Run step3 Failed!!!")
        overwrite = True

    ## step4 polish
    ## correct mode
    if args.skip_denovo:
        fa_in = SV_consensus_fasta
    else:
        fa_in = merge_fasta
    t4 = time.time()    # start
    t5 = time.time()    # 
    t6 = time.time()
    print("STEP3 cost: ", t6 - t3)
    polish_out = os.path.join(step4_dir, "racon.fasta")
    if re_run == "step4": overwrite = True
    if os.path.isfile(polish_out) and not overwrite:
        logger.info("Skipped step4")
    else:
        logger.info("Run step4 polish")
        '''minimap_cmd = ["minimap2", "-ax", "map-"+data_type, params, target_fa, query_fa, "-t", str(threads), ">", minimap_out]'''
        '''racon {input.fastq} {output.sam} {input.asm} -t {threads} > {output.fasta}'''
        logger.info("Run mapping to consensus fasta")
        minimap_out2 = os.path.join(step4_dir, "aln.sam")
        minimap2_params_ls2 = config["step4"]["minimap2_params_ls"]
        racon_params = " ".join(config["step4"]["racon_params_ls"])
        
        ## 
        try:
            run_minimap2(fa_in, args.fastq, args.data_type, args.threads, minimap_out2, minimap2_params_ls2)
        except:
            raise Exception("Run mapping to consensus fasta failed!!!")
        t5 = time.time()

        logger.info("Run racon polish")
        ''' racon: -u, --include-unpolished   output unpolished target sequences'''
        polish_cmd = ["/usr/bin/time -v racon", racon_params, args.fastq, minimap_out2, fa_in, "-u", "-t", str(args.threads), ">", polish_out]
        logger.info("Running: %s", " ".join(polish_cmd))
        try:
            subprocess.check_call(" ".join(polish_cmd), shell=True)
        except:
            raise Exception("Run racon polish failed!!!")
        overwrite = True
        t6 = time.time()
    if args.skip_denovo:
        return
    ## step5 scaffolding
    print("STEP4 cost: ", t6 - t4)
    connect_file = os.path.join(step3_dir, "connect.txt")
    scaffold_out = os.path.join(step5_dir, "final.fasta")
    if re_run == "step5": overwrite = True
    if os.path.isfile(scaffold_out) and not overwrite:
        logger.info("Skipped all!!!")
    else:
        logger.info("Run step5 scaffolding")
        scaffold_by_rec.scaffold(polish_out, connect_file, scaffold_out, config["step5"])
        logger.info("Run scaffolding Done!!!")
    t7 = time.time()

    ## time stastic
    time_ls = [t1, t2, t3, t4, t5, t6, t7]
    print("Time stastic of the pipe: ")
    print("\t".join(["mapping raw reads,", "Find candidate regions,", "SVconsensus,", "mapping,", "polish,", "scaffolding"]))
    for i in range(len(time_ls)-1):
        print("t{}-t{}: {}s".format(i + 1, i + 2, time_ls[i + 1] - time_ls[i]))

if __name__ == "__main__":
    main()
    pass
