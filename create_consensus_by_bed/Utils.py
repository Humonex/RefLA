import logging
import subprocess
import os
import sys
import shutil
import time
import gzip
from typing import Any
from pathlib import Path
# import numpy as np

class Asm_Region():
    '''改类用于记录局部组装区域的相关信息'''
    def __init__(self, chr_id, start, end) -> None:
        self.chr_id = chr_id
        self.start = start
        self.end = end
        self.all_reads = None
        self.l_reads = None
        self.r_reads = None
        self.loca = None
    def add_all_reads(self, reads):
        self.all_reads = reads
    def add_l_reads(self, reads):
        self.l_reads = reads
    def add_r_reads(self, reads):
        self.r_reads = reads
    def add_loca(self, ctg_len):
        if self.end - self.start / ctg_len > 0.8: return "whole_ctg"
        elif self.start < 10000:return "left"
        elif self.end > ctg_len - 10000:return "right"
        else: return "mid"
    def get_reads(self, part):   # l、r、all
        if part == "all":
            return self.all_reads
        elif part == "left":
            return self.l_reads
        elif part == "right":
            return self.r_reads
        else: raise ValueError
    @staticmethod
    def write_reads(file, reads_out):
        with open(file, "w") as fo:
            for i in reads_out: fo.write(i + "\n")

 
'''Depth record for later'''    
class DepthRec():
    def __init__(self, chr_id, chr_len, dp_ls, win_size, block_size, whole_dp) -> None:
        self.chr_id = chr_id
        self.chr_len = chr_len
        self.dp_ls = dp_ls
        # self.dp_ls = np.array(dp_ls)
        self.dp_block_size = block_size
        self.dp_win_size = win_size     # win_size 应该提供初始计算win_size的大小
        self.block_dp_ls = []
        self.whole_dp = whole_dp
        ## 
        block_batch = self.dp_block_size // self.dp_win_size
        block_num = self.chr_len // self.dp_block_size + 1 if chr_len % self.dp_block_size > 0 else chr_len // self.dp_block_size  # 有多少个block
        for i in range(block_num):
            if i != block_num:
                block_dp = sum(self.dp_ls[i * block_batch : (i + 1) * (block_batch)]) // block_batch
            else:
                block_dp = sum(self.dp_ls[i * block_batch : len(self.dp_ls)]) // (len(self.dp_ls) - i * block_batch)
            self.block_dp_ls.append(block_dp)
    def get_reg_dp(self, reg):  # 实例方法
        start = reg[1]  # 100
        end = reg[2]    # 300
        l = start // self.dp_win_size
        r = end // self.dp_win_size
        return sum(self.dp_ls[l:r]) / (r - l)
    def get_block_dp(self, reg):
        start = reg[1]  # 1100
        end = reg[2]    # 1300
        l = start // self.dp_block_size     # 0
        r = end // self.dp_block_size + 1       # 
        return sum(self.block_dp_ls[l:r]) / (r - l)
    def get_block_high_dp(self, reg):
        start = reg[1]  # 1100
        end = reg[2]    # 1300
        l = start // self.dp_block_size     # 0
        r = end // self.dp_block_size + 1       # 
        # max(self.block_dp_ls[l:r])
        # max_dp = 0
        # for i in range(l, r):
        #     max_dp = max(self.block_dp_ls[i], max_dp)
        return max(self.block_dp_ls[l:r])
    @staticmethod
    def read_mosdepth_dp_file(dp_file):
        dp_ls = []
        with gzip.open(dp_file, "rt") as f:
            for line in f:
                fields = line.strip().split("\t")
                dp = float(fields[3])
                dp_ls.append(dp)
        return dp_ls
    

'''Record for consensus'''
class Record():
    def __init__(self, chr_id, start, end) -> None:
        self.chr_id = chr_id
        self.start = start
        self.end = end
        self.operation = None
        self.info = None
        self.patch_id = None
        self.loca = None    # 判断处于哪个部位
        self.ctg_len = None
        self.size = self.end - self.start
        self.old_start = start
        self.old_end = end
        self.skip = False
    def add_info(self, info):   # 实例方法
        self.info = info
    def add_operation(self, operation):
        self.operation = operation
    def add_patch_id(self, patch_id):   # 用于填充的asm_id read_id
        self.patch_id = patch_id
    def add_loca(self, loca):
        self.loca = loca
    ## 
    def get_op_ls(self):
        return self.operation.split(",")
    def get_patch_ls(self):
        return self.patch_id.split(",")
    def get_info_ls(self):
        return self.info.split(",")
    def print_info(self):
        for attr, value in vars(self).items():
            print(attr + ":", value)
    def __str__(self):
        return (f"Record(chr_id={self.chr_id}, start={self.start}, end={self.end}, "
                f"operation={self.operation}, patch_id={self.patch_id}, "
                f"loca={self.loca}, ctg_len={self.ctg_len}, size={self.size}, "
                f"old_start={self.old_start}, old_end={self.old_end}, skip={self.skip})")

    def __repr__(self):
        return self.__str__()
    ## 
    @staticmethod
    def write_record(record_ls, file_out):
        with open(file_out, "w") as fo:
            fo.write("#chr\tstart\tend\toperation\tinfo\tpatch_id\told_start\told_end\n")
            for rec in record_ls:
                fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(rec.chr_id, rec.start, rec.end, rec.operation, rec.info, rec.patch_id, rec.old_start, rec.old_end))
    @staticmethod
    def write_rec_short(record_ls, file_out):
        with open(file_out, "w") as fo:
            fo.write("#chr\tstart\tend\tloca\toperation\tpatch_id\treg_len\tpatch_len\n")
            for rec in record_ls:
                if rec.operation == "asm_patch" or rec.operation == "read_patch":
                    fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(rec.chr_id, rec.start, rec.end, rec.loca, rec.operation, rec.patch_id, rec.end - rec.start, len(rec.info)))
                else:
                    fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(rec.chr_id, rec.start, rec.end, rec.loca, rec.operation, rec.patch_id, rec.end - rec.start, None))
    @staticmethod
    def read_record(file_in):
        ls =[]
        with open(file_in, "r") as fin:
            for line in fin:
                if line.startswith("#"): continue
                chr_id, start, end, operation, info, patch_id = line.strip().split("\t")[:6]
                start, end = int(start), int(end)
                rec = Record(chr_id, start, end)
                rec.add_patch_id(patch_id)
                rec.add_info(info)
                rec.add_operation(operation)
                ls.append(rec)
        return ls


'''Connect info for scaffolding'''
class Connect_info():
    def __init__(self, chr_id, connect_ls:list, gap_ls:list) -> None:
        self.chr_id = chr_id
        self.connect_ls = connect_ls
        self.gap_ls = gap_ls

    @staticmethod
    def read_connect_info(file):
        ls = []
        with open(file, "r") as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < 2:continue
                chr_id = fields[0]
                connect_ls = fields[1].split(",")
                if len(connect_ls) > 1: gap_ls = fields[2].split(",")
                else: gap_ls = []
                ls.append(Connect_info(chr_id, connect_ls, gap_ls))
        return ls

    @staticmethod
    def write_connect_info(connect_info_ls, file):
        ## tab键分割，第一列记录chr_id，第二列记录connect_ls，第三列记录gap_ls
        ## chr1 chr1_0,chr_1,chr_2  100,300
        with open(file, "w") as f:
            for connect_info in connect_info_ls:
                chr_id, connect_ls, gap_ls = connect_info.chr_id, connect_info.connect_ls, connect_info.gap_ls
                if len(connect_info.gap_ls) > 0:
                    f.write("{}\t{}\t{}\n".format(chr_id, ",".join(connect_ls), ",".join(gap_ls)))
                else:
                    f.write("{}\t{}\t{}\n".format(chr_id, ",".join(connect_ls), "."))
    def __str__(self):
        return f'Connect_info(chr_id={self.chr_id}, connect_ls={self.connect_ls}, gap_ls={self.gap_ls})'


def write_connect_info(connect_info_ls:list, file):
    ## tab键分割，第一列记录chr_id，第二列记录connect_ls，第三列记录gap_ls
    ## chr1 chr1_0,chr_1,chr_2  100,300
    with open(file, "w") as f:
        for connect_info in connect_info_ls:
            chr_id, connect_ls, gap_ls = connect_info.chr_id, connect_info.connect_ls, connect_info.gap_ls
            f.write("{}\t{}\t{}\n".format(chr_id, ",".join(connect_ls), ",".join(gap_ls)))

def read_connect_info(file):
    ls = []
    with open(file, "r") as f:
        for line in f:
            chr_id, connect_ls, gap_ls = line.strip().split("\t")
            ls.append(Connect_info(chr_id, connect_ls, gap_ls))
    return ls

def run_cmd_ls(cmd_ls):
    t0 = time.time()
    print("Run: {}".format(" ".join(cmd_ls)))
    subprocess.check_call(" ".join(cmd_ls), shell=True)
    print("Run: {} finished, cost {}s".format(" ".join(cmd_ls), time.time() - t0))

def make_dir(dir):
    if not os.path.isdir(dir):
        os.makedirs(dir)
def reg_to_id(reg):
    return reg[0] + ":" + str(reg[1]) + "-" + str(reg[2])

def id_to_reg(reg_str:str):
    ctg = reg_str.split(":")[0]
    start, end = reg_str.split(":")[1].split("-")
    start, end = int(start), int(end)
    return [ctg, start, end]


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

logger = logging.getLogger()
def _enable_logging(log_file, debug, overwrite):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    logger = logging.getLogger()
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
    # return logger

def write_rec_to_bed(bed_file, rec_ls):
    with open(bed_file, "w") as f:
        for rec in rec_ls:
            f.write("{}\t{}\t{}\n".format(rec.chr_id, rec.start, rec.end))

def write_reg_to_bed(bed_file, reg_ls):
    with open(bed_file, "w") as f:
        for reg in reg_ls:
            f.write("{}\t{}\t{}\n".format(reg[0], reg[1], reg[2]))

def get_unmapped_reads(threads, bam_in, f_out, f_type):   # f_type: fastq fasta
    t1 = time.time()
    cmd_ls = ["samtools", "view", "-h", "-f4", "-@", str(threads), bam_in, "|", "samtools", f_type, "-@", str(threads), ">", f_out]
    logger.info("Running: %s", " ".join(cmd_ls))
    subprocess.check_call(" ".join(cmd_ls), shell=True)
    print("Get unmapped reads finished, cost {}s".format(time.time() - t1))

#### Tools for denovo
def run_hifiasm_hom(fq_in, work_dir, threads, data_type):  # Assemble inbred/homozygous genomes (-l0 disables duplication purging)
    '''hifiasm -o $work_dir/out -t $threads -l0 $fq_path  # (-l0 disables duplication purging)'''
    ''' awk '/^S/{print ">"$2;print $3}' $work_dir/out.bp.p_ctg.gfa > $work_dir/out.p_ctg.fasta '''
    if data_type != "hifi": raise("hifiasm only for hifi data")   # 
    if not os.path.isdir(work_dir):os.makedirs(work_dir)
    hifiasm_fa_out = os.path.join(work_dir, "out.p_ctg.fasta")
    cmd_1 = ["hifiasm", "-o", work_dir + "/out", "-t", str(threads), "-l0", fq_in]
    cmd_2 = ["awk", "\'/^S/{print \">\"$2;print $3}\'", work_dir + "/out.bp.p_ctg.gfa", ">", hifiasm_fa_out]
    logger.info("Running: %s", " ".join(cmd_1))
    subprocess.check_call(" ".join(cmd_1), shell=True)
    logger.info("Running: %s", " ".join(cmd_2))
    subprocess.check_call(" ".join(cmd_2), shell=True)
    logger.info("hifiasm running Done!!!")
    return hifiasm_fa_out

def run_hifiasm(fq_in, work_dir, threads, data_type, hifiasm_options):  # Assemble inbred/homozygous genomes (-l0 disables duplication purging)
    '''hifiasm -o $work_dir/out -t $threads $fq_path  # '''
    '''hifiasm -o $work_dir/out -t $threads -l0 $fq_path  # (-l0 disables duplication purging)'''
    ''' awk '/^S/{print ">"$2;print $3}' $work_dir/out.bp.p_ctg.gfa > $work_dir/out.p_ctg.fasta '''
    if data_type != "hifi": raise("hifiasm only for hifi data")   # 
    if os.path.exists(work_dir):shutil.rmtree(work_dir)     # 目录存在要先删除，注意目录要给对，最好不是已有目录, hifiasm会受到原有数据的影响
    if not os.path.exists(work_dir):os.makedirs(work_dir)
    hifiasm_fa_out = os.path.join(work_dir, "out.p_ctg.fasta")
    cmd_1 = ["hifiasm", "-o", work_dir + "/out", "-t", str(threads), fq_in, hifiasm_options]  # hifiasm组装，添加
    cmd_2 = ["awk", "\'/^S/{print \">\"$2;print $3}\'", work_dir + "/out.bp.p_ctg.gfa", ">", hifiasm_fa_out]
    logger.info("Running: %s", " ".join(cmd_1))
    subprocess.check_call(" ".join(cmd_1), shell=True)
    logger.info("Running: %s", " ".join(cmd_2))
    subprocess.check_call(" ".join(cmd_2), shell=True)
    logger.info("hifiasm running Done!!!")
    return hifiasm_fa_out

def run_flye(fq_in, work_dir, threads, genome_size, data_type, flye_options):
    print("Choose flye for denovo")
    if not os.path.isdir(work_dir):os.makedirs(work_dir)
    flye_fa_out = os.path.join(work_dir, "assembly.fasta")
    if data_type == "ont":
        cmd = ["flye", "--nano-raw", fq_in, "-o", work_dir, "-t", str(threads), flye_options]  # flye组装
    elif data_type == "hifi":
        cmd = ["flye", "--pacbio-hifi", fq_in, "-o", work_dir, "-t", str(threads), flye_options]
    else:
        raise("Error data type: {}".format(data_type))
    logger.info("Running: %s", " ".join(cmd))
    subprocess.check_call(" ".join(cmd), shell=True)
    logger.info("flye running Done!!!")
    return flye_fa_out
    
def run_wtdbg2(fq_in, work_dir, threads, genome_size, data_type, wtdbg2_options):
    print("Choose wtdbg2 for denovo")
    if not os.path.isdir(work_dir):os.makedirs(work_dir)
    reads_fq = fq_in
    out_prefix = work_dir+"/out"
    wtdbg2_fa_out = out_prefix+".ctg.fa"    # 
    ## preset choose
    if data_type == "ont": ## wtdbg2选择不同的模式
        if genome_size < 1000000000:preset_type = "preset2" # ont <1G
        else:preset_type = "preset3" # ont >=1G
    elif data_type == "hifi":preset_type = "preset4" # hifi
    else:raise("wtdbg2 fastq type not support")
    ### cmd set 
    # asm_cmd1 = "time wtdbg2 -x preset2 -t " + str(threads) + " -i " + reads_fq +" -fo " + out_prefix + " > wtdbg2.log"
    # asm_cmd2 = "time wtpoa-cns -t " + str(threads) + " -i " + out_prefix+".ctg.lay.gz" + " -fo " + out_prefix+".ctg.fa" + " > wtpoa-cns.log"
    asm_cmd1 = ["time wtdbg2", "-x", preset_type, "-t", str(threads), "-i ", reads_fq, " -fo ",  out_prefix]
    asm_cmd2 = ["time wtpoa-cns", "-t", str(threads), "-i", out_prefix+".ctg.lay.gz", "-fo", out_prefix+".ctg.fa"]
    for cmd in [asm_cmd1, asm_cmd2]:
        logger.info("Run: {}".format(cmd))
        subprocess.check_call(" ".join(cmd), shell=True)
    logger.info("Run wtdbg2 Done !!!")
    return wtdbg2_fa_out

def run_shasta(fq_in, work_dir, threads, data_type, shasta_options):    # Haploid assembly
    '''shasta will first remove the dir exist'''
    ''' 30X 
    shasta-Linux-0.11.0 --input input.fasta --config Nanopore-May2022 '''
    ''' 30X 
    --config Nanopore-Human-SingleFlowcell-May2022 '''
    print("Choose shasta for denovo")
    shasta_fa_out = os.path.join(work_dir, "Assembly.fasta")
    if os.path.exists(work_dir):shutil.rmtree(work_dir)     # 目录存在要先删除，注意目录要给对，最好不是已有目录
    if not os.path.exists(work_dir):os.makedirs(work_dir)
    # cmd = ["time shasta-Linux-0.11.0", "--input", fq_in, "--config", shasta_config, "--threads", str(threads), "--assemblyDirectory", work_dir]
    cmd = ["time shasta-Linux-0.11.0", "--input", fq_in, shasta_options, "--threads", str(threads), "--assemblyDirectory", work_dir]
    logger.info("Run: {}".format(cmd))
    subprocess.check_call(" ".join(cmd), shell=True)
    logger.info(" Run shasta Done !!!")
    return shasta_fa_out

def run_raven(fq_in, work_dir, threads, data_type, raven_options):
    '''usage: raven [options ...] <sequences> [<sequences> ...]'''
    pass

def run_miniasm(fq_in, work_dir, threads, data_type, miniasm_options):
    ''' from NieFan
    minimap2 -x ava-ont -t32 all.fastq all.fastq | gzip -1 > reads.paf.gz
    miniasm -f all.fastq reads.paf.gz > $genomeName.gfa
    awk '/^S/{print ">"$2"\n"$3}' $genomeName.gfa > $genomeName.fasta
    '''
    '''
    # Overlap for PacBio reads (or use "-x ava-ont" for nanopore read overlapping)
    minimap2/minimap2 -x ava-pb -t8 pb-reads.fq pb-reads.fq | gzip -1 > reads.paf.gz
    # Layout
    miniasm/miniasm -f reads.fq reads.paf.gz > reads.gfa
    '''
    
    if not os.path.isdir(work_dir):os.makedirs(work_dir)
    # step1
    minimap2_out = os.path.join(work_dir, "reads.paf.gz")
    if data_type == "ont":
        cmd1 = ["minimap2 -x ava-ont", "-t", str(threads), fq_in, fq_in, "| gzip -1 >", minimap2_out]
        logger.info("Run: {}".format(cmd1))
        subprocess.check_call(" ".join(cmd1), shell=True)
    elif data_type == "hifi":
        cmd1 = ["minimap2 -x ava-pb", "-t", str(threads), fq_in, fq_in, "| gzip -1 >", minimap2_out]
        logger.info("Run: {}".format(cmd1))
        subprocess.check_call(" ".join(cmd1), shell=True)
    # step2
    miniasm_out = os.path.join(work_dir, "reads.gfa")
    cmd2 = ["miniasm -f", fq_in, minimap2_out, ">", miniasm_out]
    logger.info("Run: {}".format(cmd2))
    subprocess.check_call(" ".join(cmd2), shell=True)
    # step3
    # cmd_2 = ["awk", "\'/^S/{print \">\"$2;print $3}\'", work_dir + "/out.bp.p_ctg.gfa", ">", hifiasm_fa_out]
    miniasm_fa_out = os.path.join(work_dir, "asm.fa")
    cmd3 = ["awk", "\'/^S/{print \">\"$2;print $3}\'", miniasm_out, ">", miniasm_fa_out]
    logger.info("Run: {}".format(cmd3))
    subprocess.check_call(" ".join(cmd3), shell=True)
    return miniasm_fa_out


'''Run local denovo assembly'''
def Run_for_denovo(fq_in, work_dir, threads, genome_size, data_type, config):
    t0 = time.time()
    denovo_asm_params = config["denovo_asm"]
    tool = denovo_asm_params[data_type]  # 
    # tool_config = config["denovo_asm"]["tool_config"]
    tool_option_ls = denovo_asm_params.get(tool + "_option_ls", [])
    tool_options = " ".join(tool_option_ls)     # get the tool options of the specify tool
    print("Running {} options:{}".format(tool, tool_options))
    logger.info("Choose {} for local denovo".format(tool))
    try:
        # if tool == "hifiasm": denovo_asm_out = run_hifiasm_hom(fq_in, work_dir + "/hifiasm", threads, data_type)
        if tool == "hifiasm": denovo_asm_out = run_hifiasm(fq_in, work_dir + "/hifiasm", threads, data_type, tool_options)
        elif tool == "shasta": denovo_asm_out = run_shasta(fq_in, work_dir + "/shasta", threads, data_type, tool_options)
        elif tool == "wtdbg2": denovo_asm_out = run_wtdbg2(fq_in, work_dir + "/wtdbg2", threads, genome_size, data_type, tool_options)
        elif tool == "flye": denovo_asm_out = run_flye(fq_in, work_dir + "/flye", threads, genome_size, data_type, tool_options)
        elif tool == "miniasm": denovo_asm_out = run_miniasm(fq_in, work_dir + "/miniasm", threads, data_type, tool_options)
        else: raise("Error tool")
    except Exception as e:
        print("组装失败: ", e)
        failed_asm = work_dir + "failed.fa"
        failed_asm_path = Path(failed_asm)
        failed_asm_path.touch()
        return failed_asm
    print("------------------------------------Assembly cost:{}s------------------------------------".format(time.time() - t0))
    return denovo_asm_out
def Run_for_denovo2(fq_in, work_dir, threads, genome_size, data_type, config):
    denovo_asm_params = config["denovo_asm"]
    tool = denovo_asm_params[data_type]  # 
    # tool_config = config["denovo_asm"]["tool_config"]
    tool_option_ls = denovo_asm_params.get(tool + "_option_ls", [])
    tool_options = " ".join(tool_option_ls)     # get the tool options of the specify tool
    print("Running {} options:{}".format(tool, tool_options))
    logger.info("Choose {} for local denovo".format(tool))
    # if tool == "hifiasm": denovo_asm_out = run_hifiasm_hom(fq_in, work_dir + "/hifiasm", threads, data_type)
    if tool == "hifiasm": denovo_asm_out = run_hifiasm(fq_in, work_dir + "/hifiasm", threads, data_type, tool_options)
    elif tool == "shasta": denovo_asm_out = run_shasta(fq_in, work_dir + "/shasta", threads, data_type, tool_options)
    elif tool == "wtdbg2": denovo_asm_out = run_wtdbg2(fq_in, work_dir + "/wtdbg2", threads, genome_size, data_type, tool_options)
    elif tool == "flye": denovo_asm_out = run_flye(fq_in, work_dir + "/flye", threads, genome_size, data_type, tool_options)
    else: raise("Error tool")
    return denovo_asm_out

def cal_genome_siez():
    # bed = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/yeast_ont/my_pipe/step3_SV_consensus/candidate_op/candidate_op.bed"
    bed = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/yeast_ont/my_pipe/step3_SV_consensus/candidate_op/merge0/merge0.bed"
    ctg_len_dic = {}
    with open(bed, "r") as f:
        genome_size = 0
        ls = []
        for line in f:
            fields = line.strip().split("\t")
            ctg, start, end, op = fields[:4]
            start, end = int(start), int(end)
            if op.endswith("asm"):
                ls.append(end - start)
                ctg_len_dic[ctg] = ctg_len_dic.get(ctg, 0) + end - start
                # if ctg in ctg_len_dic:
                #     ctg_len_dic[ctg] += end-start
                # else:
                #     ctg_len_dic[ctg] = 0
                #     ctg_len_dic[ctg] += end - start
                if end - start < 100000:
                    genome_size += (end - start) * 2
                else:
                    genome_size += (end - start) * 1.5
    print(genome_size)
    print(sorted(ls))
    print(ctg_len_dic)
def stats_asm_patch():
    bed = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/yeast_ont/my_pipe/step3_SV_consensus/consensus.bed"
    ls = []
    with open(bed, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            ctg, start, end, op = fields[:4]
            start, end = int(start), int(end)
            if op == "asm_patch":
                ls.append(end -start)
    print(sorted(ls))

if __name__ == "__main__":
    deprec_dic = {}
    win_size = 100
    block_size = 5000
    whole_dp = 31.98
    ctg_len_dic = {'NC_000932.1':154478, 'NC_003070.9':30427671, 'NC_003071.7':19698289, 'NC_003074.8':23459830, 'NC_003075.7':18585056, 'NC_003076.8':26975502, 'NC_037304.1':367808}
    ctg_ls = ['NC_000932.1', 'NC_003070.9', 'NC_003071.7', 'NC_003074.8', 'NC_003075.7', 'NC_003076.8', 'NC_037304.1']
    dir = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/thaliana_ont/my_pipe2/step2_candidate_regions/depths"
    for ctg in ctg_ls:
        chr_len = ctg_len_dic[ctg]
        dp_ls = DepthRec.read_mosdepth_dp_file(dir + "/" + ctg + "/" + ctg + ".regions.bed.gz")
        deprec = DepthRec(ctg, chr_len, dp_ls, win_size, block_size, whole_dp)
        deprec_dic[ctg] = deprec
        # print(chr_len, deprec.block_dp_ls)
        # break
    # print(deprec_dic)
    # chr_id = "NC_003071.7"
    # chr_len = 19698289
    # deprec = DepthRec(chr_id, chr_len, dp_ls, win_size, block_size, whole_dp)
    reg_ls =  [['NC_003070.9', 0, 52200], ['NC_003070.9', 2443500, 2476500], ['NC_003070.9', 6180500, 6215500], ['NC_003070.9', 14313700, 14406000], ['NC_003070.9', 14455500, 14901300], ['NC_003070.9', 15004700, 15467100], ['NC_003070.9', 20927500, 20957000], ['NC_003070.9', 30372000, 30427671], ['NC_003071.7', 0, 88900], ['NC_003071.7', 880000, 903000], ['NC_003071.7', 3242100, 3425700], ['NC_003071.7', 3584500, 3676300], ['NC_003071.7', 5166700, 5238500], ['NC_003071.7', 19686000, 19698289], ['NC_003074.8', 0, 96400], ['NC_003074.8', 13554000, 13895300], ['NC_003074.8', 14128500, 14268700], ['NC_003074.8', 23409000, 23459830], ['NC_003075.7', 0, 12500], ['NC_003075.7', 2892700, 3164300], ['NC_003075.7', 3909500, 4048000], ['NC_003075.7', 5012700, 5111000], ['NC_003075.7', 9246000, 9346500], ['NC_003075.7', 18537500, 18585056], ['NC_003076.8', 0, 27900], ['NC_003076.8', 11150900, 11262900], ['NC_003076.8', 11638700, 12091100], ['NC_003076.8', 12770700, 12908500], ['NC_003076.8', 14389000, 14538100], ['NC_003076.8', 15625700, 15716000], ['NC_003076.8', 26953400, 26975502], ['NC_037304.1', 0, 88500], ['NC_037304.1', 104500, 291000], ['NC_037304.1', 351000, 367808], ['NC_000932.1', 0, 154478]]
    # print(deprec.get_block_dp(["NC_003071.7", 0, 88900]))
    for reg in reg_ls:
        if deprec_dic[reg[0]].get_block_dp(reg) / whole_dp > 2:
            print("{}dp:{}".format(reg, str(deprec_dic[reg[0]].get_block_dp(reg))))
    print("\nMax dp")
    norm_ls = []
    for reg in reg_ls:
        if deprec_dic[reg[0]].get_block_high_dp(reg) / whole_dp > 2:
            print("{}dp:{}".format(reg, str(deprec_dic[reg[0]].get_block_high_dp(reg))))
        norm_ls.append(deprec_dic[reg[0]].get_block_high_dp(reg) / whole_dp)
    print(sorted(norm_ls))
    # print(os.path.dirname("/public/home/hpc214712170/Test/tests/yeast_ont/my_pipe3/step3_SV_consensus"))
    # fq_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/yeast_ont/my_pipe/step3_SV_consensus/denovo_asm/denovo.fastq"
    # work_dir = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/yeast_ont/my_pipe/step3_SV_consensus/denovo_asm/shasta"
    # fq_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_ont/my_pipe/step3_SV_consensus/denovo_asm/denovo.fastq"
    # work_dir = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_ont/my_pipe/step3_SV_consensus/denovo_asm/shasta"
    # threads = 10
    # data_type = "ont"
    # shasta_config = "Nanopore-May2022"
    # run_shasta(fq_in, work_dir, threads, data_type, shasta_config)
    # cal_genome_siez()
    # stats_asm_patch()
    # op = " ".join([])
    # cmd = ["ni", op, "nia"]
    # print(" ".join(cmd))
    # rec = Record("1", 0, 100)
    # rec.add_operation("patch")
    # rec.add_info("AGCT")
    
    pass