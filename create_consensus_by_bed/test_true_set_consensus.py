import fasta_parser
import os
import pysam
from multiprocessing import Pool
from collections import defaultdict
from get_fasta_consensus2 import SV_consensus_on_ref
class Record():
    def __init__(self, chr_id, start, end) -> None:
        self.chr_id = chr_id
        self.start = start
        self.end = end
        self.operation = None
        self.info = None
        self.patch_id = None
    def add_info(self, info):
        self.info = info
    def add_operation(self, operation):
        self.operation = operation
    def add_patch_id(self, patch_id):   # 用于填充的asm_id read_id
        self.patch_id = patch_id
def get_candidate_op_dic():
    candidate_op_dic = defaultdict(list)
    bed_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi/my_pipe/step3_SV_consensus/candidate_op/candidate_op.bed"
    with open(bed_in, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            ctg, start, end, op, info, patch_id = fields[0], int(fields[1]), int(fields[2]), fields[3], fields[4], fields[5]
            rec = Record(ctg, start, end)
            rec.add_info(info)
            rec.add_operation(op)
            rec.add_patch_id(patch_id)
            candidate_op_dic[ctg].append(rec)
    return candidate_op_dic
def test():
    candidate_op_dic = get_candidate_op_dic()
    threads = 20
    Nfill_size = 20
    bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi/my_pipe/step1_mapping/aln.sorted.bam"
    reference_fn = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/chm13/REF/ref1/new_GRCH38.fna"

    SV_consensus_dir = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/create_consensus_by_bed/test_consensus_with_true_set"
    if not os.path.exists(SV_consensus_dir):
        os.makedirs(SV_consensus_dir)
    denovo_fa = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/chm13/REF/ref2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"    # denovo asm out
    SVconsensus_bed_out = SV_consensus_dir + "/consensus.bed"
    consensus_fasta_out = SV_consensus_dir + "/consensus.fasta"
    asm_to_ref_bam = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/chm13/REF/ref2/bam/chm13_to_GRCH38.sort.bam"
    # ref_dic = parse_fasta(reference_fn)   # 
    # asm_fa_dic = parse_fasta(denovo_fa)
    ref_dic = fasta_parser.read_sequence_dict(reference_fn)
    asm_fa_dic = fasta_parser.read_sequence_dict(denovo_fa)
    consensus_fa_dic = {}
    print(" Parse fasta Done !!!")
    bam_reader = pysam.AlignmentFile(bam_in, "rb")
    all_chrs = bam_reader.references
    # run
    print("Run SV_consensus_on_ref")
    pool = Pool(processes=threads)
    results = [pool.apply_async(SV_consensus_on_ref, args=(candidate_op_dic[ctg], ref_dic[ctg], asm_fa_dic, asm_to_ref_bam, Nfill_size)) for ctg in all_chrs]
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    with open(consensus_fasta_out, "w") as fo, open(SVconsensus_bed_out, "w") as fo2:     # 非并行写文件
        for i, res in enumerate(results):
            final_rec_ls, new_seq_ls = res.get()
            ctg = all_chrs[i]
            ##
            if len(new_seq_ls) == 1:
                # fasta_write(ctg, new_seq_ls[0], fo)
                consensus_fa_dic[ctg] = new_seq_ls[0]
            else:
                print("multiple seqs of {} !!!".format(ctg))
                for i in range(len(new_seq_ls)):
                    # fasta_write(ctg+str(i), new_seq_ls[i], fo)
                    consensus_fa_dic[ctg+str(i)] = new_seq_ls[i]
            ## write reccord bed
            for rec in final_rec_ls:
                fo2.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(rec.chr_id, rec.start, rec.end, rec.operation, rec.info, rec.patch_id))
            # break
    fasta_parser.write_fasta_dict(consensus_fa_dic, consensus_fasta_out)
test()