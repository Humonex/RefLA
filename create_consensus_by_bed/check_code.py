# bed_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi_ctgs/NC_19/my_pipe/step3_SV_consensus/consensus.bed"
# with open(bed_in, "r")as f:
#     for line in f:
#         fields = line.strip().split("\t")
#         # print(fields[5])
#         if fields[3].startswith("reads"):
#             print(fields[4], fields[5])
import get_fasta_consensus2
import subprocess
fa_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi_ctgs/NC_19/my_pipe/step3_SV_consensus/denovo_asm/out.ctg.fa"
rec_bed = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi_ctgs/NC_19/my_pipe/step3_SV_consensus/consensus.bed"
# res = get_fasta_consensus2.get_frag_reads(fa_in, rec_bed)
cmd_2 = ["awk", "\'/^S/{print \">\"$2;print $3}\'", "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/create_consensus_by_bed/test_consensus/denovo_asm/hifiasm" + "/out.bp.p_ctg.gfa", ">", "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/create_consensus_by_bed/test_consensus/denovo_asm/hifiasm/out.p_ctg.fasta"]  # 注意使用转义字符
print(" ".join(cmd_2))
subprocess.check_call(" ".join(cmd_2), shell=True)