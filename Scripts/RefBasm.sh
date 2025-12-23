#!/bin/bash
source activate assembly
script_dir="/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/pipe_main.py"
quast_route="/public/home/hpc214712170/shixf/tools/quast-5.0.2"
work_dir=/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/bench/chm13_hifi/my_pipe
ref1=/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/chm13/REF/ref1/GCF_000001405.40_GRCh38.p14_genomic.fna  # ref for assembly
ref2=/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/chm13/REF/ref2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna  # ref for eval
fq_path=/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/chm13/hifi/chm13_hifi.fastq  # fastq in
large=yes   # no or yes
data_type=hifi	# data type: ont hifi
genome_size=3G	# 10m 100m 1G  1000b
ctg_ls=NC_000001.11,NC_000002.12,NC_000003.12,NC_000004.12,NC_000005.10,NC_000006.12,NC_000007.14,NC_000008.11,NC_000009.12,NC_000010.11,NC_000011.10,NC_000012.12,NC_000013.11,NC_000014.9,NC_000015.10,NC_000016.10,NC_000017.11,NC_000018.10,NC_000019.10,NC_000020.11,NC_000021.9,NC_000022.11,NC_000023.11,NC_000024.10     # ctg_ls kept to perform assembly
asm=$work_dir/step4_polish/racon.fasta     # 
threads=40    # 
config=Config.yaml      # config file path
# /usr/bin/time -v python $script_dir -t $threads -a --work-dir $work_dir --ref $ref1 --fastq $fq_path --data-type $data_type -g $genome_size --config $config
/usr/bin/time -v python $script_dir -t $threads --ctg-ls $ctg_ls --work-dir $work_dir --ref $ref1 --fastq $fq_path --data-type $data_type -g $genome_size --config $config

## 大于100m选large==yes
# bash $work_dir/run_quast.sh $asm $Reference $work_dir/quast_eval $threads $large
sh $work_dir/run_quast.sh

