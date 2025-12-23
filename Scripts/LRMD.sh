#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate assembly
script_dir=""	# path to the pipe_main.py
work_dir=my_pipe
ref=  # ref for assembly
fq_path=  # fastq in
data_type=hifi	# data type: ont hifi
bam=
genome_size=3G	# 10m 100m 1G  1000b
ctg_ls=     # ctg_ls kept to perform assembly
threads=    # 
config=Config.yaml      # config file path
# time python $script_dir -t $threads -a --work-dir $work_dir --ref $ref1 --fastq $fq_path --data-type $data_type -g $genome_size --config $config
/usr/bin/time -v python $script_dir -t $threads -a --work-dir $work_dir --ref $ref --fastq $fq_path --data-type $data_type -g $genome_size --config $config --find_mis --bam $bam

