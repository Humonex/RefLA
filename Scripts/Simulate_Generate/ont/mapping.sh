#!/bin/bash

ref=template_simu.fasta
# asm=
fq=fq/merge.fastq
work_dir=`pwd`
prefix=aln2simu
threads=40

# 1、asm mapping
# minimap2 -ax asm20 -t $threads $ref $asm | samtools sort -@ $threads -O BAM -o $work_dir/prefix.sort.bam
# 2、fastq mapping
minimap2 -ax map-ont ${ref} ${fq} -t $threads | samtools sort -@ $threads > $work_dir/$prefix.sort.bam && echo "mapping done!!!"
samtools index -@ $threads $work_dir/$prefix.sort.bam

