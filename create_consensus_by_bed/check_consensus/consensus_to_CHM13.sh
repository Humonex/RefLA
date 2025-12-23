#!/bin/bash
# map query_asm to target_asm
out_prefix="consensus_to_CHM13" # asm to asm
query_asm="/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/create_consensus_by_bed/NC_060930.1_hap1.fasta"
target_asm="/public/home/hpc214712170/shixf/projects/ref-guided-assembly/T2T_CHM13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
threads=48
##
time minimap2 -ax asm20 $target_asm $query_asm -t $threads > $out_prefix.sam && echo "** minimap2 done **"
# split-prefix for multi-part index
#time minimap2 -ax asm20 --split-prefix $out_prefix $target_asm $query_asm -t $threads > $out_prefix.sam && echo "** minimap2 done **"

samtools view -bS -@ $threads $out_prefix.sam > $out_prefix.bam
samtools sort -@ $threads $out_prefix.bam -o $out_prefix.sort.bam
samtools index -@ $threads $out_prefix.sort.bam

rm $out_prefix.sam $out_prefix.bam