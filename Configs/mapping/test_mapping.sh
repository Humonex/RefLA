#!/bin/bash

ref=/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/yeast/REF/ref1/Y12.fasta
# fq=/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/yeast/ont/CRR198362/NC_001133.9.fastq
fq=/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/yeast/ont/CRR198362/CRR198362_40X.fastq.gz
out_dir=mapping
time minimap2 -t 40 -ax map-ont $ref $fq > $out_dir/aln1.sam
time minimap2 -t 40 -ax map-ont -k 25 -w 16 -f 0.0004 $ref $fq > $out_dir/aln2.sam
time minimap2 -t 40 -a -k 25 -w 16 -f 0.0004 $ref $fq > $out_dir/aln3.sam
time minimap2 -t 40 -k 25 -w 16 -f 0.0004 $ref $fq > $out_dir/aln4.sam
time minimap2 -t 40 -x map-ont -k 25 -w 16 -f 0.0004 $ref $fq > $out_dir/aln5.sam
time minimap2 -t 40 -x map-ont -k 25 -w 16 -f 0.0004 $ref $fq > $out_dir/aln5.sam
-w30 -k19 -I 50g, -K 50G