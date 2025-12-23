#!/bin/bash

out_dir=/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/thaliana_ont/tests/my_pipe4/step3_SV_consensus/endpoint_reads
reg_file=/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/create_consensus_by_bed/asm_reg/reg.bed
bam=/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/thaliana_ont/tests/my_pipe4/step1_mapping/aln.sorted.bam

while read line
do
    array=(`echo $line | tr '\t' ' '`)
    for i in "${!array[@]}"; do
        echo "$i=>${array[i]}"
    done
    chr=${array[0]}
    left=${array[1]}
    right=${array[2]}
    left2=$((left+2))
    right2=$((right+2))
    echo $chr $left2 $right2
    samtools view -F 256 $bam ${chr}:${left}-${left2} | cut -f1 | sort | uniq > $out_dir/${chr}:${left}-${right}.left.bed
    samtools view -F 256 $bam ${chr}:${right}-${right2} | cut -f1 | sort | uniq > $out_dir/${chr}:${left}-${right}.right.bed
done < $reg_file
