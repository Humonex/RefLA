#!/bin/bash
# simulate ont reads
eval "$(conda shell.bash hook)"
conda activate pbsim3

fa_in=../Ref/template_ref.fasta
out_dir=fq
coverage=40
model=/public/home/hpc214712170/shixf/tools/pbsim3/data/QSHMM-ONT-HQ.model # QSHMM-ONT-HQ.model/QSHMM-ONT.model
prefix=$out_dir/ont_${coverage}x

# WGS simulation
echo "-------------------------Run pbsim3 for ${fa_in}-------------------------"
mkdir -p $out_dir
pbsim --strategy wgs --method qshmm --qshmm $model --depth $coverage --genome $fa_in --prefix $prefix
echo "---------------Run pbsim3 done---------------"

# 
echo "---------------Merge fastq---------------"
cat $prefix*fastq > $out_dir/merge.fastq
rm $prefix*fastq
