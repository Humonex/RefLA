#!/bin/bash
# simulate hifi reads

eval "$(conda shell.bash hook)"
conda activate pbsim3

fa_in=../Ref/template_ref.fasta
out_dir=fq
coverage=40
model=/public/home/hpc214712170/shixf/tools/pbsim3/data/QSHMM-RSII.model    # QSHMM-RSII.model
prefix=$out_dir/hifi_${coverage}x
pass_num=10     # default:10
threads=48

# ---------------WGS simulation---------------
echo "-------------------------Run pbsim3 for ${fa_in}, simulate hifi reads-------------------------"
mkdir -p $out_dir
# 1、run pbsim
pbsim --strategy wgs \
      --method qshmm \
      --qshmm $model \
      --depth $coverage \
      --genome $fa_in \
      --prefix $prefix \
      --pass-num $pass_num
echo "---------------Run pbsim3 done---------------"
# 
sam_ls=()   # ls
for file in $out_dir/*sam; do
    sam=$(echo $file | sed 's/.*\///; s/\..*//')
    # echo "$sam"
    sam_ls[${#sam_ls[*]}]=$sam     # 数组添加
done
echo ${sam_ls[*]}   # 输出数组所有

# 2、samtools view -bS sd_0001.sam > sd_0001.bam
echo "---------------sam->bam---------------"
time parallel --joblog samtools.log -j$threads "samtools view -bS $out_dir/{1}.sam > $out_dir/{1}.bam && rm $out_dir/{1}.sam" ::: ${sam_ls[@]}

# 3、ccs sd_0001.bam sd_0001.fastq.gz -j $threads
echo "---------------Run ccs---------------"
time parallel --joblog ccs.log -j$threads "ccs $out_dir/{1}.bam $out_dir/{1}.fastq.gz -j 2 && rm $out_dir/{1}.bam" ::: ${sam_ls[@]}

echo "---------------Merge fastq---------------"
cat $prefix*fastq.gz > $out_dir/merge.fastq
rm $prefix*fastq.gz
