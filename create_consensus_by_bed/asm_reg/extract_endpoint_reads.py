
import argparse
import os
import sys
import pysam


def read_reg_file(file):
    ls = []
    with open(file, "r") as f:
        for line in f:
            chr, start, end = line.strip().split()
            ls.append([chr, start, end])
    return ls

# def 
def fun(bam_in, reg_file, out_dir):
    '''samtools view -F 256 $bam ${chr}:${left}-${left2} | cut -f1 | sort | uniq > $out_dir/${chr}:${left}-${right}.left.bed'''
    '''
    out: 
    read_id     +/-
    '''
    bam_reader = pysam.AlignmentFile(bam_in, "rb")
    reg_ls = read_reg_file(reg_file)
    for reg in reg_ls:
        chr, start, end = reg
        start, end = int(start), int(end)
        out_file1 = os.path.join(out_dir, chr+":"+str(start)+"-"+str(end)+".left.bed")  # {chr}:${left}-${right}.left.bed
        out_file2 = os.path.join(out_dir, chr+":"+str(start)+"-"+str(end)+".right.bed")
        with open(out_file1, "w") as fo1, open(out_file2, "w") as fo2:
            read_set_l = set()
            read_set_r = set()
            l_0 = start
            l_1 = start + 2
            r_0 = end
            r_1 = end + 2
            for read in bam_reader.fetch(chr, start=l_0, end=l_1):
                # filter
                if read.is_secondary: continue
                
                if read.is_reverse: # 反链
                    read_set_l.add((read.query_name, "-"))
                else:
                    read_set_l.add((read.query_name, "+"))
            for read in bam_reader.fetch(chr, start=r_0, end=r_1):
                # filter
                if read.is_secondary: continue
                
                if read.is_reverse: # 反链
                    read_set_r.add((read.query_name, "-"))
                else:
                    read_set_r.add((read.query_name, "+"))
            ## 
            for read_id, strand in read_set_l:
                fo1.write("{}\t{}\n".format(read_id, strand))
            for read_id, strand in read_set_r:
                fo2.write("{}\t{}\n".format(read_id, strand))

if __name__ == "__main__":
    
    bam_in = sys.argv[1]
    reg_file = sys.argv[2]
    out_dir = os.path.abspath(sys.argv[3])
    fun(bam_in, reg_file, out_dir)
