import gzip
def select_reads_from_names(fastq_in, fastq_out, readset:set):
    # reads_names = get_reads_names(bam)
    # readset = set(reads_names) # 去重一下
    if fastq_in.endswith("gz"):
        fin = gzip.open(fastq_in, "rt")
    else:
        fin = open(fastq_in, "r")
    fout = open(fastq_out,'w')
    for line in fin:
        if line.startswith('@'):
            if line[1:].strip().split()[0] in readset:
                fout.write(line)
                fout.write(fin.readline())
                fout.write(fin.readline())
                fout.write(fin.readline())
                readset.remove(line[1:].strip().split()[0])     # 从中去除
            else:
                fin.readline()
                fin.readline()
                fin.readline()
    fin.close()
    fout.close()