import sys
import pysam


aln = sys.argv[1]
fai = sys.argv[2]
ctg = sys.argv[3]
pos1 = int(sys.argv[4])
pos2 = int(sys.argv[5])


def convert_reference_pos_to_raw_pos(raw_to_reference_file, candidate_pos):
    candidate_pos = set(candidate_pos)
    raw_ref_pos_map={}
    samfile = pysam.AlignmentFile(raw_to_reference_file, "rb")
    read_length = 0
    for read in samfile.fetch():
        ref_pos = read.reference_start
        query_pos = 0
        read_length = read.query_length
        for (ct,cl) in read.cigartuples:
            if ct==0:
                for i in range(cl):
                    query_pos+=1
                    ref_pos+=1
                    if ref_pos in candidate_pos:
                        raw_ref_pos_map[ref_pos] = query_pos
            elif ct==1:
                query_pos+=cl
            elif ct==2:
                for i in range(cl):
                    ref_pos+=1
                    if ref_pos in candidate_pos:
                        raw_ref_pos_map[ref_pos] = query_pos
            elif ct==4:
                query_pos+=cl
            else:
                continue
    return raw_ref_pos_map,read_length


def annotate(bam,ctgname,pos1,pos2, read_lengths):
    snp_sites = [pos1,pos2]
    samfile = pysam.AlignmentFile(bam, "rb")
    for read in samfile.fetch(ctgname):
        rname = read.query_name
        ref_pos = read.reference_start  # 0-based
        query_pos = 0
        read_length = read_lengths[rname]
        strand = read.is_reverse
        read_snp_sites = [] # 0-based
        if not strand:
            for (ct,cl) in read.cigartuples:
                if ct==0:
                    for i in range(cl):
                        query_pos+=1
                        ref_pos+=1
                        if ref_pos in snp_sites:
                            # read_snp_sites.append(str(ref_pos)+':'+str(query_pos))
                            read_snp_sites.append(str(query_pos))
                elif ct==1:
                    query_pos+=cl
                elif ct==2:
                    for i in range(cl):
                        ref_pos+=1
                        if ref_pos in snp_sites:
                            # read_snp_sites.append(str(ref_pos)+':'+str(query_pos))
                            read_snp_sites.append(str(query_pos))
                elif ct==4 or ct==5:
                    query_pos+=cl
                else:
                    continue
        else:
            for (ct,cl) in read.cigartuples:
                if ct==0:
                    for i in range(cl):
                        query_pos+=1
                        ref_pos+=1
                        if ref_pos in snp_sites:
                            # read_snp_sites.append(str(ref_pos)+':'+str(read_length-query_pos-1))
                            read_snp_sites.append(str(read_length-query_pos-1))
                elif ct==1:
                    query_pos+=cl
                elif ct==2:
                    for i in range(cl):
                        ref_pos+=1
                        if ref_pos in snp_sites:
                            # read_snp_sites.append(str(ref_pos)+':'+str(read_length-query_pos-1))
                            read_snp_sites.append(str(read_length-query_pos-1))
                elif ct==4 or ct==5:
                    query_pos+=cl
                else:
                    continue
        if len(read_snp_sites)==0:
            continue
        print(rname+":"+str(ref_pos)+'\t'+','.join(read_snp_sites)+'\n')

def get_reference_length(fai_file):
    read_legnths = {}
    with open(fai_file) as f:
        for line in f:
            line = line.strip().split()
            read_legnths[line[0]] = int(line[1])
    return read_legnths

read_lengths = get_reference_length(fai)
# print(read_lengths)
annotate(aln,ctg,pos1,pos2,read_lengths)