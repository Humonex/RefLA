import pysam
import re
from create_consensus_by_bed.Utils import DepthRec
def convert_reference_pos_to_raw_pos2(read, candidate_pos):  
    # 输入read信息和要求的参考上的位置，将参考的位置转换为read_querysequence上的位置
    '''注意:这个位置并不一定是读数原来的位置'''
    candidate_pos = set(candidate_pos)
    raw_ref_pos_map={}
    ref_pos = read.reference_start
    query_pos = 0
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
    # return raw_ref_pos_map,read_length
    return raw_ref_pos_map
def solve_by_reads(bam_in, ctg, start, end, Depthrec:DepthRec, config):
    '''
    1、引入identity选择最优读数。计算方式：%identity based on NM 
    2、引入对区域depth的判断，若存在过高的depth，则False'''
    print("*****{}:{}-{}*******".format(ctg, start, end))
    # dp_upper_bound = config["dp_upper_bound"]
    block_high_dp_bound = config["block_high_dp_bound"]
    if Depthrec.get_block_high_dp([ctg, start, end]) / Depthrec.whole_dp > block_high_dp_bound:
        return (False, "", "")
    MIN_SUPPORT_READS = 3
    MIN_MAPPING_QUALITY = 40
    MIN_ALIGN_LENGTH = 10000
    MIN_ALIGN_RATE = 0.95
    MIN_CLIP_LEN = 500  # 对ont不知是否要开大点
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    # support_reads_ls = []
    support_read_sim_ls = []
    for read in bam_reader.fetch(ctg, start, end):
        ## filter1
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < MIN_MAPPING_QUALITY:
            continue    # 只保留primary
        ## filter2
        if read.query_alignment_length < MIN_ALIGN_LENGTH or (read.query_alignment_length / read.query_length) < MIN_ALIGN_RATE:
            # print(read.query_name, read.query_alignment_length, read.query_alignment_length / read.query_length)
            continue
        ## filter3
        if read.reference_start > start-1000 or read.reference_end < end + 1000:continue    # 跨度不够长，保留足够长的。暂不考虑端粒区域，默认组装不了
        cigar = read.cigarstring
        tokens = re.findall("[\d]{0,}[A-Z]{1}", cigar)
        left, right = tokens[0], tokens[-1]
        if left[-1] in "HS" and int(left[:-1]) > MIN_CLIP_LEN: continue     # 保留没有大段的剪切的读数
        if right[-1] in "HS" and int(right[:-1]) > MIN_CLIP_LEN: continue
        ##
        # support_reads_ls.append(read)
        support_read_sim_ls.append([read, 1-read.get_tag("NM")/read.query_length])
        # if read.reference_start < start-1000 and read.reference_end > end + 1000:
        #     support_reads_ls.append(read)
        # print("{}:{}-{}, {}".format(read.query_name, read.reference_start, read.reference_end, 1-read.get_tag("NM")/read.query_length))
    support_read_sim_ls.sort(key=lambda x:x[1], reverse=True)
    # for read,sim in support_read_sim_ls:print("{}:{}-{}, {}".format(read.query_name, read.reference_start, read.reference_end, sim))
    # for read,sim in support_read_sim_ls:print("{}".format(sim), end="\t")
    # print("\n")
    if len(support_read_sim_ls) > MIN_SUPPORT_READS:
        print("{}:{}-{}, solve with read: {}, similarity:{}".format(ctg, start, end, support_read_sim_ls[0][0].query_name, support_read_sim_ls[0][1]))
        best_read = support_read_sim_ls[0][0]
        ref_to_read = convert_reference_pos_to_raw_pos2(best_read, [start, end])
        query_start = ref_to_read.get(start, -1)
        query_end = ref_to_read.get(end, -1)
        if query_start < 0 or query_end < 0 or query_start > query_end:
            print("{}:{}-{} ERROR!!!".format(ctg, start, end))
            return (False, "", "")
        solve_seq = best_read.query_sequence[query_start:query_end]
        return (True, solve_seq, best_read.query_name)
    return (False, "", "")
def main():
    bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi_ctgs/NC_19/my_pipe/step1_mapping/aln.sorted.bam"
    
    # reg = "NC_000019.10"
    # reg = "NC_000019.10:27369501-27371000"
    # reg = "NC_000019.10:24558501-24564500"
    reg = "NC_000019.10:564501-564600"
    ctg = reg.split(":")[0]
    start, end = reg.split(":")[1].split("-")
    start, end = int(start), int(end)
    print(reg)
    res = solve_by_reads(bam_in, ctg, start, end)
    print(res)
    pass
if __name__ == "__main__":
    main()