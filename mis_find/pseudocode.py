
# def win_check(ctg, start, end, ctg_len, bam_reader, min_span_num, min_supp_portion, MIN_MAPPING_QUALITY, MIN_ALIGN_LENGTH, MIN_ALIGN_RATE, ins_threshold, del_threshold, min_clip_len):

#     span_ls = []
#     for read in bam_reader.fetch(ctg, start, end):
#         '''if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < MIN_MAPPING_QUALITY or read.query_alignment_length < MIN_ALIGN_LENGTH or (read.query_alignment_length / read.query_length) < MIN_ALIGN_RATE:
#             continue
#         cigar = read.cigartuples
#         left = cigar[0]
#         if left[0] == 4 or left[0] == 5:
#             if left[1] > min_clip_len:
#                 continue
#         right = cigar[-1]
#         if right[0] == 4 or right[0] == 5:
#             if right[1] > min_clip_len:
#                 continue
#         if read.reference_start <= max(0, start - 500) and read.reference_end >= min(ctg_len, end + 500):
#             span_ls.append(read)'''
#         if check_read(read):
#             span_ls.append(read)
#     if len(span_ls) < min_span_num and start >= 5000 and end <= ctg_len - 5000:
#         return True, "low_supp_mis"
#     else:
#         min_supp_num = max(min_span_num, min_supp_portion*len(span_ls))
#         # cal indel of the span read
#         Inss_ls = []
#         Dels_ls = []
#         for read in span_ls:   # 收集span_ls读数的indels
#             # read_l = max(start - 5000, read.reference_start + 100)
#             # read_r = min(end + 5000, read.reference_end - 100)
#             indels = cal_idels(read, read_l, read_r)
#             Inss_ls.append(indels[0])
#             Dels_ls.append(indels[1])
#         # 
#         Inss_ls.sort()
#         Dels_ls.sort()
#         # check indels
        
#         # avg_ins = sum(Inss_ls[:min_supp_num]) // min_supp_num
#         # avg_del = sum(Dels_ls[:min_supp_num]) // min_supp_num
#         avg_ins = get_avgins(Inss_ls, min_supp_num)
#         avg_del = get_avgdel(Dels_ls, min_supp_num)
#         if avg_ins < ins_threshold and avg_del < del_threshold: # no mis
#             return False, "No_mis"
#         else:
#             return True, "Reads_mis"

'''
输入：候选区间位置ctg, start, end
输入：BAM比对信息bam_reader，读数质量检查参数参数MIN_MAPPING_QUALITY, MIN_ALIGN_LENGTH, MIN_ALIGN_RATE, min_clip_len
输入：支持读数参数min_span_num, min_supp_portion, ins_threshold, del_threshold
输出：该窗口是否为读数支持窗口
span_ls = list()     # 初始化span_ls列表
FOR read IN bam_reader.fetch(ctg, start, end)
    IF check_read(read) THEN  # 调用读数质量检查函数
       APPEND read TO span_ls   # 将读数加入span_ls列表中
    ENDIF
ENDFOR

IF LEN(span_ls) < min_span_num AND THEN
    RETURN True, "low_supp_mis"  # 若跨越区域较少，读数支持不足
ELSE
    SET min_supp_num TO MAX(min_span_num, min_supp_portion * LEN(span_ls)) # 更新最小跨越数目
    Inss_ls = list() # 插入列表初始化
    Dels_ls = list() # 删除列表初始化
    FOR read IN span_ls
        SET indels TO cal_idels(read, read_l, read_r)  # 从span_ls读取中计算indels
        APPEND indels[0] TO Inss_ls # 插入列表追加indels中的插入事件
        APPEND indels[1] TO Dels_ls # 删除列表追加indels中的删除事件
    ENDFOR
    CALL SORT(Inss_ls) # 将插入列表排序
    CALL SORT(Dels_ls) # 将删除列表排序

    SET avg_ins TO Cal_avg_indel(Inss_ls, min_supp_num)  # 计算平均插入大小
    SET avg_del TO Cal_avg_indel(Dels_ls, min_supp_num)  # 计算平均删除大小

    IF avg_ins < ins_threshold AND avg_del < del_threshold THEN  # 检查indels是否在阈值之内
        RETURN False, "No_mis"  # 若indels小于阈值，该窗口为读数支持窗口，无组装错误
    ELSE
        RETURN True, "Reads_mis" # 否则，该窗口为读数不支持窗口，存在组装错误
    ENDIF
ENDIF


'''