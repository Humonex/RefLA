
'''
'''
import gc
import logging
import cigar
import pysam

dic_starnd = {1: '+', 2: '-'}
signal = {1 << 2: 0, \
            1 >> 1: 1, \
            1 << 4: 2, \
            1 << 11: 3, \
            1 << 4 | 1 << 11: 4}
'''
    1 >> 1 means normal_foward read
    1 << 2 means unmapped read
    1 << 4 means reverse_complement read
    1 << 11 means supplementary alignment read
    1 << 4 | 1 << 11 means supplementary alignment with reverse_complement read
'''
def detect_flag(Flag):
    back_sig = signal[Flag] if Flag in signal else 0
    return back_sig

def analysis_inv(ele_1, ele_2, read_name, candidate, SV_size):
    ''' diff strand
    ele: 
    [304, 1721, 221967, 223351, 'NC_001133.9', '-']
    [read_start, read_end, ref_start, ref_end, chr, strand]

    Return candidate:
    ['++', 71637, 73397, 'S7_4563', 'INV', 'chrVII', 2], read1与ref同向 head-to-head
    ['--', 71640, 73407, 'S7_675', 'INV', 'chrVII', 3]， read1与ref反向 tail-to-tail
    '''
    '''
    abs(ele_1[3] - ele_2[3]): inv_size
    1、inv_size > min_size  是否再加一个max_size
    2、read2_start + 1/2(inv_size) >= read1_end; 小于的话说明两个比对交叠太多了/或者是说明这个inv太小，为假信号。过滤较少
    '''
    # INV = []
    # inv_size = abs(ele_1[3] - ele_2[3])
    if ele_1[5] == '+':
        # +-
        if ele_1[3] - ele_2[3] >= SV_size:
            if ele_2[0] + 0.5 * (ele_1[3] - ele_2[3]) >= ele_1[1]:  # 该条件用来过滤一部分
                candidate.append(["++", 
                                    ele_2[3], 
                                    ele_1[3], 
                                    read_name,
                                    "INV",
                                    ele_1[4]])
                # INV.append(["++", 
                #                     ele_2[3], 
                #                     ele_1[3], 
                #                     read_name,
                #                     "INV",
                #                     ele_1[4], 1])
                # head-to-head
                # 5'->5'
            # else:
            #     print("Failed INV:", ele_1, ele_2)
        if ele_2[3] - ele_1[3] >= SV_size:
            if ele_2[0] + 0.5 * (ele_2[3] - ele_1[3]) >= ele_1[1]:
                candidate.append(["++", 
                                    ele_1[3], 
                                    ele_2[3], 
                                    read_name,
                                    "INV",
                                    ele_1[4]])
                # INV.append(["++", 
                #                     ele_1[3], 
                #                     ele_2[3], 
                #                     read_name,
                #                     "INV",
                #                     ele_1[4], 2])
                # head-to-head
                # 5'->5'
            # else:
            #     print("Failed INV:", ele_1, ele_2)
    else:
        # -+
        if ele_2[2] - ele_1[2] >= SV_size:
            if ele_2[0] + 0.5 * (ele_2[2] - ele_1[2]) >= ele_1[1]:
                candidate.append(["--", 
                                    ele_1[2], 
                                    ele_2[2], 
                                    read_name,
                                    "INV",
                                    ele_1[4]])
                # INV.append(["--", 
                #                     ele_1[2], 
                #                     ele_2[2], 
                #                     read_name,
                #                     "INV",
                #                     ele_1[4], 3])
                # tail-to-tail
                # 3'->3'
            # else:
            #     print("Failed INV:", ele_1, ele_2)
        if ele_1[2] - ele_2[2] >= SV_size:
            if ele_2[0] + 0.5 * (ele_1[2] - ele_2[2]) >= ele_1[1]:
                candidate.append(["--", 
                                    ele_2[2], 
                                    ele_1[2], 
                                    read_name,
                                    "INV",
                                    ele_1[4]])
                # INV.append(["--", 
                #                     ele_2[2], 
                #                     ele_1[2], 
                #                     read_name,
                #                     "INV",
                #                     ele_1[4], 4])
                # tail-to-tail
                # 3'->3'
            # else:
            #     print("Failed INV:", ele_1, ele_2)
    # print("INV:", ele_1, ele_2, "->", INV, "size:", inv_size)

def analysis_bnd(ele_1, ele_2, read_name, candidate):
    '''
    *********Description*********
    *	TYPE A:		N[chr:pos[	*
    *	TYPE B:		N]chr:pos]	*
    *	TYPE C:		[chr:pos[N	*
    *	TYPE D:		]chr:pos]N	*
    *****************************
    '''
    if ele_2[0] - ele_1[1] <= 100:
        if ele_1[5] == '+':
            if ele_2[5] == '+':
                # +&+
                if ele_1[4] < ele_2[4]:
                    candidate.append(['A', 
                                        ele_1[3], 
                                        ele_2[4], 
                                        ele_2[2], 
                                        read_name,
                                        "TRA",
                                        ele_1[4]])
                    # N[chr:pos[
                else:
                    candidate.append(['D', 
                                        ele_2[2], 
                                        ele_1[4], 
                                        ele_1[3], 
                                        read_name,
                                        "TRA",
                                        ele_2[4]])
                    # ]chr:pos]N
            else:
                # +&-
                if ele_1[4] < ele_2[4]:
                    candidate.append(['B', 
                                        ele_1[3], 
                                        ele_2[4], 
                                        ele_2[3], 
                                        read_name,
                                        "TRA",
                                        ele_1[4]])
                    # N]chr:pos]
                else:
                    candidate.append(['B', 
                                        ele_2[3], 
                                        ele_1[4], 
                                        ele_1[3], 
                                        read_name,
                                        "TRA",
                                        ele_2[4]])
                    # N]chr:pos]
        else:
            if ele_2[5] == '+':
                # -&+
                if ele_1[4] < ele_2[4]:
                    candidate.append(['C', 
                                        ele_1[2], 
                                        ele_2[4], 
                                        ele_2[2], 
                                        read_name,
                                        "TRA",
                                        ele_1[4]])
                    # [chr:pos[N
                else:
                    candidate.append(['C', 
                                        ele_2[2], 
                                        ele_1[4], 
                                        ele_1[2], 
                                        read_name,
                                        "TRA",
                                        ele_2[4]])
                    # [chr:pos[N
            else:
                # -&-
                if ele_1[4] < ele_2[4]:
                    candidate.append(['D', 
                                        ele_1[2], 
                                        ele_2[4], 
                                        ele_2[3], 
                                        read_name,
                                        "TRA",
                                        ele_1[4]])
                    # ]chr:pos]N
                else:
                    candidate.append(['A', 
                                        ele_2[3], 
                                        ele_1[4], 
                                        ele_1[2], 
                                        read_name,
                                        "TRA",
                                        ele_2[4]])
                    # N[chr:pos[

def analysis_split_read(split_read, SV_size, RLength, read_name, candidate, MaxSize, query):
    '''split_read
    read_start	read_end	ref_start	ref_end	chr	strand
    #0			#1			#2			#3		#4	#5
    '''
    
    '''
    split_read: [[304, 1721, 221967, 223351, 'NC_001133.9', '-'], [1797, 2815, 1524211, 1525231, 'NC_001136.10', '-']]
    SV_size:min_size
    RLength:query_length
    query:read.query_sequence
    '''
    '''
    1、对于一类读数：比对到同一条染色体但相距比较远的位置，应该视为什么signal。有关relocations???
    '''
    SP_list = sorted(split_read, key = lambda x:x[0])   # read based sort
    # print(SP_list)
    if len(SP_list) < 2:    # 可改进之处？？，0/1的话，说明primary比对太差了??
        # print(SP_list)
        return

    # detect INS involoved in a translocation   ????
    trigger_INS_TRA = 0	

    '''sig = False'''
    # Store Strands of INV
    if len(SP_list) == 2:   # [[304, 1721, 221967, 223351, 'NC_001133.9', '-'], [1797, 2815, 1524211, 1525231, 'NC_001136.10', '-']]
        ele_1 = SP_list[0]
        ele_2 = SP_list[1]
        if ele_1[4] == ele_2[4]:    # same chr:dup、ins、del、inv
            '''if ele_2[2] - ele_1[3] > 10000:
                sig = True
                print("Relocation:", ele_1, ele_2, read_name)
            else:
                sig = False'''
            
            if ele_1[5] != ele_2[5]:    # diff strand: inv
                analysis_inv(ele_1, 
                                ele_2, 
                                read_name, 
                                candidate,
                                SV_size)

            else:   # same strand: dup & ins & del ...
                # dup & ins & del
                a = 0
                if ele_1[5] == '-': # ele_1 ele_2颠换顺序
                    # print(ele_1, "->", [RLength-SP_list[a+1][1], RLength-SP_list[a+1][0]]+SP_list[a+1][2:], ele_2, "->", [RLength-SP_list[a][1], RLength-SP_list[a][0]]+SP_list[a][2:])
                    ele_1 = [RLength-SP_list[a+1][1], RLength-SP_list[a+1][0]]+SP_list[a+1][2:]
                    ele_2 = [RLength-SP_list[a][1], RLength-SP_list[a][0]]+SP_list[a][2:]
                    query = query[::-1] # 逆序query_sequence
                
                # print(ele_1, ele_2)
                if ele_1[3] - ele_2[2] >= SV_size:  # 第一段的末尾在第二段的起始的后面。典型的DUP比对
                    '''
                    INS和DUP都是序列插入，但是DUP由于是复制前面一段序列，所以第二段比对与前一段有交叠
                    条件设置
                    1、read gap size大于ref上的gap size
                    '''
                    # if ele_2[1] - ele_1[1] >= ele_1[3] - ele_2[2]:
                    if ele_2[0] - ele_1[1] >= ele_1[3] - ele_2[2]:  # read gap size大于ref上的gap size
                        candidate.append([int((ele_1[3]+ele_2[2])/2), 
                                        ele_2[0]+ele_1[3]-ele_2[2]-ele_1[1], 
                                        read_name,
                                        str(query[ele_1[1]+int((ele_1[3]-ele_2[2])/2):ele_2[0]-int((ele_1[3]-ele_2[2])/2)]),
                                        "INS",
                                        ele_2[4]])
                    else:
                        candidate.append([ele_2[2], 
                                            ele_1[3], 
                                            read_name,
                                            "DUP",
                                            ele_2[4]])
                    # if ele_2[2] < ele_1[2]:
                    #     if len(candidate[-1]) == 6:
                    #         print(ele_1, ele_2, "->candidate:", candidate[-1][:3] + candidate[-1][4:])
                    #     else:
                    #         print(ele_1, ele_2, "->candidate:", candidate[-1])

                '''
                read_gap: ele_2[0] - ele_1[1]
                ref_gap: ele_2[2] - ele_1[3]
                1、read gap > ref gap: INS
                2、read gap < ref gap: DEL
                '''
                delta_length = ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1]    # (ele_2[0] - ele_1[1]) - (ele_2[2] - ele_1[3]) : read_gap - ref_gap
                if ele_1[3] - ele_2[2] < max(SV_size, delta_length/5) and delta_length >= SV_size:  # 过滤掉第一段末尾在第二段起始之后的读数、过滤掉del
                    if ele_2[2] - ele_1[3] <= max(100, delta_length/5) and (delta_length <= MaxSize or MaxSize == -1):
                        candidate.append([(ele_2[2]+ele_1[3])/2, 
                                            delta_length, 
                                            read_name,
                                            str(query[ele_1[1]+int((ele_2[2]-ele_1[3])/2):ele_2[0]-int((ele_2[2]-ele_1[3])/2)]),
                                            "INS",
                                            ele_2[4]])
                delta_length = ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3]    # (ele_2[2] - ele_1[3]) - (ele_2[0] - ele_1[1]) : ref_gap - read_gap
                if ele_1[3] - ele_2[2] < max(SV_size, delta_length/5) and delta_length >= SV_size:
                    if ele_2[0] - ele_1[1] <= max(100, delta_length/5) and (delta_length <= MaxSize or MaxSize == -1):
                        candidate.append([ele_1[3], 
                                            delta_length, 
                                            read_name,
                                            "DEL",
                                            ele_2[4]])
        else:   # diff chr -> tra
            trigger_INS_TRA = 1
            analysis_bnd(ele_1, ele_2, read_name, candidate)

    else:
        
        # over three splits
        for a in range(len(SP_list[1:-1])):     # = len(SP_list)-2
            ele_1 = SP_list[a]
            ele_2 = SP_list[a+1]
            ele_3 = SP_list[a+2]

            if ele_1[4] == ele_2[4]:    # same chr:dup、ins、del、inv
                if ele_2[4] == ele_3[4]:
                    if ele_1[5] == ele_3[5] and ele_1[5] != ele_2[5]:
                        if ele_2[5] == '-':
                            # +-+
                            if ele_2[0] + 0.5 * (ele_3[2] - ele_1[3]) >= ele_1[1] and ele_3[0] + 0.5 * (ele_3[2] - ele_1[3]) >= ele_2[1]:
                                # No overlaps in split reads

                                if ele_2[2] >= ele_1[3] and ele_3[2] >= ele_2[3]:
                                    candidate.append(["++", 
                                                        ele_1[3], 
                                                        ele_2[3], 
                                                        read_name,
                                                        "INV",
                                                        ele_1[4]])
                                    # head-to-head
                                    # 5'->5'
                                    candidate.append(["--", 
                                                        ele_2[2], 
                                                        ele_3[2], 
                                                        read_name,
                                                        "INV",
                                                        ele_1[4]])
                                    # tail-to-tail
                                    # 3'->3'
                        else:
                            # -+-
                            if ele_1[1] <= ele_2[0] + 0.5 * (ele_1[2] - ele_3[3]) and ele_3[0] + 0.5 * (ele_1[2] - ele_3[3]) >= ele_2[1]:
                                # No overlaps in split reads

                                if ele_2[2] - ele_3[3] >= -50 and ele_1[2] - ele_2[3] >= -50:
                                    candidate.append(["++", 
                                                        ele_3[3], 
                                                        ele_2[3], 
                                                        read_name,
                                                        "INV",
                                                        ele_1[4]])
                                    # head-to-head
                                    # 5'->5'
                                    candidate.append(["--", 
                                                        ele_2[2], 
                                                        ele_1[2], 
                                                        read_name,
                                                        "INV",
                                                        ele_1[4]])
                                    # tail-to-tail
                                    # 3'->3'	

                    if len(SP_list) - 3 == a:
                        if ele_1[5] != ele_3[5]:
                            if ele_2[5] == ele_1[5]:
                                # ++-/--+
                                analysis_inv(ele_2, 
                                                ele_3, 
                                                read_name, 
                                                candidate, 
                                                SV_size)
                            else:
                                # +--/-++
                                analysis_inv(ele_1, 
                                                ele_2, 
                                                read_name, 
                                                candidate, 
                                                SV_size)

                    if ele_1[5] == ele_3[5] and ele_1[5] == ele_2[5]:
                        # dup & ins & del 
                        if ele_1[5] == '-':
                            ele_1 = [RLength-SP_list[a+2][1], RLength-SP_list[a+2][0]]+SP_list[a+2][2:]
                            ele_2 = [RLength-SP_list[a+1][1], RLength-SP_list[a+1][0]]+SP_list[a+1][2:]
                            ele_3 = [RLength-SP_list[a][1], RLength-SP_list[a][0]]+SP_list[a][2:]
                            query = query[::-1]

                        if ele_2[3] - ele_3[2] >= SV_size and ele_2[2] < ele_3[3]:
                            candidate.append([ele_3[2], 
                                                ele_2[3], 
                                                read_name,
                                                "DUP",
                                                ele_2[4]])

                        if a == 0:
                            if ele_1[3] - ele_2[2] >= SV_size:
                                candidate.append([ele_2[2], 
                                                    ele_1[3], 
                                                    read_name,
                                                    "DUP",
                                                    ele_2[4]])

                        delta_length = ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1]
                        if ele_1[3] - ele_2[2] < max(SV_size, delta_length/5) and delta_length >= SV_size:
                            if ele_2[2] - ele_1[3] <= max(100, delta_length/5) and (delta_length <= MaxSize or MaxSize == -1):
                                if ele_3[2] >= ele_2[3]:
                                    candidate.append([(ele_2[2]+ele_1[3])/2, 
                                                        delta_length, 
                                                        read_name,
                                                        str(query[ele_1[1]+int((ele_2[2]-ele_1[3])/2):ele_2[0]-int((ele_2[2]-ele_1[3])/2)]),
                                                        "INS",
                                                        ele_2[4]])
                        delta_length = ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3]
                        if ele_1[3] - ele_2[2] < max(SV_size, delta_length/5) and delta_length >= SV_size:
                            if ele_2[0] - ele_1[1] <= max(100, delta_length/5) and (delta_length <= MaxSize or MaxSize == -1):
                                if ele_3[2] >= ele_2[3]:
                                    candidate.append([ele_1[3], 
                                                        delta_length, 
                                                        read_name,
                                                        "DEL",
                                                        ele_2[4]])
                        
                        if len(SP_list) - 3 == a:
                            ele_1 = ele_2
                            ele_2 = ele_3

                            delta_length = ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1]
                            if ele_1[3] - ele_2[2] < max(SV_size, delta_length/5) and delta_length >= SV_size:
                                if ele_2[2] - ele_1[3] <= max(100, delta_length/5) and (delta_length <= MaxSize or MaxSize == -1):
                                    candidate.append([(ele_2[2]+ele_1[3])/2, 
                                                        delta_length, 
                                                        read_name,
                                                        str(query[ele_1[1]+int((ele_2[2]-ele_1[3])/2):ele_2[0]-int((ele_2[2]-ele_1[3])/2)]),
                                                        "INS",
                                                        ele_2[4]])

                            delta_length = ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3]
                            if ele_1[3] - ele_2[2] < max(SV_size, delta_length/5) and ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3] >= SV_size:
                                if ele_2[0] - ele_1[1] <= max(100, delta_length/5) and (delta_length <= MaxSize or MaxSize == -1):
                                    candidate.append([ele_1[3], 
                                                        delta_length, 
                                                        read_name,
                                                        "DEL",
                                                        ele_2[4]])

                    if len(SP_list) - 3 == a and ele_1[5] != ele_2[5] and ele_2[5] == ele_3[5]:
                        ele_1 = ele_2
                        ele_2 = ele_3
                        ele_3 = None
                    if ele_3 == None or (ele_1[5] == ele_2[5] and ele_2[5] != ele_3[5]):
                        if ele_1[5] == '-':
                            ele_1 = [RLength-SP_list[a+2][1], RLength-SP_list[a+2][0]]+SP_list[a+2][2:]
                            ele_2 = [RLength-SP_list[a+1][1], RLength-SP_list[a+1][0]]+SP_list[a+1][2:]
                            query = query[::-1]
                        delta_length = ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1]
                        if ele_1[3] - ele_2[2] < max(SV_size, delta_length/5) and delta_length >= SV_size:
                            if ele_2[2] - ele_1[3] <= max(100, delta_length/5) and (delta_length <= MaxSize or MaxSize == -1):
                                candidate.append([(ele_2[2]+ele_1[3])/2, 
                                                    delta_length, 
                                                    read_name,
                                                    str(query[ele_1[1]+int((ele_2[2]-ele_1[3])/2):ele_2[0]-int((ele_2[2]-ele_1[3])/2)]),
                                                    "INS",
                                                    ele_2[4]])

                        delta_length = ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3]
                        if ele_1[3] - ele_2[2] < max(SV_size, delta_length/5) and delta_length >= SV_size:
                            if ele_2[0] - ele_1[1] <= max(100, delta_length/5) and (delta_length <= MaxSize or MaxSize == -1):
                                candidate.append([ele_1[3], 
                                                    delta_length, 
                                                    read_name,
                                                    "DEL",
                                                    ele_2[4]])

            else:
                trigger_INS_TRA = 1
                analysis_bnd(ele_1, ele_2, read_name, candidate)

                if len(SP_list) - 3 == a:
                    if ele_2[4] != ele_3[4]:
                        analysis_bnd(ele_2, ele_3, read_name, candidate)

    if len(SP_list) >= 3 and trigger_INS_TRA == 1:
        if SP_list[0][4] == SP_list[-1][4]:
            # print(SP_list[0])
            # print(SP_list[-1])
            if SP_list[0][5] != SP_list[-1][5]:
                pass
            else:
                if SP_list[0][5] == '+':
                    ele_1 = SP_list[0]
                    ele_2 = SP_list[-1]
                else:
                    ele_1 = [RLength-SP_list[-1][1], RLength-SP_list[-1][0]]+SP_list[-1][2:]
                    ele_2 = [RLength-SP_list[0][1],RLength-SP_list[0][0]]+SP_list[0][2:]
                    query = query[::-1]
                # print(ele_1)
                # print(ele_2)
                dis_ref = ele_2[2] - ele_1[3]
                dis_read = ele_2[0] - ele_1[1]
                if dis_ref < 100 and dis_read - dis_ref >= SV_size and (dis_read - dis_ref <= MaxSize or MaxSize == -1):
                    # print(min(ele_2[2], ele_1[3]), dis_read - dis_ref, read_name)
                    candidate.append([min(ele_2[2], ele_1[3]), 
                                        dis_read - dis_ref, 
                                        read_name,
                                        str(query[ele_1[1]+int(dis_ref/2):ele_2[0]-int(dis_ref/2)]),
                                        "INS",
                                        ele_2[4]])	

                if dis_ref <= -SV_size:
                    candidate.append([ele_2[2], 
                                        ele_1[3], 
                                        read_name,
                                        "DUP",
                                        ele_2[4]])
    '''if sig:
        print(candidate[-1])'''
def acquire_clip_pos(deal_cigar):
    '''
    acquire the pos of clip in the query sequence
    Return: [left_clip, right_clip, bias]
    '''
    seq = list(cigar.Cigar(deal_cigar).items()) # [(6636, 'S'), (1515, 'M'), (165, 'I'), (9, 'S')]
    
    if seq[0][1] == 'S':
        first_pos = seq[0][0]
    else:
        first_pos = 0
    if seq[-1][1] == 'S':
        last_pos = seq[-1][0]
    else:
        last_pos = 0

    bias = 0
    for i in seq:
        if i[1] == 'M' or i[1] == 'D' or i[1] == '=' or i[1] == 'X':
            bias += i[0]
    # print(seq, "->", [first_pos, last_pos, bias]) # [(1141, 'S'), (3804, 'M'), (595, 'I'), (8569, 'S')] [1141, 8569, 3804]
    return [first_pos, last_pos, bias]

def organize_split_signal(primary_info, Supplementary_info, total_L, SV_size, 
    min_mapq, max_split_parts, read_name, candidate, MaxSize, query):
    '''
    +:primary_info:[softclip_left, read.query_length-softclip_right, pos_start, pos_end, Chr_name, strand]
    -:primary_info:[softclip_right, read.query_length-softclip_left, pos_start, pos_end, Chr_name, strand]
    Supplementary_info: from SA tag
    ['NC_001136.10,,-,15237541512M4I6002S,25,276', 'NC_001140.6,551510,-,6562S883M70I3S,12,251']
    [ref, ref_start, strand, cigar, mapq, NM]   NM:
    total_L:query_length
    query:query_sequence
    '''
    # print(primary_info, Supplementary_info, total_L, SV_size, min_mapq, max_split_parts, read_name, candidate, MaxSize)
    # if len(primary_info) == 0:
    #     print(Supplementary_info)
    # print(primary_info, Supplementary_info)
    # print("NM:", Supplementary_info[0].split(',')[-1])
    split_read = list()
    if len(primary_info) > 0:   # if read.mapq < min_mapq -> primary_info = []
        split_read.append(primary_info)
        min_mapq = 0

    for i in Supplementary_info:
        seq = i.split(',')  # [NC_001136.10, 1523754, -, 1512M4I6002S, 25, 276]: [ref, ref_start, strand, cigar, mapq, NM]
        local_chr = seq[0]
        local_start = int(seq[1])
        local_cigar = seq[3]
        local_strand = seq[2]
        local_mapq = int(seq[4])
        if local_mapq >= min_mapq:
        # if local_mapq >= 0:	
            local_set = acquire_clip_pos(local_cigar)   # local_set: [left_clip, right_clip, bias]
            if local_strand == '+':
                 split_read.append([local_set[0], total_L-local_set[1], local_start, 
                     local_start+local_set[2], local_chr, local_strand])    #  local_set[0]:相当于是该条比对在queryseq中的比对起始位置，total_L-local_set[1]是比对结束位置
            else:   #  - strand要颠倒位置，因为我们存的是比对起始位置，要将读数作为正向为标准。
                try:
                    split_read.append([local_set[1], total_L-local_set[0], local_start, 
                        local_start+local_set[2], local_chr, local_strand])
                except:
                    # print("Error")
                    pass
    if len(split_read) <= max_split_parts or max_split_parts == -1:
        # print("split_read:", split_read)
        # if len(split_read) < 2:
        #     print(split_read)
        analysis_split_read(split_read, SV_size, total_L, read_name, candidate, MaxSize, query)

def generate_combine_sigs(sigs, Chr_name, read_name, svtype, candidate, merge_dis):
    '''
    sigs:
    1、ins: [sig_start, oplen, str(read.query_sequence[shift_ins_read-oplen:shift_ins_read])]   sig_start:在ref上的sigs对应左端位置
    2、del: [sig_start, oplen].
    
    candidate: None。元素如下
    1、ins: [sig_start, sv_len, read_name, temp_sig[2], svtype, Chr_name]
    2、del: [sig_start, sv_len, read_name, svtype, Chr_name]
    '''
    # for i in sigs:
    # 	print(svtype,i, len(sigs))
    # if len(sigs) > 0:
    #     print(sigs, Chr_name, read_name, svtype, candidate, merge_dis)
    
    if len(sigs) == 0:
        pass
    elif len(sigs) == 1:
        if svtype == 'INS':
            candidate.append([sigs[0][0], 
                                            sigs[0][1], 
                                            read_name,
                                            sigs[0][2],
                                            svtype,
                                            Chr_name])
        else:   # svtype='DEL'
            candidate.append([sigs[0][0], 
                                            sigs[0][1], 
                                            read_name,
                                            svtype,
                                            Chr_name])
    else:   # 同一个读数中的多个sigs，merge
        temp_sig = sigs[0]
        # print(sigs[0], temp_sig)
        if svtype == "INS":
            temp_sig += [sigs[0][0]]
            for i in sigs[1:]:
                if i[0] - temp_sig[3] <= merge_dis:
                    temp_sig[1] += i[1]
                    temp_sig[2] += i[2]
                    temp_sig[3] = i[0]  # 新的sigs start，主要作用是用于Merge
                    # print(temp_sig)
                else:
                    # if temp_sig[0] !=temp_sig[3]:
                    #     print("merge:", temp_sig[0], "->", temp_sig[3])
                    candidate.append([temp_sig[0], 
                                                        temp_sig[1], 
                                                        read_name,
                                                        temp_sig[2],
                                                        svtype,
                                                        Chr_name])
                    temp_sig = i
                    temp_sig.append(i[0])
            candidate.append([temp_sig[0], 
                                                temp_sig[1], 
                                                read_name,
                                                temp_sig[2],
                                                svtype,
                                                Chr_name])
        else:   # del
            # print(sigs) # [[223127, 13], [224020, 10], [226913, 72]]
            temp_sig += [sum(sigs[0])]
            # merge_dis_bias = max([i[1]] for i in sigs)
            for i in sigs[1:]:
                if i[0] - temp_sig[2] <= merge_dis:
                    temp_sig[1] += i[1]
                    temp_sig[2] = sum(i)
                else: 
                    candidate.append([temp_sig[0], 
                                                        temp_sig[1], 
                                                        read_name,
                                                        svtype,
                                                        Chr_name])
                    temp_sig = i
                    temp_sig.append(i[0])
            candidate.append([temp_sig[0], 
                                                temp_sig[1], 
                                                read_name,
                                                svtype,
                                                Chr_name])
OPLIST=[
    pysam.CBACK,
    pysam.CDEL,
    pysam.CDIFF,
    pysam.CEQUAL,
    pysam.CHARD_CLIP,
    pysam.CINS,
    pysam.CMATCH,
    pysam.CPAD,
    pysam.CREF_SKIP,
    pysam.CSOFT_CLIP
]
RefChangeOp=set([0,2,7,8])

#QUERY CHANGE, REF CHANGE
CHANGETABLE={
    pysam.CMATCH:     (True,True),
    pysam.CINS:       (True,False),
    pysam.CDEL:       (False,True),
    pysam.CREF_SKIP:  (False,True),
    pysam.CPAD:       (False,False),
    pysam.CEQUAL:     (True,True),
    pysam.CDIFF:      (True,True)
}

CHANGEOP=[CHANGETABLE[i] if i in CHANGETABLE.keys() else (False,False) for i in range(max(OPLIST)+1)]
REFCHANGEOP=[CHANGETABLE[i][1] if i in CHANGETABLE.keys() else False for i in range(max(OPLIST)+1)]
INDELOP=[(i==pysam.CDEL or i==pysam.CINS) for i in range(max(OPLIST)+1)]

def parse_read(read, Chr_name, SV_size, min_mapq, max_split_parts, min_read_len, min_siglength, merge_del_threshold, merge_ins_threshold, MaxSize):
    # read = pysam.AlignedSegment()

    if read.query_length < min_read_len:
        return []
    all_candidate = list()
    cigar_candidate = list()
    seg_candidate = list()
    candidate = list()
    Combine_sig_in_same_read_ins = list()
    Combine_sig_in_same_read_del = list()

    process_signal = detect_flag(read.flag)
    if read.mapq >= min_mapq:
        pos_start = read.reference_start # 0-based
        pos_end = read.reference_end
        # shift_del = 0
        # shift_ins = 0
        sig_start = pos_start
        softclip_left = 0
        softclip_right = 0
        hardclip_left = 0
        hardclip_right = 0
        shift_ins_read = 0
        if read.cigar[0][0] == 4:
            softclip_left = read.cigar[0][1]
        if read.cigar[0][0] == 5:
            hardclip_left = read.cigar[0][1]

        '''for element in read.cigar:
            if element[0] in [0, 7 ,8]:
                shift_del += element[1]
            if element[0] == 2 and element[1] < min_siglength: ## changed SV_size to min_siglength
                shift_del += element[1]
            if element[0] == 2 and element[1] >= min_siglength: ## changed SV_size to min_siglength
                Combine_sig_in_same_read_del.append([pos_start+shift_del, element[1]])
                shift_del += element[1]

            # calculate offset of an ins sig in read
            if element[0] != 2:
                shift_ins_read += element[1]

            if element[0] in [0, 2, 7, 8]:
                shift_ins += element[1]
            if element[0] == 1 and element[1] >= min_siglength: ## changed SV_size to min_siglength
                Combine_sig_in_same_read_ins.append([pos_start+shift_ins, element[1],
                    str(read.query_sequence[shift_ins_read-element[1]-hardclip_left:shift_ins_read-hardclip_left])])'''

        for op, oplen in read.cigartuples:
            # calculate offset of an ins sig in read
            if op != 2:#might be fixed later
                shift_ins_read += oplen
            if oplen >= min_siglength and INDELOP[op]:
                if op==2:
                    Combine_sig_in_same_read_del.append([sig_start, oplen])
                    sig_start += oplen
                else:
                    Combine_sig_in_same_read_ins.append([sig_start, oplen,
                        str(read.query_sequence[shift_ins_read-oplen:shift_ins_read])])
            else:
                # if op in RefChangeOp:
                if REFCHANGEOP[op]:
                    sig_start += oplen
        
        if read.cigar[-1][0] == 4:
            softclip_right = read.cigar[-1][1]
        if read.cigar[-1][0] == 5:
            hardclip_right = read.cigar[-1][1]

        if hardclip_left != 0:
            softclip_left = hardclip_left
        if hardclip_right != 0:
            softclip_right = hardclip_right

    # ************Combine signals in same read********************
    generate_combine_sigs(Combine_sig_in_same_read_ins, Chr_name, read.query_name, "INS", cigar_candidate, merge_ins_threshold)
    generate_combine_sigs(Combine_sig_in_same_read_del, Chr_name, read.query_name, "DEL", cigar_candidate, merge_del_threshold)

    if process_signal == 1 or process_signal == 2: # 0 / 16     # primary
        Tags = read.get_tags()
        if read.mapq >= min_mapq:
            if process_signal == 1:
                primary_info = [softclip_left, read.query_length-softclip_right, pos_start, 
                pos_end, Chr_name, dic_starnd[process_signal]]
            else:
                primary_info = [softclip_right, read.query_length-softclip_left, pos_start, 
                pos_end, Chr_name, dic_starnd[process_signal]]
            # print("primary_info:", primary_info, softclip_left, softclip_right)
        else:
            primary_info = []

        for i in Tags:
            if i[0] == 'SA':    # ('SA', 'NC_001136.10,1523754,-,1512M4I6002S,25,276;NC_001140.6,551510,-,6562S883M70I3S,12,251;')
                # print("SA:", i)
                Supplementary_info = i[1].split(';')[:-1]   # 
                # print("Supplementary_info:", Supplementary_info)    # ['NC_001136.10,1523754,-,1512M4I6002S,25,276', 'NC_001140.6,551510,-,6562S883M70I3S,12,251']
                organize_split_signal(primary_info, Supplementary_info, read.query_length, 
                    SV_size, min_mapq, max_split_parts, read.query_name, seg_candidate, MaxSize, read.query_sequence)
    # if len(candidate) > 1:
    #     print("candidate:", candidate)  # candidate: [[222345, 13, 'S12_2760', 'ACCGACTGCAATC', 'INS', 'NC_001133.9'], [223127, 11, 'S12_2760', 'DEL', 'NC_001133.9']]
    for ele in cigar_candidate:
        ele[-2] = ele[-2] + "-cigar"
        all_candidate.append(ele)
    for ele in seg_candidate:
        ele[-2] = ele[-2] + "-seg"
        all_candidate.append(ele)
    return all_candidate
    # return candidate

def single_pipe(sam_path, min_length, min_mapq, max_split_parts, min_read_len, temp_dir, 
                task, min_siglength, merge_del_threshold, merge_ins_threshold, MaxSize, bed_regions):
    print("Start for:{}:{}-{}".format(task[0], task[1], task[2]))
    candidate = list()
    reads_info_list = list()
    Chr_name = task[0]
    samfile = pysam.AlignmentFile(sam_path)

    for read in samfile.fetch(Chr_name, task[1], task[2]):
        # print("flag:", read.flag, read.is_secondary, read.is_supplementary)
        if read.flag == 256 or read.flag == 272:    # filter secondary reads
            continue
        pos_start = read.reference_start # 0-based
        pos_end = read.reference_end
        in_bed = False
        if bed_regions != None:
            for bed_region in bed_regions:
                if pos_end <= bed_region[0] or pos_start >= bed_region[1]:
                    continue
                else:
                    in_bed = True
                    break
        else:
            in_bed = True
        if read.reference_start >= task[1] and in_bed:
            read_candidate = parse_read(read, Chr_name, min_length, min_mapq, max_split_parts, 
                                    min_read_len, min_siglength, merge_del_threshold, 
                                    merge_ins_threshold, MaxSize)
            candidate.extend(read_candidate)
            if read.mapq >= min_mapq:
                is_primary = 0
                if read.flag in [0, 16]:
                    is_primary = 1
                reads_info_list.append([pos_start, pos_end, is_primary, read.query_name])
    samfile.close()
    # print('finish %s:%d-%d in %f seconds.'%(task[0], task[1], len(reads_info_list), time.time() - start_time))
   
    # if len(candidate) == 0:
    #     logging.info("Skip %s:%d-%d."%(Chr_name, task[1], task[2]))
    #     return

    output = "%ssignatures/_%s_%d_%d.bed"%(temp_dir, Chr_name, task[1], task[2])
    file = open(output, 'w')
    for ele in candidate:
        if len(ele) == 5:
            file.write("%s\t%s\t%d\t%d\t%s\n"%(ele[-2], ele[-1], ele[0], ele[1], ele[2]))
        elif len(ele) == 7:
            file.write("%s\t%s\t%s\t%d\t%s\t%d\t%s\n"%(ele[-2], ele[-1], ele[0], 
                ele[1], ele[2], ele[3], ele[4]))
        elif len(ele) == 6:
            try:
                file.write("%s\t%s\t%s\t%d\t%d\t%s\n"%(ele[-2], ele[-1], ele[0], ele[1], ele[2], ele[3]))
                # INV chr strand pos1 pos2 read_ID
            except:
                file.write("%s\t%s\t%d\t%d\t%s\t%s\n"%(ele[-2], ele[-1], ele[0], ele[1], ele[2], ele[3]))
                # INS chr pos len read_ID seq
    file.close()
    reads_output = "%ssignatures/_%s_%d_%d.reads"%(temp_dir, Chr_name, task[1], task[2])
    reads_file = open(reads_output, 'w')
    for ele in reads_info_list: # [pos_start, pos_end, is_primary, read.query_name]
        reads_file.write("%s\t%d\t%d\t%d\t%s\n"%(Chr_name, ele[0], ele[1], ele[2], ele[3]))
    reads_file.close()
    logging.info("Finished %s:%d-%d."%(Chr_name, task[1], task[2]))	
    gc.collect()

if __name__ == '__main__':
    # run(sys.argv[1:])
    '''
    args.input, 
    args.min_size, 
    args.min_mapq, 
    args.max_split_parts, 
    args.min_read_len, 
    temporary_dir, 
    Task_list[i],   # [i[0], pos, local_ref_len]
    args.min_siglength, 
    args.merge_del_threshold, 
    args.merge_ins_threshold, 
    args.max_size,
    '''
    print("Test")
    ## Test1
    # sam_path = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/aln2alt.sort.bam"
    # min_length = 30
    # min_mapq = 20
    # max_split_parts = 7
    # min_read_len = 500
    # temp_dir = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/cute_SV2/"
    # task = ["NC_001133.9", 0, 230218]
    # min_siglength = 10
    # merge_del_threshold = 0     # 500   # sigs的合并距离
    # merge_ins_threshold = 100   # 500
    # MaxSize = 100000
    # bed_regions = None
    # Test2
    sam_path = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast/simu/aln2trim.sort.bam"
    min_length = 30
    min_mapq = 20
    max_split_parts = 7
    min_read_len = 500
    temp_dir = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast/cute_SV/"
    task = ["chrVII", 0, 1123738]
    min_siglength = 10
    merge_del_threshold = 0     # 500   # sigs的合并距离
    merge_ins_threshold = 100   # 500
    MaxSize = 100000
    bed_regions = None
    single_pipe(sam_path, min_length, min_mapq, max_split_parts, min_read_len, temp_dir, 
                task, min_siglength, merge_del_threshold, merge_ins_threshold, MaxSize, bed_regions)
    # Test3
    # bam_reader = pysam.AlignmentFile(sam_path)
    # task_ls = []
    # for chr in bam_reader.references:
    #     task_ls.append([chr, 0, bam_reader.get_reference_length(chr)])
    # for task in task_ls:
    #     single_pipe(sam_path, min_length, min_mapq, max_split_parts, min_read_len, temp_dir, 
    #             task, min_siglength, merge_del_threshold, merge_ins_threshold, MaxSize, bed_regions)
    pass