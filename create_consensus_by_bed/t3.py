import sys
import fasta_parser
import re
'''

'''
"ctg1_s_e-ctg2_s_e"
# ctg1_1_1000-ctg2_1_1000
contig_sep = "~"
in_sep = "_"
fa1 = sys.argv[1]
fa2 = sys.argv[2]

# NC_000008.11:244300-6221362_NC_000008.11:8246062-11966644_NC_000008.11:12638544-43907689
fa_dic = fasta_parser.read_sequence_dict(fa1)
new_dic = {}
# str1 = "NC_000008.11:244300-6221362_NC_000008.11:8246062-11966644_NC_000008.11:12638544-43907689"
# str2 = "NC_000008.11:244300-6221362"
# str3 = "284_fragemant"
# str_ls = [str1, str2, str3]
for scaff_id, seq in fa_dic.items():
# for str in str_ls:
    try:
        if "NC" in scaff_id or "CM" in scaff_id:
            if "NC" in scaff_id:
                split_s = re.split("(NC)", scaff_id)[1:] # CM
            elif "CM" in scaff_id:
                split_s = re.split("(CM)", scaff_id)[1:] # CM
            split_s = [i+j for i,j in zip(split_s[::2], split_s[1::2])]
            split_s = [s.strip("_") for s in split_s]
            # print(split_s)
            new_id = ""
            for ctg_id in split_s:
                ctg, pos_ls = ctg_id.split(":")
                start, end = pos_ls.split("-")
                # print(ctg, start, end)
                ctg_id = in_sep.join([ctg, start, end])
                new_id += contig_sep + ctg_id
            new_id = new_id[1:]
            # if len(split_s) > 0:

        else:
            split_s = scaff_id # CM
            # print(split_s)
            new_id = split_s
        print(new_id)
        new_dic[new_id] = seq
    except:
        print(ValueError)
fasta_parser.write_fasta_dict(new_dic, fa2)
    
