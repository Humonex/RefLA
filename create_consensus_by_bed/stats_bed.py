import numpy as np
import sys

def stats(bed_file):
    with open(bed_file, "r")as f:
        ls = []
        for line in f:
            if not line or line.startswith("#"): continue
            ctg, start, end = line.strip().split("\t")[:3]
            start, end = int(start), int(end)
            ls.append(end - start)
        
        ## stats
        print(sorted(ls))
        print("reg num: ", len(ls))
        print("mean: ", np.mean(ls))
        print("median: ", np.median(ls))
        step = 0.01
        quantile_ls = [i * step for i in range(int(1 / step))]
        quantile_num_ls = np.quantile(ls, quantile_ls)
        # if reverse_flag:
        #     quantile_num_ls = np.flip(quantile_num_ls)
        for i in range(len(quantile_ls)):
            print("{:.3f}:{:.0f}".format(quantile_ls[i], quantile_num_ls[i]))

if __name__ == "__main__":
    bed_file = sys.argv[1]
    stats(bed_file)
