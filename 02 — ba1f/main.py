#!/home/vvi/anaconda3/bin/python
import numpy as np

encoder = {"A" : 0, "G" : 1, "C" : -1, "T" : 0}

with open("data.txt", "r") as infile,\
        open("answers.txt", "w") as outfile:

    prefix_sums = np.cumsum([0] + [encoder[c] for c in infile.readline().strip()])  
    min_indices = np.ravel(np.where(prefix_sums == np.min(prefix_sums)))
    answer = ' '.join([str(idx) for idx in min_indices])
    outfile.write("{}\n".format(answer))



