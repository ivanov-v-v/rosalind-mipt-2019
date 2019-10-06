#!/home/vvi/anaconda3/bin/python

from collections import defaultdict
import numpy as np

lines = list(open("data.txt", "r"))
curr_line = 0

with open("answers.txt", "w") as outfile:
    while curr_line < len(lines):
        genome = lines[curr_line]
        (k, L, t) = [int(s) for s in lines[curr_line + 1].split()]
        curr_line += 2
#        print(genome, k, L, t)
        kmer2pos = defaultdict(list)
        for i in range(len(genome) - k + 1):
            kmer2pos[genome[i : i + k]].append(i)

#        print(kmer2pos)
        lt_clumps = []
        for kmer, positions in kmer2pos.items():
            for i in range(len(positions) - t + 1):
                if positions[i + t - 1] - positions[i] + 1 <= L:
                    lt_clumps.append(kmer)
                    break
#        print(lt_clumps)
        outfile.write("{}\n".format(' '.join(lt_clumps)))





