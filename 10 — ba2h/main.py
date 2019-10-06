#!/home/vvi/anaconda3/bin/python

from collections import Counter
import itertools

import numpy as np


def get_kmers(dna, k):
    return [np.array(list(dna[i : i + k]))
            for i in range(len(dna) - k)]


def kmer_based_distance(pattern, reads):
    pattern = np.array(list(pattern))
    k = pattern.size
    return np.sum([
        np.min([
            (pattern != kmer).sum() 
            for kmer in get_kmers(dna, k)
        ]).astype(np.int32)
        for dna in reads
    ])
        


lines = list(open("data.txt", "r"))

with open("answers.txt", "w") as outfile:
    outfile.write(
        str(kmer_based_distance(
            lines[0].strip(), 
            lines[1].strip().split()
        ))
    )

