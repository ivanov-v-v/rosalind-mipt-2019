#!/home/vvi/anaconda3/bin/python

from collections import Counter
import numpy as np
import itertools

def hamming_distance(arr1, arr2):
    mismatch_mask = ~(arr1 == arr2)
    return np.count_nonzero(mismatch_mask)


def find_fuzzy_matches(pattern, genome, d):
    fuzzy_match_positions = []
    for i in range(genome.size - pattern.size + 1): 
        if hamming_distance(pattern, genome[i : i + pattern.size]) <= d:
            fuzzy_match_positions.append(i)
    return fuzzy_match_positions


def generate_kmers(k):
    return itertools.product("ACGT", repeat=k)


def reverse_complement(pattern):
    base2complement = {"A" : "T", "C" : "G", "G" : "C", "T" : "A"}
    return np.array([base2complement[base] for base in pattern[::-1]])


lines = list(open("data.txt", "r"))
curr_line_id = 0

with open("answers.txt", "w") as outfile:
    while curr_line_id < len(lines):
        genome = np.array(list(lines[curr_line_id].strip()))
        (k, d) = [int(num) for num in lines[curr_line_id + 1].strip().split()]
        curr_line_id += 2

        kmer_freq = Counter()
        for kmer in generate_kmers(k):
            pattern = np.array(kmer)
            fuzzy_match_positions = find_fuzzy_matches(pattern, genome, d)
            revcomp_fuzzy_match_positions = find_fuzzy_matches(
                reverse_complement(pattern), genome, d
            )
            kmer_freq[''.join(kmer)] = (
                len(fuzzy_match_positions) 
                + len(revcomp_fuzzy_match_positions)
            )
        top_freq = kmer_freq.most_common(1)[0][1]
        top_kmers = []
        for kmer, freq in kmer_freq.most_common():
            if freq < top_freq:
                break
            top_kmers.append(kmer)
        answers = ' '.join(top_kmers)
        outfile.write("{}\n".format(answers))

