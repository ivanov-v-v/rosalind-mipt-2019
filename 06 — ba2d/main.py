#!/home/vvi/anaconda3/bin/python

from collections import Counter
import itertools

import numpy as np
from scipy.spatial import distance


def get_kmers(dna, k):
    return [np.array(list(dna[i : i + k]))
            for i in range(len(dna) - k)]


def get_base_counts(dna):
    bases = ['A', 'T', 'G', 'C']
    return np.array([
        np.sum(dna == base) 
        for base in bases
    ])


base_to_num = {
    'A' : 0, 'T' : 1,
    'G' : 2, 'C' : 3
}

num_to_base = {
    val : key
    for key, val in base_to_num.items()
}


def get_profile(motif_mx):
    t, k = motif_mx.shape
    counts = np.column_stack([
        get_base_counts(col)
        for col in motif_mx.T
    ])
    assert counts.shape == (4, k)
    return counts / counts.sum(axis=0)


def eval_kmer_proba(kmer, profile):
    return np.prod([
        profile[base_to_num[kmer[i]], i]    
        for i in range(len(kmer))
    ])


def find_most_probable_motif(dna, k, profile):
    kmers_with_probas = np.vstack([
        (kmer, eval_kmer_proba(kmer, profile))
        for kmer in get_kmers(dna, k)
    ])
    top_row = np.argmax(kmers_with_probas[:, 1])
    return kmers_with_probas[top_row, 0]


def get_consensus(profile):
    return np.array([
        num_to_base[np.argmax(col)]
        for col in profile.T
    ])


def get_score(motifs):
    motif_mx = np.vstack(motifs)
    profile = get_profile(motif_mx)
    consensus = get_consensus(profile)
    return np.sum(
        distance.hamming(motif, consensus)
        for motif in motif_mx
    )


lines = list(open("data.txt", "r"))
k, t = (int(num) for num in lines[0].split())
best_motifs = [
    np.array(list(dna[:k]))
    for dna in lines[1:]
]

with open("answers.txt", "w") as outfile:
    ref_dna = lines[1]
    for ref_motif in get_kmers(ref_dna, k):
        motifs = [ref_motif]
        for dna in lines[2:]:
            dna = np.array(list(dna))
            motif_mx = np.vstack(motifs)
            profile = get_profile(motif_mx)
            likely_motif = find_most_probable_motif(dna, k, profile)
            motifs.append(likely_motif)
        if get_score(motifs) < get_score(best_motifs):
            best_motifs = motifs

    answers = '\n'.join([''.join(motif) for motif in best_motifs])
    outfile.write(answers)

