#!/home/vvi/anaconda3/bin/python

from collections import Counter
import itertools
from tqdm import tqdm

import numpy as np
import scipy.stats as sps


def get_kmers(dna, k):
    return [np.array(list(dna[i : i + k]))
            for i in range(len(dna) - k)]


def get_base_counts(dna):
    bases = ['A', 'T', 'G', 'C']
    return np.array([
        np.sum(dna == base) 
        for base in bases
    ])


def get_base_pseudocounts(dna):
    return get_base_counts(dna) + 1


base_to_num = {
    'A' : 0, 'T' : 1,
    'G' : 2, 'C' : 3
}

num_to_base = {
    val : key
    for key, val in base_to_num.items()
}


def get_profile(motifs):
    motif_mx = np.vstack(motifs)
    t, k = motif_mx.shape
    counts = np.column_stack([
        get_base_pseudocounts(col)
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
        (motif != consensus).sum()
        for motif in motif_mx
    )


def randomized_motif_search(reads, k, t):
    n = len(reads[0])
    motifs = [
        np.array(list(reads[i][start : start + k]))
        for i, start in enumerate(np.random.randint(low=0, high=n-k, size=t))
    ]
    best_motifs = motifs
    while True:
        profile = get_profile(motifs)
        motifs = [
            find_most_probable_motif(dna, k, profile)
            for dna in reads
        ]
        if get_score(motifs) < get_score(best_motifs):
            best_motifs = motifs
        else: 
            return best_motifs


def sample_from_profile(dna, profile):
    k = profile.shape[1]
    kmers_with_probas = np.vstack([
        (kmer, eval_kmer_proba(kmer, profile))
        for kmer in get_kmers(dna, k)
    ])
    kmers = kmers_with_probas[:, 0]
    probas = kmers_with_probas[:, 1]
    probas /= probas.sum()
    kmer_dist = sps.rv_discrete(values=(np.arange(len(kmers)), probas))
    return kmers[kmer_dist.rvs()]
    

def gibbs_sampler_motif_search(reads, k, t, n_restarts=100):
    n = len(reads[0])
    motifs = [
        np.array(list(reads[i][start : start + k]))
        for i, start in enumerate(np.random.randint(low=0, high=n-k, size=t))
    ]
    best_motifs = motifs
    for _ in range(n_restarts):
        i = np.random.randint(t)
        profile = get_profile(motifs[:i] + motifs[i+1:])
        missing_motif = sample_from_profile(reads[i], profile)
        motifs[i] = missing_motif
        if get_score(motifs) < get_score(best_motifs):
            best_motifs = motifs
    return best_motifs


lines = [line.strip() for line in list(open("data.txt", "r"))]
target_motifs = [np.array(list(line.strip())) for line in list(open("target.txt", "r"))]
k, t, n_restarts = (int(num) for num in lines[0].split())
reads = lines[1:]
MC_ITERS = 20
best_motifs = [ 
    np.array(list(dna[:k]))
    for dna in reads
]  
best_score = get_score(best_motifs)
outfile = open("answers.txt", "w")
for _ in tqdm(range(MC_ITERS), desc="gibbs sampler iterations passed"):
    candidate_motifs = gibbs_sampler_motif_search(reads, k, t, n_restarts)
    candidate_score = get_score(candidate_motifs)
    if candidate_score < best_score:
        best_score = candidate_score
        best_motifs = candidate_motifs
        print("CURRENT", best_score)
        print("TARGET: ", get_score(target_motifs))
        outfile.seek(0)        
        answers = '\n'.join([''.join(motif) for motif in best_motifs])
        outfile.write(answers)
outfile.close()
