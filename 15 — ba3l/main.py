#!/home/vvi/anaconda3/bin/python
import numpy as np

def str_from_kmers(kmers):
    return kmers[0] + ''.join([kmer[-1] for kmer in kmers[1:]])

def reconstruct_from_paired_reads(paired_reads, k, d):
    pref = str_from_kmers(paired_reads[:, 0])
    suff = str_from_kmers(paired_reads[:, 1])
    for i in range(k + d + 1, len(pref)):
        if pref[i] != suff[i - k - d]:
            return ''
    return pref + suff[-(k + d):]

lines = [s.strip() for s in list(open("data.txt", "r"))]
k, d = [int(val) for val in lines[0].split()]
paired_reads = np.vstack([pair.split('|') for pair in lines[1:]])

reconstruction = reconstruct_from_paired_reads(paired_reads, k, d)
with  open("answers.txt", "w") as outfile:
    answer = ''.join(reconstruction)
    outfile.write(answer)
