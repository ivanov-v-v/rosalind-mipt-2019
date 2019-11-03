#!/home/vvi/anaconda3/bin/python
from collections import defaultdict
import numpy as np

def dna_to_rna(dna_str):
    return dna_str.replace("T", "U")


def rna_to_dna(rna_str):
    return rna_str.replace("U", "T")


def dna_revcomp(dna_str):
    base2complement = {
        "A" : "T", "T" : "A",
        "C" : "G", "G" : "C"
    }
    return ''.join([base2complement[base] for base in dna_str])[::-1]


def rna_revcomp(rna_str):
    return dna_to_rna(dna_revcomp(rna_to_dna(rna_str)))


codon2protein = {
    #AXX
    'AUG' : 'M', 'AUA' : 'I', 'AUC' : 'I',
    'AUU' : 'I', 'AGG' : 'R', 'AGA' : 'R',
    'AGC' : 'S', 'AGU' : 'S', 'ACG' : 'T',
    'ACA' : 'T', 'ACC' : 'T', 'ACU' : 'T',
    'AAG' : 'K', 'AAA' : 'K', 'AAC' : 'N', 'AAU' : 'N',
    #CXX
    'CAU' : 'H', 'CAC' : 'H', 
    'CAA' : 'Q', 'CAG' : 'Q', 
    'CCU' : 'P', 'CCC' : 'P', 'CCA' : 'P', 'CCG' : 'P',
    'CGU' : 'R', 'CGC' : 'R', 'CGA' : 'R', 'CGG' : 'R',
    'CUU' : 'L', 'CUC' : 'L', 'CUA' : 'L', 'CUG' : 'L',
    #GXX
    'GAU' : 'D', 'GAC' : 'D',
    'GAA' : 'E', 'GAG' : 'E',
    'GCU' : 'A', 'GCC' : 'A', 'GCA' : 'A', 'GCG' : 'A',
    'GGU' : 'G', 'GGC' : 'G', 'GGA' : 'G', 'GGG' : 'G',
    'GUU' : 'V', 'GUC' : 'V', 'GUA' : 'V', 'GUG' : 'V',
    #UXX
    'UAU' : 'Y', 'UAC' : 'Y',
    'UAA' : '*', 'UAG' : '*',
    'UCU' : 'S', 'UCC' : 'S', 'UCA' : 'S', 'UCG' : 'S',
    'UGU' : 'C', 'UGC' : 'C',
    'UGA' : '*',
    'UGG' : 'W',
    'UUU' : 'F', 'UUC' : 'F',
    'UUA' : 'L', 'UUG' : 'L'
}

protein2codon = defaultdict(list)
for key, val in codon2protein.items():
    protein2codon[val].append(key) 

lines = [s.strip() for s in list(open("data.txt", "r"))]
dna_str = lines[0]
protein_lst = list(lines[1])

n_codons = len(protein_lst)
enc_lst = []
for i in range(len(dna_str) - 3 * n_codons + 1):
    precodon_lst = [
        dna_str[i + 3 * j : i + 3 * (j + 1)] 
        for j in range(n_codons)
    ]
    revcomp_precodon_lst = [
        dna_revcomp(precodon) 
        for precodon in precodon_lst
    ][::-1]
    enc_at_i_lst = []
    for j in range(len(protein_lst)):
        precodon, revcomp_precodon = precodon_lst[j], revcomp_precodon_lst[j]
        protein = protein_lst[j]
        codon, revcomp_codon = dna_to_rna(precodon), dna_to_rna(revcomp_precodon)
        if (codon2protein[codon] != protein) and (codon2protein[revcomp_codon] != protein):
            break
        enc_at_i_lst.append(precodon)
    if len(enc_at_i_lst) == len(protein_lst):
        enc_lst.append(''.join(enc_at_i_lst))
    
with  open("answers.txt", "w") as outfile:
    answer = '\n'.join(enc_lst)
    outfile.write(answer)