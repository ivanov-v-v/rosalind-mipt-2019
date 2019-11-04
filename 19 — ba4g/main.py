from collections import Counter, defaultdict
import numpy as np
import numba
from tqdm import tqdm

aminoacid_lst = np.array([
    57, 71, 87, 97, 99, 101, 103, 
    113, 114, 115, 128, 129, 131, 
    137, 147, 156, 163, 186 
])  


#@numba.jit(nopython=True)
def expand(leaderboard):
    expanded_leaderboard = []
    for peptide in leaderboard:
        for aminoacid in aminoacid_lst:
            expanded_leaderboard.append(peptide + [aminoacid])
    return expanded_leaderboard


#@numba.jit(nopython=True)
def mass(peptide):
    return sum(peptide)


#@numba.jit(nopython=True)
def parent_mass(spectrum):
    return spectrum[-1]


#@numba.jit(nopython=True)
def cyclospectrum(peptide):
    selfconcat = peptide + peptide
    spectrum = [0]
    for l in range(1, len(peptide)):
        for i in range(len(peptide)):
            spectrum.append(mass(selfconcat[i : i + l]))
    spectrum.append(mass(peptide))
    return sorted(spectrum)


def linspectrum(peptide):
    spectrum = [0] 
    for l in range(1, len(peptide)):
        for i in range(len(peptide) - l + 1):
            spectrum.append(mass(peptide[i : i + l]))
    spectrum.append(mass(peptide))
    return sorted(spectrum)


def score(peptide, spectrum):
    peptide_counts = Counter(cyclospectrum(peptide))
    theor_counts = Counter(spectrum)
    return np.sum([
        min(peptide_counts[base], theor_counts[base])
        for base in peptide_counts.keys()
    ])


def cut(leaderboard, spectrum, n):
    return sorted(
        leaderboard, 
        key=lambda peptide: score(peptide, spectrum), 
        reverse=True
    )[:n]

def peptide_repr(peptide):
    return '-'.join([str(x) for x in peptide])


#@numba.jit(nopython=True)
used = defaultdict(bool)
def leaderboard_cyclopeptide_sequencing(spectrum, n):
    leaderboard = [[]]
    leader_peptide = []
    while leaderboard:
#        print(len(leaderboard))
        leaderboard = expand(leaderboard)
        candidate_lst = []
        for peptide in tqdm(leaderboard):
            used[peptide_repr] = True
#            print(mass(peptide))
            if mass(peptide) == parent_mass(spectrum):
#                print(peptide, cyclospectrum(peptide))
                if score(peptide, spectrum) > score(leader_peptide, spectrum):
                    leader_peptide = peptide
                    print(leader_peptide)
            elif mass(peptide) < parent_mass(spectrum) and not used[peptide_repr(peptide)]:
                candidate_lst.append(peptide)
        leaderboard = cut(candidate_lst, spectrum, n)
        #print(leaderboard)
    return leader_peptide
    
lines = [s.strip() for s in list(open("data.txt", "r"))]
n = int(lines[0])
spectrum = [int(c) for c in lines[1].split()]

with  open("answers.txt", "w") as outfile:
    answer = '-'.join([str(mass) for mass in leaderboard_cyclopeptide_sequencing(spectrum, n)])
    print(answer)
    outfile.write(answer)
