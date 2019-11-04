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
def expand(leaderboard, allowed_aminoacids):
    expanded_leaderboard = []
    for peptide in leaderboard:
        for aminoacid in allowed_aminoacids:
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
    scores = [score(peptide, spectrum) for peptide in leaderboard]
    candidates = sorted(
        leaderboard, 
        key=lambda peptide: score(peptide, spectrum), 
        reverse=True
    )[:n]
    for i in range(len(leaderboard)):
        if scores[i] == score(candidates[-1], spectrum):
            if leaderboard[i] not in candidates:
                candidates.append(leaderboard[i])
    return candidates


def peptide_repr(peptide):
    return '-'.join([str(x) for x in peptide])


def convolve(spectrum):
    return sorted(
        Counter([
            abs(spectrum[i] - spectrum[j]) 
            for i in range(len(spectrum) - 1) 
            for j in range(i + 1, len(spectrum))
        ]).items(), key=lambda p: p[1], reverse=True
    )


#@numba.jit(nopython=True)
def leaderboard_cyclopeptide_sequencing(spectrum, n, allowed_aminoacids):
    leaderboard = [[]]
    leader_peptide = []
    used = defaultdict(bool)
    while leaderboard:
        leaderboard = expand(leaderboard, allowed_aminoacids)
        candidate_lst = []
        for peptide in leaderboard:
            used[peptide_repr] = True
            if mass(peptide) == parent_mass(spectrum):
                if score(peptide, spectrum) > score(leader_peptide, spectrum):
                    leader_peptide = peptide
                    print(peptide_repr(leader_peptide))
            elif mass(peptide) < parent_mass(spectrum) and not used[peptide_repr(peptide)]:
                candidate_lst.append(peptide)
        leaderboard = cut(candidate_lst, spectrum, n)
    return leader_peptide


def convolution_cyclopeptide_sequencing(spectrum, n, m):
    spectrum_conv = convolve(spectrum)
    i = 0
    allowed_aminoacids = set()
    while True:
        mass = spectrum_conv[i][0]
        if 57 <= mass  <= 200:
            allowed_aminoacids.add(mass)
        if len(allowed_aminoacids) == m:
            allowed_aminoacids.update([
                counter_item[0]
                for counter_item in spectrum_conv[i + 1:]
                if (counter_item[1] == spectrum_conv[i][1] 
                    and 57 <= counter_item[0] <= 200)
            ])
            break
        i += 1
    return leaderboard_cyclopeptide_sequencing(
        spectrum, n, allowed_aminoacids
    )

    
lines = [s.strip() for s in list(open("data.txt", "r"))]
m = int(lines[0])
n = int(lines[1])
spectrum = sorted([int(c) for c in lines[2].split()])

with  open("answers.txt", "w") as outfile:
    leader_peptide = convolution_cyclopeptide_sequencing(spectrum, n, m)
    print(
        peptide_repr(leader_peptide), 
        score(leader_peptide, spectrum)
    )
    outfile.write(peptide_repr(leader_peptide))
