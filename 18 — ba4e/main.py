#!/home/vvi/anaconda3/bin/python
import numpy as np
import numba
from tqdm import tqdm

aminoacid_lst = np.array([
    57, 71, 87, 97, 99, 101, 103, 
    113, 114, 115, 128, 129, 131, 
    137, 147, 156, 163, 186 
])  



#@numba.jit(nopython=True)
def expand(peptide_lst):
    expanded_peptide_lst = []
    for peptide in peptide_lst:
        for aminoacid in aminoacid_lst:
            expanded_peptide_lst.append(peptide + aminoacid)
    return expanded_peptide_lst


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
            spectrum.append(mass(selfconcat[i : i + l + 1]))
    spectrum.append(mass(peptide))
    return sorted(spectrum)


#@numba.jit(nopython=True)
def cyclopeptide_sequencing(spectrum):
    peptide_lst = [[]]
    result_lst = []
    curr_iter = 0
    while peptide_lst:
        print(len(peptide_lst))
        peptide_lst = expand(peptide_lst)
        peptide_to_try_lst = []
        for peptide in tqdm(peptide_lst, desc="iter {}".format(curr_iter + 1)):
            if mass(peptide) == parent_mass(spectrum):
                if cyclospectrum(peptide) == spectrum:
                    result_lst.append('-'.join(peptide))
                    print(result_lst[-1])
            elif np.all(np.isin(cyclospectrum(peptide), spectrum)):
                peptide_to_try_lst.append(peptide)
        peptide_lst = peptide_to_try_lst
        curr_iter += 1
    return result_lst
    
lines = [s.strip() for s in list(open("data.txt", "r"))]
spectrum = [int(c) for c in lines[0].split()]

with  open("answers.txt", "w") as outfile:
    answer = ' '.join(cyclopeptide_sequencing(spectrum))
    print(answer)
    outfile.write(answer)
