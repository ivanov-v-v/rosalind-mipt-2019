#!/home/vvi/anaconda3/bin/python
import numpy as np
import numba

@numba.jit(nopython=True)
def n_proteins_by_mass(mass):
    aminoacid_masses = np.array([
        57, 71, 87, 97, 99, 101, 103, 
        113, 114, 115, 128, 129, 131, 
        137, 147, 156, 163, 186
    ])
    dp = np.zeros(mass + 1, dtype=np.int64)
    dp[0] = 1
    for i in range(mass + 1):
        for am_mass in aminoacid_masses:
            if am_mass > i: 
                break 
            dp[i] += dp[i - am_mass]
    return dp[mass]

lines = [s.strip() for s in list(open("data.txt", "r"))]
mass = int(lines[0])

with  open("answers.txt", "w") as outfile:
    answer = str(n_proteins_by_mass(mass))
    outfile.write(answer)
