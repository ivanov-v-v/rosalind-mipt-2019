#!/home/vvi/anaconda3/bin/python

import numpy as np

lines = list(open("data.txt", "r"))
curr_line_id = 0

with open("answers.txt", "w") as outfile:
    while curr_line_id < len(lines):
        pattern = np.array(list(lines[curr_line_id].strip()))
        genome = np.array(list(lines[curr_line_id + 1].strip()))
        d = int(lines[curr_line_id + 2])
        curr_line_id += 3
        
#        print(pattern, genome, d)

        fuzzy_match_positions = []
        for i in range(genome.size - pattern.size + 1):
            mismatch_mask = ~(pattern == genome[i : i + pattern.size])
#            print(i, np.sum(mismatch_mask))
            if np.count_nonzero(mismatch_mask) <= d:
                fuzzy_match_positions.append(i)
        answers = ' '.join([str(i) for i in fuzzy_match_positions])
#        print(answers)
        outfile.write("{}\n".format(answers))

