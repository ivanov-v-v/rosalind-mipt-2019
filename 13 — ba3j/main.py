#!/home/vvi/anaconda3/bin/python
from collections import defaultdict
import numpy as np

lines = [s.strip() for s in list(open("data.txt", "r"))]
k, d = [int(val) for val in lines[0].split()]
paired_reads = [s.split('|') for s in lines[1:]]

graph = defaultdict(list)
nodes = []
source = None
balance = defaultdict(int)
for l, r in paired_reads:
    pref = (l[:-1], r[:-1])
    suff = (l[1:], r[1:])
    nodes.extend([pref, suff])
    graph[pref].append(suff)
    balance[pref] -= 1
    balance[suff] += 1

balance = np.array([balance[v] for v in nodes])
source = nodes[np.where(balance < 0)[0][0]]

curr_path, ecycle = [source], []
curr_node = source
while curr_path:
    if len(graph[curr_node]):
        curr_path.append(curr_node)
        curr_node = graph[curr_node].pop()
    else:
        ecycle.append(curr_node)
        curr_node = curr_path.pop()
ecycle.reverse()
ecycle = np.vstack(ecycle)

prefs = ecycle[0, 0] + ''.join([s[k-2:] for s in ecycle[1:, 0]])
suffs = ecycle[0, 1] + ''.join([s[k-2:] for s in ecycle[1:, 1]])
answer = prefs[:d + k] + suffs
with open("answers.txt", "w") as outfile:
    outfile.write(answer)
