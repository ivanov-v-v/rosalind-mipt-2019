#!/home/vvi/anaconda3/bin/python
from collections import defaultdict
import numpy as np

patterns = [s.strip() for s in list(open("data.txt", "r"))]
print(len(patterns) - np.unique(patterns).size)
k = len(patterns[0])
print(patterns)

graph = defaultdict(list)
indeg, outdeg = defaultdict(int), defaultdict(int)

nodes = []
for kmer in patterns:
    pref, suff = kmer[:-1], kmer[1:]
    outdeg[pref] += 1
    indeg[suff] += 1
    nodes.extend([pref, suff])
    graph[pref].append(suff)
is_ok = {v : indeg[v] == 1 and outdeg[v] == 1 for v in nodes}

#print(nodes)

contigs = set()
visited = defaultdict(int)
for u in nodes:
    if not is_ok[u]:
        visited[u] = True
        if outdeg[u] > 0:
            for w in graph[u]:
                path = [u, w[-1]]
                v = w
                while is_ok[v]:
                    visited[v] = True
                    v = graph[v][0]
                    path.append(v[-1])
                contigs.add(''.join(path))
#print(contigs)
for u in nodes:
    if visited[u] > 0:
        continue
    assert indeg[u] == 1 and outdeg[u] == 1
    visited[u] = True
    cycle = [u]
    v = u
    while True:
        v = graph[v][0]
        assert indeg[v] == 1 and outdeg[v] == 1, print(u, v)
        visited[v] = True
        if v == u:
            break
        cycle.append(v[-1])
    contigs.add(''.join(cycle))

answer = ' '.join(sorted(contigs))
print(answer)
with open("my_answers.txt", "w") as outfile:
    outfile.write(answer)
