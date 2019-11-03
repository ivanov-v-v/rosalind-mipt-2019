import numpy as np

def deBruijn_collection(pattern,\
                        head=lambda kmer: kmer[0:-1],
                        tail=lambda kmer: kmer[1:]):
    graph={}
    k=len(pattern[0])
    for kmer in pattern:
        if not head(kmer) in graph:
            graph[head(kmer)]=[]
        graph[head(kmer)].append(tail(kmer))
    for kmer in graph.keys():
        graph[kmer].sort()
    return graph

def non_branching_paths(graph):
    def add_counts(graph):
        counts={}
        for node,outs in graph.items():
            counts[node]=(0,len(outs))
        for node,outs in graph.items():
            for out in outs:
                if not out in counts:
                    counts[out]=(0,0)
                x,y=counts[out]
                counts[out]=x+1,y
        return counts
    def make_cycle(start,graph):
        def standardize(cycle):
            mm=min(cycle)
            ii=cycle.index(mm)
            return cycle[ii:]+cycle[:ii]+[mm]
        cycle=[start]
        node=start

        while node in graph and len(graph[node])==1:
            succ=graph[node][0]
            if succ==start:
                return standardize(cycle)
            else:
                cycle.append(succ)
                node=succ

        return []
    def isolated_cycles(nodes,graph):
        cycles=[]
        for node in graph.keys():
            cycle=make_cycle(node,graph)
            if len(cycle)>0 and not cycle in cycles:
                cycles.append(cycle)
        return cycles

    paths=[]
    nodes=add_counts(graph)

    for v in nodes:
        (ins,outs)=nodes[v]
        if ins!=1 or outs!=1:
            if outs>0:
                for w in graph[v]:
                    nbp=[v,w]
                    w_in,w_out=nodes[w]
                    while w_in==1 and w_out==1:
                        u=graph[w][0]
                        nbp.append(u)
                        w=u
                        w_in,w_out=nodes[w]
                    paths.append(nbp)
    return isolated_cycles(nodes,graph)+paths

def create_contigs(patterns):  
    contigs=[]
    for path in non_branching_paths(deBruijn_collection(patterns)):
        contig=path[0]
        for p in path[1:]:
            contig=contig+p[-1]
        contigs.append(contig)   
    return contigs

patterns = [s.strip() for s in list(open("data.txt", "r"))]
k = len(patterns[0])

contigs = sorted(create_contigs(patterns))
print(contigs)
answer = ' '.join(contigs)
with open("answers.txt", "w") as outfile:
    outfile.write(answer)
