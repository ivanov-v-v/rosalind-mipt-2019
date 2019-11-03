#!/home/vvi/anaconda3/bin/python

def str_from_kmers(kmers):
    k = len(kmers[0])
    kmers = sorted(kmers)
    output = kmers[0]
    while len(kmers) > 1:
        for i, kmer in enumerate(kmers):
            start_len = len(output)
            if output[-(k-1):] == kmer[:-1]:
                output += kmer[-1]
            elif kmer[1:] == output[:k-1]:
                output = kmer[0] + output
            if start_len != len(output): 
                kmers.pop(i)
                break
    return output

lines = [s.strip() for s in list(open("data.txt", "r"))]
k = int(lines[0])
kmers = lines[1:]

with open("answers.txt", "w") as outfile:
    outfile.write(str_from_kmers(kmers))
