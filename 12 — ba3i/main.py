#!/home/vvi/anaconda3/bin/python

def k_universal_circular_str(k):
    """
    Using '0 is better than 1' rule,
    constructs binary DeBrujin sequence
    containing all k-mers, and then converts
    it to circular one by stripping the tail of 1s.
    """
    seen_kmers = set()
    output = ''.join(['1' for _ in range(k)])
    seen_kmers.add(output)
    while True:
        tail = output[-(k-1):]
        if tail + '0' in seen_kmers:
            if tail + '1' in seen_kmers:
                break
            output = output + '1'
            seen_kmers.add(tail + '1')
        else:
            output = output + '0'
            seen_kmers.add(tail + '0')
    assert len(seen_kmers) == 2**k
    return output[:-(k-1)]

k = int(list(open("data.txt", "r"))[0])
with open("answers.txt", "w") as outfile:
    answer = k_universal_circular_str(k)
    outfile.write(answer)
