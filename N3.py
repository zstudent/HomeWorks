from collections import defaultdict
def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()
            seq = fh.readline().rstrip()
            fh.readline()
            qual = fh.readline().rstrip()
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


def overlap(a, b, min_length=3):
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0

        if b.startswith(a[start:]):
            return len(a) - start
        start += 1


def overlap_graph(reads, k):

    index = defaultdict(set)
    for read in reads:
        for i in range(len(read) - k + 1):
            index[read[i:i + k]].add(read)


    graph = defaultdict(set)
    for r in reads:
        for o in index[r[-k:]]:
            if r != o:
                if overlap(r, o, k):
                    graph[r].add(o)

    edges = 0
    for read in graph:
        edges += len(graph[read])
    return (edges, len(graph))


seqs, quals = readFastq('ERR266411_1.for_asm.fastq')
edges, suffixes = overlap_graph(seqs, 30)
print(edges)