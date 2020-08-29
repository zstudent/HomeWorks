# Question 3

# Say we are concerned only with overlaps that (a) are exact matches (no differences allowed), and (b) are at least
# k\verb|k|k bases long. To make an overlap graph, we could call overlap(a, b, min_length=k)\verb|
# overlap(a, b, min_length=k)|overlap(a, b, min_length=k) on every possible pair of reads from the dataset.
# Unfortunately, that will be very slow! Consider this: Say we are using k=6, and we have a read a\verb|a|a whose
# length-6 suffix is GTCCTA\verb|GTCCTA|GTCCTA. Say GTCCTA\verb|GTCCTA|GTCCTA does not occur in any other read
# in the dataset. In other words, the 6-mer GTCCTA\verb|GTCCTA|GTCCTA occurs at the end of read a\verb|a|a and
# nowhere else. It follows that a\verb|a|a's suffix cannot possibly overlap the prefix of any other read by 6 or
# more characters. Put another way, if we want to find the overlaps involving a suffix of read a\verb|a|a and a prefix
# of some other read, we can ignore any reads that don't contain the length-k suffix of a\verb|a|a. This is good news
# because it can save us a lot of work! Here is a suggestion for how to implement this idea. You don't have to do it
# this way, but this might help you. Let every k-mer in the dataset have an associated Python set\verb|set|set object,
# which starts out empty. We use a Python dictionary to associate each k-mer with its corresponding set\verb|set|set.
# (1) For every k-mer in a read, we add the read to the set\verb|set|set object corresponding to that k-mer.
# If our read is GATTA\verb|GATTA|GATTA and k=3, we would add GATTA\verb|GATTA|GATTA to the set\verb|set|set objects
# for GAT\verb|GAT|GAT, ATT\verb|ATT|ATT and TTA\verb|TTA|TTA. We do this for every read so that, at the end, each
# set\verb|set|set contains all reads containing the corresponding k-mer. (2) Now, for each read a\verb|a|a, we find
# all overlaps involving a suffix of a\verb|a|a. To do this, we take a\verb|a|a's length-k suffix,
# find all reads containing that k-mer (obtained from the corresponding set\verb|set|set) and call
# overlap(a, b, min_length=k)\verb|overlap(a, b, min_length=k)|overlap(a, b, min_length=k) for each.
# The most important point is that we do not call overlap(a, b, min_length=k)\verb|overlap(a, b, min_length=k)|
# overlap(a, b, min_length=k) if b\verb|b|b does not contain the length-k suffix of a\verb|a|a.
# Download and parse the read sequences from the provided Phi-X FASTQ file. We'll just use their base sequences,
# so you can ignore read names and base qualities. Also, no two reads in the FASTQ have the same sequence of bases.
# This makes things simpler. find all pairs of reads with an exact suffix / prefix match of length at least 30.
# Don't overlap a read with itself; if a read has a suffix/prefix match to itself, ignore that match.
# Ignore reverse complements.

from collections import defaultdict


def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip()  # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1  # move just past previous match


def overlap_graph(reads, k):
    # Make index
    index = defaultdict(set)
    for read in reads:
        for i in range(len(read) - k + 1):
            index[read[i:i + k]].add(read)

    # Make graph
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
