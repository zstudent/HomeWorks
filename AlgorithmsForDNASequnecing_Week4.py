#!/usr/bin/env python
# coding: utf-8

# In[31]:


import itertools
def scss(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    shortest_cnt = 1
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
            shortest_cnt = 1
        elif len(sup) == len(shortest_sup):
            shortest_cnt += 1
    return shortest_sup, shortest_cnt  # return shortest  


# In[30]:


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
            return len(a)-start
        start += 1  # move just past previous match


# In[19]:


def scs_with_counts(ss):
    """ Returns shortest common superstring of given strings,
        assuming no string is a strict substring of another """
    all_scs = []
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]
        for i in range(len(ss)-1):
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) <= len(shortest_sup):
            shortest_sup = sup
            all_scs.append(sup)
    return shortest_sup,all_scs
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


# In[20]:


def overlap_graph(reads, k):
    # Make index
    index = defaultdict(set)
    for read in reads:
        for i in range(len(read) - k + 1):
            index[read[i:i+k]].add(read)

    # Make graph
    graph = defaultdict(set)
    for r in reads:
        for o in index[r[-k:]]:
            if r != o:
                if overlap(r, o, k):
                    graph[r].add(o)

    return graph
def pick_maximal_overlap(reads, k):
    """ Return a pair of reads from the list with a
        maximal suffix/prefix overlap >= k.  Returns
        overlap length 0 if there are no such overlaps."""
    reada, readb = None, None
    best_olen = 0

    # Make index
    index = defaultdict(set)
    for read in reads:
        for i in range(len(read) - k + 1):
            index[read[i:i+k]].add(read)

    for r in reads:
        for o in index[r[-k:]]:
            if r != o:
                olen = overlap(r, o, k)
                if olen > best_olen:
                    reada, readb = r, o
                    best_olen = olen

    return reada, readb, best_olen


# In[21]:


def greedy_scs(reads, k):
    """ Greedy shortest-common-superstring merge.
        Repeat until no edges (overlaps of length >= k)
        remain. """
    read_a, read_b, olen = pick_maximal_overlap(reads, k)
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a, read_b, olen = pick_maximal_overlap(reads, k)
    return ''.join(reads)


# # 1)

# In[32]:


shortest, cnt = scss(['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT'])
print(len(shortest))


# # 2)

# In[35]:


shortest, scs = scs_with_counts(['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT'])
print(len(scs))


# # 3)

# In[28]:


sequences, qualities = readFastq('ads1_week4_reads.fq')
from collections import defaultdict
for k in range(100, 1, -1):
        genome = greedy_scs(sequences, k)
        if len(genome) == 15894:
            print(genome.count('A'))
            break


# # 4)

# In[29]:


print(genome.count('T'))

