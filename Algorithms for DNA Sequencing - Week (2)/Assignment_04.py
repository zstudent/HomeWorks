# Question 4
# Index-assisted approximate matching. In practicals, we built a Python class called Index\verb|Index|Index
# implementing an ordered-list version of the k-mer index. The Index\verb|Index|Index class is copied below.
# We also implemented the pigeonhole principle using Boyer-Moore as our exact matching algorithm.
# Implement the pigeonhole principle using Index\verb|Index|Index to find exact matches for the partitions.
# Assume P always has length 24, and that we are looking for approximate matches with up to 2 mismatches (substitutions).
# We will use an 8-mer index.
# Download the Python module for building a k-mer index.
# https://d28rh4a8wq0iu5.cloudfront.net/ads1/code/kmer_index.py
# Write a function that, given a length-24 pattern P and given an Index\verb|Index|Index object built on 8-mers,
# finds all approximate occurrences of P within T with up to 2 mismatches. Insertions and deletions are not allowed.
# Don't consider any reverse complements.
# How many times does the string GGCGCGGTGGCTCACGCCTGTAAT\verb|GGCGCGGTGGCTCACGCCTGTAAT|GGCGCGGTGGCTCACGCCTGTAAT,
# which is derived from a human Alu sequence, occur with up to 2 substitutions in the excerpt of human chromosome 1?
# (Don't consider reverse complements here.)
# Hint 1: Multiple index hits might direct you to the same match multiple times, but be careful not to count a match more than once.
# Hint 2: You can check your work by comparing the output of your new function to that of the naive_2mm\verb|naive_2mm|naive_2mm
# function implemented in the previous module.

import bisect

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

class Index(object):
    def __init__(self, t, k):
        ''' Create index from all substrings of size 'length' '''
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i + k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer

    def query(self, p):
        ''' Return index hits for first k-mer of P '''
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits


def approximate_match(p, t, n):
    segment_length = int(round(len(p) / (n+1)))
    all_matches = set()
    p_idx = Index(t, segment_length)
    idx_hits = 0
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        matches = p_idx.query(p[start:end])

        # Extend matching segments to see if whole p matches
        for m in matches:
            idx_hits += 1
            if m < start or m-start+len(p) > len(t):
                continue

            mismatches = 0

            for j in range(0, start):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            for j in range(end, len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break

            if mismatches <= n:
                all_matches.add(m - start)
    return list(all_matches), idx_hits


p = "GGCGCGGTGGCTCACGCCTGTAAT"
chr1 = readGenome('chr1.GRCh38.excerpt.fasta')

print(len(approximate_match(p, chr1, 2)[0]))
