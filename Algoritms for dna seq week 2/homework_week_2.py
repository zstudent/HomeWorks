# Question 1
# How many alignments does the naive exact matching algorithm try when matching the string
# GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG\verb|GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG|
# GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to the excerpt of human chromosome 1?
# (Don't consider reverse complements.)

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def naive(p, t):
    occurrences = []
    comparisons = 0
    alignments = 0
    for i in range(len(t) - len(p) + 1):
        alignments += 1
        match = True
        for j in range(len(p)):
            comparisons += 1
            if t[i+j] != p[j]:
                match = False
                break
        if match:
            occurrences.append(i)
    return occurrences, comparisons, alignments


p = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
chr1 = readGenome("chr1.GRCh38.excerpt.fasta")
print("Question 1")
print(naive(p, chr1)[2])

# Question 2
# How many character comparisons does the naive exact matching algorithm try when matching the string
# GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG\verb|GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG|
# GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to the excerpt of human chromosome 1?
# (Don't consider reverse complements.)




p = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
chr1 = readGenome("chr1.GRCh38.excerpt.fasta")

print("Question 2")
print(naive(p, chr1)[1])


# Question 3
# How many alignments does Boyer-Moore try when matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG\verb|
# GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG|GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG
# (derived from human Alu sequences) to the excerpt of human chromosome 1? (Don't consider reverse complements.)

from bm_preproc import BoyerMoore



def boyer_moore(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    comparisons = 0
    alignments = 0
    while i < len(t) - len(p) + 1:
        alignments += 1
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            comparisons += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences, comparisons, alignments


p = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
chr1 = readGenome("chr1.GRCh38.excerpt.fasta")
p_bm = BoyerMoore(p)
print("Question 3")
print(boyer_moore(p, p_bm, chr1)[2])


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



class Index(object):
    def __init__(self, t, k):
        self.k = k
        self.index = []
        for i in range(len(t) - k + 1):
            self.index.append((t[i:i + k], i))
        self.index.sort()

    def query(self, p):
        kmer = p[:self.k]
        i = bisect.bisect_left(self.index, (kmer, -1))
        hits = []
        while i < len(self.index):
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
print("Question 4")
print(len(approximate_match(p, chr1, 2)[0]))


# Question 5
# Using the instructions given in Question 4, how many total index hits are there when searching for occurrences of
# GGCGCGGTGGCTCACGCCTGTAAT\verb|GGCGCGGTGGCTCACGCCTGTAAT|GGCGCGGTGGCTCACGCCTGTAAT with up to 2 substitutions in the
# excerpt of human chromosome 1?
# (Don't consider reverse complements.)
# Hint: You should be able to use the boyer_moore\verb|boyer_moore|boyer_moore function
# (or the slower naive\verb|naive|naive function) to double-check your answer.

import bisect


p = "GGCGCGGTGGCTCACGCCTGTAAT"
chr1 = readGenome('chr1.GRCh38.excerpt.fasta')
print("Question 5")
print(approximate_match(p, chr1, 2)[1])


# Question 6
# Write a function that, given a length-24 pattern P and given a SubseqIndex\verb|SubseqIndex|SubseqIndex object built
# with k = 8 and ival = 3, finds all approximate occurrences of P within T with up to 2 mismatches.
# When using this function, how many total index hits are there when searching for
# GGCGCGGTGGCTCACGCCTGTAAT\verb|GGCGCGGTGGCTCACGCCTGTAAT|GGCGCGGTGGCTCACGCCTGTAAT with up to 2 substitutions in the
# excerpt of human chromosome 1? (Again, don't consider reverse complements.)
import bisect


class SubseqIndex(object):

    def __init__(self, t, k, ival):
        self.k = k
        self.ival = ival
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):
            self.index.append((t[i:i + self.span:ival], i))
        self.index.sort()

    def query(self, p):
        subseq = p[:self.span:self.ival]
        i = bisect.bisect_left(self.index, (subseq, -1))
        hits = []
        while i < len(self.index):
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

def approximate_match_subseq(p, t, n, ival):
    segment_length = int(round(len(p) / (n+1)))
    all_matches = set()
    p_idx = SubseqIndex(t, segment_length, ival)
    idx_hits = 0
    for i in range(n+1):
        start = i
        matches = p_idx.query(p[start:])

        for m in matches:
            idx_hits += 1
            if m < start or m-start+len(p) > len(t):
                continue

            mismatches = 0

            for j in range(0, len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break

            if mismatches <= n:
                all_matches.add(m - start)
    return list(all_matches), idx_hits


p = "GGCGCGGTGGCTCACGCCTGTAAT"
chr1 = readGenome('chr1.GRCh38.excerpt.fasta')
print("Question 6")
print(approximate_match_subseq(p, chr1, 2, 3)[1])
