# Question 1
# How many alignments does the naive exact matching algorithm try when matching the string
# GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG\verb|GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG|
# GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to the excerpt of human chromosome 1?
# (Don't consider reverse complements.)

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def naive(p, t):
    occurrences = []
    comparisons = 0
    alignments = 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        alignments += 1
        match = True
        for j in range(len(p)):  # loop over characters
            comparisons += 1
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences, comparisons, alignments


p = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
chr1 = readGenome("chr1.GRCh38.excerpt.fasta")

print(naive(p, chr1)[2])