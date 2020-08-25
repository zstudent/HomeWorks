#1
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

print(naive(p, chr1)[2])