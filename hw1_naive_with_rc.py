
def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[ i +j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        mismatches = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                mismatches += 1
                if mismatches > 2:
                    break
        if mismatches <= 2:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def phred33ToQ(qual):
    return ord(qual) - 33


genome = readGenome('lambda_virus.fa')
sequences, qualities = readFastq('ERR037900_1.first1000.fastq')
############## Problem #1 ####################
p = 'AGGT'
rcp = reverseComplement(p)
matches = naive(p, genome)
matches.extend(naive(rcp, genome))
print("%s or %s appears in genome %d times" % (p, rcp, len(matches)))

############## Problem #2 ####################
p = 'TTAA'
matches = naive(p, genome)
print("%s appears in genome %d times" % (p, len(matches)))


############## Problem #3 ####################
p = 'ACTAAGT'
rcp = reverseComplement(p)
matches = naive(p, genome)
matches.extend(naive(rcp, genome))
leftmost = sorted(matches)[0]
print("Leftmost appearance of %s or %s is %d" % (p, rcp, leftmost))

############## Problem #4 ####################
p = 'AGTCGA'
rcp = reverseComplement(p)
matches = naive(p, genome)
matches.extend(naive(rcp, genome))
leftmost = sorted(matches)[0]
print("Leftmost appearance of %s or %s is %d" % (p, rcp, leftmost))

############## Problem #5 ####################
p = 'TTCAAGCC'
matches = naive_2mm(p, genome)
print("%s appears in genome %d times" % (p, len(matches)))

############## Problem #6 ####################
p = 'TTCAAGCC'
matches = naive_2mm(p, genome)
leftmost = sorted(matches)[0]
print("Leftmost appearance of %s is %d" % (p, leftmost))

############## Problem #7 ####################
def findPCByPosition(qs):
    totals = [0] * len(qs[0])
    for q in qs:
        for i in range(len(q)):
            totals[i] += phred33ToQ(q[i])
    return totals.index(min(totals))

cycleIndex = findPCByPosition(qualities)
print("Sequencing cycle %d has the problem" % cycleIndex)


