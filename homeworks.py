def dna(filename):
    genome = ''
    with open(filename, 'r') as myfile:
        for line in myfile:

            if not line[0] == '>':
                genome += line.rstrip()
    return genome


def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t


def naive_algorithm (p, t):
    p_rev = reverseComplement(p)
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        match = True
        for j in range(len(p)):
            if t[i+j] != p[j]:
                match = False
                break
        if not match:
            match = True
            for j in range(len(p)):
                if t[i + j] != p_rev[j]:
                    match = False
                    break
        if match:
            occurrences.append(i)
    return occurrences
viraldna = 'lambda_virus.fa'
genome = dna(viraldna)
print(len(naive_algorithm('AGGT', genome)))
print(len(naive_algorithm('TTAA', genome)))
print(min(naive_algorithm('ACTAAGT', genome)))
print(min(naive_algorithm('AGTCGA', genome)))
def naiv_algorithm2(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        count_mismatch = 0
        for j in range(len(p)):
            if t[i+j] != p[j]:
                count_mismatch += 1
        if count_mismatch <= 2:
            occurrences.append(i)
    return occurrences
print(len(naiv_algorithm2('TTCAAGCC',genome )))
print(min(naive_algorithm('AGGAGGTT', genome)))
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


def lowest_quality_base(qs):
    total = [0] * len(qs[0])
    for q in qs:
        for i, phred in enumerate(q):
            total[i] += phred33ToQ(phred)
    return total.index(min(total))


def phred33ToQ(qual):
    return ord(qual) - 33


sequences, qualities = readFastq('ERR037900_1.first1000.fastq')
print(lowest_quality_base(qualities))
