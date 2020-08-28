## Question 1

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

print("Question #1")
print(readGenome("lambda_virus.fa").count("AGGT") + readGenome("lambda_virus.fa").count(reverseComplement("AGGT")))

## Question 2

print("Question #2")
print(readGenome("lambda_virus.fa").count("TTAA"))

## Question 3

EDS = readGenome("lambda_virus.fa").find("ACTAAGT")
Glucose = readGenome("lambda_virus.fa").find(reverseComplement("ACTAAGT"))

print("Question #3")
if EDS >= Glucose:
    print(Glucose)
else:
    print(EDS)

## Question 4

ROE = readGenome("lambda_virus.fa").find("AGTCGA")
SOE = readGenome("lambda_virus.fa").find(reverseComplement("AGTCGA"))
print("Question #4")
print(ROE)
print(SOE)

# Question 5

def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        mismatches = 0
        for CS in range(len(p)):
            if t[i+CS] != p[CS]:
                mismatches += 1
                if mismatches > 2:
                    break
        if mismatches <= 2:
            occurrences.append(i)
    return occurrences

p = "TTCAAGCC"

print("Question #5")
print(naive_2mm(p, readGenome("lambda_virus.fa")))

##Question 6


def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        mismatches = 0
        for j in range(len(p)):
            if t[i + j] != p[j]:
                mismatches += 1
                if mismatches > 2:
                    break
        if mismatches <= 2:
            occurrences.append(i)
    return occurrences


p = "AGGAGGTT"
print("Question #6")
print(naive_2mm(p, readGenome("lambda_virus.fa"))[0])


# Question 7

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
print("Question #7")
print(lowest_quality_base(qualities))



