



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
        count = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                count += 1
                if count >= 3:
                    break
        if count <= 2:
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


virus = readGenome('lambda_virus.fa')
sequences, qualities = readFastq('ERR037900_1.first1000.fastq')


matches = naive('AGGT', virus)
matches.extend(naive(reverseComplement('ACCT'), virus))
print(len(matches))


matches = naive('TTAA', virus)
print(len(matches))


matches = naive('ACTAAGT', virus)
matches.extend(naive(reverseComplement('ACTAAGT'), virus))
n = matches[0]
print(n)


matches = naive('AGTCGA', virus)
matches.extend(naive(reverseComplement('AGTCGA'), virus))
print(matches[0])


matches = naive_2mm('TTCAAGCC', virus)
print(len(matches))

matches = naive_2mm('TTCAAGCC', virus)
print(matches[0])


totals = [0] * len(qualities[0])
for q in qualities:
    for i,n in enumerate(q):
        totals[i] += phred33ToQ(n)

print(totals.index(min(totals)))


