# Question 3
# What is the offset of the leftmost occurrence of \verb|ACTAAGT|ACTAAGT or
# its reverse complement in the Lambda virus genome? E.g. if the leftmost occurrence of
# \verb|ACTAAGT|ACTAAGT is at offset 40 (0-based) and the leftmost occurrence of its reverse complement
# \verb|ACTTAGT|ACTTAGT is at offset 29, then report 29.

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

a = readGenome("lambda_virus.fa").find("ACTAAGT")
b = readGenome("lambda_virus.fa").find(reverseComplement("ACTAAGT"))

if a >= b:
    print(b)
else:
    print(a)

