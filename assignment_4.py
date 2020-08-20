# Question 4

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

seq = readGenome("lambda_virus.fa").find("AGTCGA")
revSeq = readGenome("lambda_virus.fa").find(reverseComplement("AGTCGA"))

print(seq)
print(revSeq)


