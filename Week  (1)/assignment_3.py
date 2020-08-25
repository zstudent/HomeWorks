#3
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

a = readGenome("lambda_virus.fa").find("ACTAAGT")
b = readGenome("lambda_virus.fa").find(reverseComplement("ACTAAGT"))

if a >= b:
    print(b)
else:
    print(a)

