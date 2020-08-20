# Question 4
# What is the offset of the leftmost occurrence of \verb|AGTCGA|AGTCGA or
# its reverse complement in the Lambda virus genome?

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

a = readGenome("lambda_virus.fa").find("AGTCGA")
b = readGenome("lambda_virus.fa").find(reverseComplement("AGTCGA"))

print(a)
print(b)


