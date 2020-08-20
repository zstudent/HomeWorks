# Question 2
# How many times does \verb|TTAA|TTAA or its reverse complement occur in the lambda virus genome?
#
# Hint: \verb|TTAA|TTAA and its reverse complement are equal, so remember not to double count.



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

print(readGenome("lambda_virus.fa").count("TTAA"))

