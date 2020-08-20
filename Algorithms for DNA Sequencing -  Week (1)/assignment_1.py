# Question 1
# How many times does \verb|AGGT|AGGT or its reverse complement (\verb|ACCT|ACCT) occur in the lambda
# virus genome? E.g. if \verb|AGGT|AGGT occurs 10 times and \verb|ACCT|ACCT occurs 12 times, you should report 22.


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


print(readGenome("lambda_virus.fa").count("AGGT") + readGenome("lambda_virus.fa").count(reverseComplement("AGGT")))
