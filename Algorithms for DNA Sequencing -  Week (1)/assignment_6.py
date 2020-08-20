# Question 6
# What is the offset of the leftmost occurrence of
# \verb|AGGAGGTT|AGGAGGTT in the Lambda virus genome when allowing up to 2 mismatches?

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


# print(readGenome("lambda_virus.fa").count("TTAA"))

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t


def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        mismatches = 0
        for j in range(len(p)):  # loop over characters
            if t[i + j] != p[j]:  # compare characters
                mismatches += 1
                if mismatches > 2:
                    break
        if mismatches <= 2:
            occurrences.append(i)
            break                   # all chars matched; record
    return occurrences


p = "AGGAGGTT"

print(naive_2mm(p, readGenome("lambda_virus.fa")))
