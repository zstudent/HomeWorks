# Question 5
# As we will discuss, sometimes we would like to find approximate matches for P in T.
# That is, we want to find occurrences with one or more differences.

# For Questions 5 and 6, make a new version of the \verb|naive|naive function called \verb|naive_2mm|naive_2mm that allows up to 2 mismatches per occurrence.
# Unlike for the previous questions, do not consider the reverse complement here. We're looking for approximate matches for P itself, not its reverse complement.
# For example, \verb|ACTTTA|ACTTTA occurs twice in \verb|ACTTACTTGATAAAGT|ACTTACTTGATAAAGT, once at offset 0 with 2 mismatches, and once at offset 4 with 1 mismatch.
# So \verb|naive_2mm('ACTTTA', 'ACTTACTTGATAAAGT')|naive_2mm(’ACTTTA’, ’ACTTACTTGATAAAGT’) should return the list \verb|[0, 4]|[0, 4].

# Hint: See this notebook for a few examples you can use to test your \verb|naive_2mm|naive_2mm function.
# How many times does \verb|TTCAAGCC|TTCAAGCC occur in the Lambda virus genome when allowing up to 2 mismatches?

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


def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        mismatches = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                mismatches += 1
                if mismatches > 2:
                    break
        if mismatches <= 2:
            occurrences.append(i)  # all chars matched; record
    return occurrences

p = "TTCAAGCC"

print(naive_2mm(p, readGenome("lambda_virus.fa")))