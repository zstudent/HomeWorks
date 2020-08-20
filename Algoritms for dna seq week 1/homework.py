# Question 1
# How many times does \verb|AGGT|AGGT or its reverse complement (\verb|ACCT|ACCT) occur in the lambda
# virus genome? E.g. if \verb|AGGT|AGGT occurs 10 times and \verb|ACCT|ACCT occurs 12 times, you should report 22.


def readGenome(file):
    genome = ''
    with open(file, 'r') as f:
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

print("Question 1 ")
print(readGenome("lambda_virus.fa").count("AGGT") + readGenome("lambda_virus.fa").count(reverseComplement("AGGT")))


# Question 2
# How many times does \verb|TTAA|TTAA or its reverse complement occur in the lambda virus genome?
#
# Hint: \verb|TTAA|TTAA and its reverse complement are equal, so remember not to double count.

print("Question 2 ")
print(readGenome("lambda_virus.fa").count("TTAA"))
# Question 3
# What is the offset of the leftmost occurrence of \verb|ACTAAGT|ACTAAGT or
# its reverse complement in the Lambda virus genome? E.g. if the leftmost occurrence of
# \verb|ACTAAGT|ACTAAGT is at offset 40 (0-based) and the leftmost occurrence of its reverse complement
# \verb|ACTTAGT|ACTTAGT is at offset 29, then report 29.
print("Question 3")

x = readGenome("lambda_virus.fa").find("ACTAAGT")
y = readGenome("lambda_virus.fa").find(reverseComplement("ACTAAGT"))

if x >= y:
    print(y)
else:
    print(x)
# Question 4
# What is the offset of the leftmost occurrence of \verb|AGTCGA|AGTCGA or
# its reverse complement in the Lambda virus genome?



x = readGenome("lambda_virus.fa").find("AGTCGA")
y = readGenome("lambda_virus.fa").find(reverseComplement("AGTCGA"))
print("Question 4")
print(x)
print(y)

# Question 5
# As we will discuss, sometimes we would like to find approximate matches for P in T.
# That is, we want to find occurrences with one or more differences.

# For Questions 5 and 6, make a new version of the \verb|naive|naive function called \verb|naive_2mm|naive_2mm that allows up to 2 mismatches per occurrence.
# Unlike for the previous questions, do not consider the reverse complement here. We're looking for approximate matches for P itself, not its reverse complement.
# For example, \verb|ACTTTA|ACTTTA occurs twice in \verb|ACTTACTTGATAAAGT|ACTTACTTGATAAAGT, once at offset 0 with 2 mismatches, and once at offset 4 with 1 mismatch.
# So \verb|naive_2mm('ACTTTA', 'ACTTACTTGATAAAGT')|naive_2mm(’ACTTTA’, ’ACTTACTTGATAAAGT’) should return the list \verb|[0, 4]|[0, 4].

# Hint: See this notebook for a few examples you can use to test your \verb|naive_2mm|naive_2mm function.
# How many times does \verb|TTCAAGCC|TTCAAGCC occur in the Lambda virus genome when allowing up to 2 mismatches?


def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        mis_score = 0
        for j in range(len(p)):
            if t[i+j] != p[j]:
                mis_score += 1
                if mis_score > 2:
                    break
        if mis_score <= 2:
            occurrences.append(i)
    return occurrences

p = "TTCAAGCC"
print("Question 5")
print(naive_2mm(p, readGenome("lambda_virus.fa")))
# Question 6
# What is the offset of the leftmost occurrence of
# \verb|AGGAGGTT|AGGAGGTT in the Lambda virus genome when allowing up to 2 mismatches?






def naive_2mm_first_occ(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        mis_score = 0
        for j in range(len(p)):
            if t[i + j] != p[j]:
                mis_score += 1
                if mis_score > 2:
                    break
        if mis_score <= 2:
            occurrences.append(i)
            break
    return occurrences


p = "AGGAGGTT"
print("Question 6")
print(naive_2mm_first_occ(p, readGenome("lambda_virus.fa")))
# Question 7
# Finally, download and parse the provided FASTQ file containing real DNA sequencing reads derived from a human:
#
# https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ERR037900_1.first1000.fastq
#
# Note that the file has many reads in it and you should examine all of them together when answering this question. The reads are taken from this study:
#
# Ajay, S. S., Parker, S. C., Abaan, H. O., Fajardo, K. V. F., & Margulies, E. H. (2011). Accurate
#
# and comprehensive sequencing of personal genomes. Genome research, 21(9), 1498-1505.
#
# This dataset has something wrong with it; one of the sequencing cycles is poor quality.
#
# Report which sequencing cycle has the problem. Remember that a sequencing cycle corresponds to a particular offset in all the reads. For example, if the leftmost read position seems to have a problem consistently across reads, report 0. If the fourth position from the left has the problem, report 3. Do whatever analysis you think is needed to identify the bad cycle. It might help to review the "Analyzing reads by position" video.

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

print("Question 7")
sequences, qualities = readFastq('ERR037900_1.first1000.fastq')
print(lowest_quality_base(qualities))