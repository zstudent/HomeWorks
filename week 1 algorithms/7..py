# Question 7
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


sequences, qualities = readFastq('ERR037900_1.first1000.fastq')
print(lowest_quality_base(qualities))