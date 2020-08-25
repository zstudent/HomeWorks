# Question 7

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()
            sequence = fh.readline().rstrip()
            fh.readline()
            quality = fh.readline().rstrip()
            if len(sequence) == 0:
                break
            sequences.append(sequence)
            qualities.append(quality)
    return sequences, qualities


def lowest_quality_base(qual):
    total = [0] * len(qual[0])
    for qu in qual:
        for i, phred in enumerate(qu):
            total[i] += phred33ToQ(phred)
    return total.index(min(total))


def phred33ToQ(qual):
    return ord(qual) - 33


sequences, qualities = readFastq('ERR037900_1.first1000.fastq')
print(lowest_quality_base(qualities))
