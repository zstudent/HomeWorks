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


sequences, qualities = readFastq('ERR037900_1.first1000.fastq')
print(lowest_quality_base(qualities))