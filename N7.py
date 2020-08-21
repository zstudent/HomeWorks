import matplotlib.pyplot as plt

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def phred33ToQ(qual):
    return ord(qual) - 33

def createHist (qualities):
    hist=[0]*len(qualities[0])
    for qual in qualities:
        for phred in qual:
            q=phred33ToQ(phred)
            hist[q]+=1
    return hist

sequences, qualities = readFastq('ERR037900_1.first1000.fastq')
h=createHist(qualities)
plt.plot(range (len(h)),h)
plt.show()
# 2 is the numeric value of #. the plot shows that rhis is dubios




