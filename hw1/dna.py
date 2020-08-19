import time

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

# p-pattern, t-text
def naive_with_rc(p, t):
    rc = reverseComplement(p)
    if rc == p:
        return naive(p, t)
    return naive(p, t)+naive(rc, t)

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


def main():
    phix_genome = readGenome('phix.fa')
    occurrences = naive_with_rc('ATTA', phix_genome)
    print('offset of leftmost occurrence: %d' % min(occurrences)) #should be 62
    print('# occurrences: %d' % len(occurrences)) #should be 60

if __name__ == "__main__":
    start = int(round(time.time() * 1000))
    main()
    finish = int(round(time.time() * 1000))
    overall = finish - start
    print("took: "+ str(overall) + " milliseconds")
