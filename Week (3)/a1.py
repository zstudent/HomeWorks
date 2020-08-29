# Question 1

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:

            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def approximate_match(p, t):

    D = []
    for i in range(len(p)+1):
        D.append([0]*(len(t)+1))


    for i in range(len(p)+1):
        D[i][0] = i
    for i in range(len(t)+1):
        D[0][i] = 0


    for i in range(1, len(p)+1):
        for j in range(1, len(t)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if p[i-1] == t[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)


    return min(D[-1])

chr1 = readGenome('chr1.GRCh38.excerpt.fasta')
p = "GCTGATCGATCGTACG"
print(approximate_match(p, chr1))
