#!/usr/bin/env python
# coding: utf-8

# In[1]:


def editDistance(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = i
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    return D[-1][-1]


# In[3]:


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


# In[13]:


def overlap_pairs(reads, k, map={}):
    def get_kmers(read, k):
        res = set()
        for i in range(0, len(read)-k+1):
            res.add(read[i:i+k])
        return res
    for read in reads:
        kmers = get_kmers(read, k)
        for kmer in kmers:
            if not kmer in map.keys():
                map[kmer] = set()
            map[kmer].add(read)
    pairs = []
    for head in reads:
        kmer = head[-k:]
        candidates = map[kmer]
        for tail in candidates:
            if (not head == tail and overlap(head, tail, k)):
                pairs.append((head, tail))
    return pairs


# In[12]:


def readFAST_Q(filename):
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
reads, _ = readFAST_Q('ERR266411_1.for_asm.fastq')


# # 1)

# In[8]:


chromosome = readGenome('chr1.GRCh38.excerpt.fasta')
print(editDistance('GCTGATCGATCGTACG', chromosome))


# # 2)

# In[9]:


print(editDistance('GATTTACCAGATTGAG', chromosome))


# # 3)

# In[15]:


pairs = overlap_pairs(reads, 30)
print(len(pairs))


# # 4)

# In[16]:


print(len(set(pair[0] for pair in pairs)))

