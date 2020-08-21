#!/usr/bin/env python
# coding: utf-8

# In[9]:


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


# In[5]:


def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t


# In[1]:


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


# In[21]:


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


# In[14]:


def leftMost(s,genome):
    rev_comp_s = reverseComplement(s)
    firstIndex = naive(s,genome)[0]
    scndIndex = naive(rev_comp_s,genome)[0] 
    if firstIndex < scndIndex:
        return (firstIndex)
    elif firstIndex > scndIndex:
        return (scndIndex)
    else:
        return (firstIndex)


# In[22]:


def naive_2mm(p,t):
    arr = []
    for i in range(len(t) - len(p) + 1):
        match = True
        mismatch_count = 0
        for j in range(len(p)):
            if t[i+j] != p[j]:
                mismatch_count += 1
            if mismatch_count > 2:
                match = False
                break
        if match: 
            arr.append(i)
    return arr


# In[33]:


def pos(qualities):
    hist=[0]*len(qualities[0])
    for i in qualities:
        for j,k in enumerate(i):
            hist[j]+=phred33toQ(k)
    return hist.index(min(hist))


# In[23]:


def phred33toQ(phred):
    return ord(phred)-33


# In[3]:


genome=readGenome('lambda_virus.fa')


# # 1)

# In[6]:


print(genome.count('AGGT')+genome.count(reverseComplement('AGGT')))


# # 2)

# In[7]:


print(genome.count('TTAA'))


# # 3)

# In[15]:


print(leftMost('ACTAAGT',genome))


# # 4)

# In[16]:


print(leftMost('AGTCGA',genome))


# # 5)

# In[18]:


print(naive_2mm('ACT','ACTACTT'))


# # 6)

# In[20]:


print(naive_2mm('AGGAGGTT',genome)[0])


# # 7)

# In[34]:


sequences, qualities = readFastq('ERR037900_1.first1000.fastq')
print(pos(qualities))

