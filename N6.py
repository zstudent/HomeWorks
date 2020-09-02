import bisect

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:

            if not line[0] == '>':
                genome += line.rstrip()
    return genome

class SubseqIndex(object):

    def __init__(self, t, k, ival):

        self.k = k
        self.ival = ival
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):
            self.index.append((t[i:i + self.span:ival], i))
        self.index.sort()

    def query(self, p):
        subseq = p[:self.span:self.ival]
        i = bisect.bisect_left(self.index, (subseq, -1))
        hits = []
        while i < len(self.index):
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

def approximate_match_subseq(p, t, n, ival):
    segment_length = int(round(len(p) / (n+1)))
    all_matches = set()
    p_idx = SubseqIndex(t, segment_length, ival)
    idx_hits = 0
    for i in range(n+1):
        start = i
        matches = p_idx.query(p[start:])


        for m in matches:
            idx_hits += 1
            if m < start or m-start+len(p) > len(t):
                continue

            mismatches = 0

            for j in range(0, len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break

            if mismatches <= n:
                all_matches.add(m - start)
    return list(all_matches), idx_hits


p = "GGCGCGGTGGCTCACGCCTGTAAT"
chr1 = readGenome('chr1.GRCh38.excerpt.fasta')

print(approximate_match_subseq(p, chr1, 2, 3)[1])