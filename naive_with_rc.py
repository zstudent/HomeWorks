def read_genome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


def reverse_complement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for letter in s:
        t  = complement[letter] + t
    return t


def naive_with_rc(p, t):
    occurrences = []
    p_reverse_complement = reverse_complement(p)
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match, match_reverse_complement = True, True
        for j in range(len(p)):  # loop over characters
            if t[i + j] != p[j]:  # compare characters
                match = False
                break

        if not match:
            for j in range(len(p)):
                if t[i + j] != p_reverse_complement[j]:  # compare characters
                    match_reverse_complement = False
                    break

        if match or match_reverse_complement:
            occurrences.append(i)  # all chars matched; record
    return occurrences


if __name__ == '__main__':
    genome = read_genome('lambda_virus.fa')
    pattern = 'ACT'
    occurences = naive_with_rc(pattern, genome)
    print(occurences)
    print(genome[195 : 198])



