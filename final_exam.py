import operator
import time

def count_records(contents):
    return contents.count(">")


def create_dict(contents):
    sequences_dict = {}
    for line in contents:
        if line.startswith('>'):
            header = line
            sequences_dict[header] = ""
        else:
            sequences_dict[header] += line
    return sequences_dict

def longest_shortest(sequences_dict):
    lengths = [len(i) for i in sequences_dict.values()]

    minimum, maximum = max(lengths), min(lengths)

    minimum_list, maximum_list = [], []

    for key, value in sequences_dict.items():
        if len(value) == minimum:
            minimum_list.append(key)
        if len(value) == maximum:
            maximum_list.append(key)

def find_orfs(sequences_dict, frame, reversed):
    start_codone = 'ATG'
    stop_codones = ['TAA', 'TAG', 'TGA']
    orf_dict = {}
    for header, sequence in sequences_dict.items():
        if reversed:
            sequence = sequence[: : -1]
        orf = ''
        found = False
        orf_list = []
        for i in range(frame - 1, len(sequence), 3):
            triplet = sequence[i : i + 3]
            if triplet in stop_codones:
                if len(orf) > 0:
                    orf += triplet
                    orf_list.append(orf)
                    orf = ''
                    found = False
                else:
                    found = False
                    orf = ''

            if found:
                orf += triplet

            elif start_codone == triplet:
                orf += start_codone
                found = True

        orf_dict[header] = orf_list
    return orf_dict

def longest_orf(orf_dict):
    maximum_orfs = {}
    maximum_orf = ''
    absolute_longest = 0
    for header, orfs in orf_dict.items():
        longest = 0
        for orf in orfs:
            length = len(orf)
            if length > longest:
                longest = length
                if longest > absolute_longest:
                    absolute_longest = longest
                    maximum_orf = orf

        maximum_orfs[header] = longest
    return sorted(maximum_orfs.items(), key=operator.itemgetter(1))[-1], maximum_orf, maximum_orfs

def find_index(sequences_dict, maximum_orf, longest_orf_tuple):
    return sequences_dict[longest_orf_tuple[0]].find(maximum_orf)

def find_repeats(sequences_dict, repeat_length):
    repeats_dict = {}
    for sequence in sequences_dict.values():
        for i in range(len(sequence) - repeat_length + 1):
            repeat = sequence[i : i + repeat_length]
            if repeat in repeats_dict.keys():
                repeats_dict[repeat] += 1
            else:
                repeats_dict[repeat] = 1
    # most_repeat = max(repeats_dict.values())
    most_repeat = sorted(repeats_dict, key = operator.itemgetter(1))[-1]
    return repeats_dict, most_repeat


if __name__ == "__main__":
    file_name = "dna.example.fasta"
    with open(file_name, "r") as file:
        contents = file.read().split("\n")
        sequences_dict = create_dict(contents)
        longest_shortest(sequences_dict)
        orf_dict = find_orfs(sequences_dict, 1, False)
        longest_orf_tuple, maximum_orf, maximum_orfs = longest_orf(orf_dict)
        print(maximum_orf)
        print(longest_orf_tuple)
        print(find_index(sequences_dict, maximum_orf, longest_orf_tuple))
        print(find_repeats(sequences_dict, 10))